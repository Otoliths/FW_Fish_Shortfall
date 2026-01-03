library(furrr)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(sf)
library(ggplot2)
library(cli)

# Suppress all warnings
options(warn = -1)
sf::sf_use_s2(F) 
# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
path_species_dir   <- "input/processed/basin_sp/"
path_occ_dir       <- "input/processed/occ/"
path_output_dir    <- "input/data_prep/basin_preserved_specimen/"

# Ensure output directory exists
if (!dir.exists(path_output_dir)) {
  dir.create(path_output_dir, recursive = TRUE)
  cli_alert_success("Created directory: {path_output_dir}")
}

# ------------------------------------------------------------
# Load and prepare basin polygons
# ------------------------------------------------------------
inland <- readRDS("input/raw/basin/basin_sf_v1.rds") %>%
  st_transform(4326) %>%
  suppressWarnings(st_make_valid()) %>%
  # Wrap around the dateline to avoid invalid geometries across 180°
  st_wrap_dateline(
    options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
    quiet   = TRUE
  ) %>%
  st_collection_extract("POLYGON", warn = FALSE) %>%
  st_as_sf()

# Extract basin IDs once
basin_ids <- inland %>%
  st_drop_geometry() %>%
  pull(basin_id)


# Create a global progress bar
cli_progress_bar("Processing basins", total = length(basin_ids))

# ------------------------------------------------------------
# Helper: does a basisOfRecord entry represent a preserved specimen?
# We accept multiple variants (case differences and formatting differences)
# ------------------------------------------------------------
is_preserved_specimen <- function(x) {
  # normalize to lowercase and remove spaces/underscores
  norm <- x %>%
    tolower() %>%
    gsub("_", "", ., fixed = TRUE) %>%
    gsub(" ", "", ., fixed = TRUE)
  
  norm %in% c(
    "PRESERVED_SPECIMEN",   # standard Darwin Core
    "preservedspecimen",   # underscore form
    "PreservedSpecimen"    # (kept multiple on purpose, harmless)
  )
}

# ------------------------------------------------------------
# Core worker: process a single basin
# ------------------------------------------------------------
process_basin <- function(basin_id_value) {
  
  cli_h2(paste("Processing basin:", basin_id_value))
  
  # ------------------------------------------------------------
  # 1. Read species table for this basin
  # ------------------------------------------------------------
  spp_file <- paste0(path_species_dir, basin_id_value, ".csv")
  
  # Skip basin if no species file
  if (!file.exists(spp_file)) {
    cli_alert_warning(paste("Species file not found:", spp_file, "- skipping basin"))
    return(NULL)
  }
  
  # Keep all columns from spp so we can "add" new fields on top
  spp <- read.csv(spp_file)[,c(1,3)]
  
  # Basic validation
  if (!is.data.frame(spp)) {
    cli_alert_warning("Species table is not a data frame - skipping this basin")
    return(NULL)
  }
  if (!"valid_name" %in% colnames(spp)) {
    cli_alert_warning("'valid_name' column missing - skipping this basin")
    return(NULL)
  }
  
  # ------------------------------------------------------------
  # 2. Container for results:
  #    start with an empty data frame with same columns as spp
  # ------------------------------------------------------------
  all_results <- spp[0, , drop = FALSE]
  
  # ------------------------------------------------------------
  # 3. Loop over all species in this basin
  # ------------------------------------------------------------
  for (species in unique(spp$valid_name)) {
    
    # Subset the basin geometry for this basin_id
    basin <- inland %>% dplyr::filter(basin_id %in% basin_id_value)
    
    # Path to species occurrence file (.rds of point records)
    occ_file <- paste0(
      path_occ_dir,
      gsub(" ", "_", species),
      ".rds"
    )
    
    # If species occurrence file is missing: keep spp rows, add NA fields
    if (!file.exists(occ_file)) {
      cli_alert_warning(
        paste("Occurrence file not found:", occ_file, "- setting year / preserved_specimen = NA")
      )
      
      species_rows <- spp %>%
        dplyr::filter(valid_name == species) %>%
        dplyr::mutate(
          year               = NA_integer_,
          preserved_specimen = NA_integer_
        )
      
      all_results <- dplyr::bind_rows(all_results, species_rows)
      next
    }
    
    cli_alert_info(paste("Reading occurrence data for species:", species))
    
    # Load occurrences and convert to sf in WGS84
    occ_sf <- readRDS(occ_file) %>%
      sf::st_as_sf(
        coords = c("decimalLongitude", "decimalLatitude"),
        crs    = 4326
      )
    
    # Intersect species records with the basin polygon
    occ_in_basin <- sf::st_intersection(occ_sf, basin)
    
    # If no points fall inside the basin: keep spp rows, add NA fields
    if (nrow(occ_in_basin) == 0) {
      cli_alert_warning(
        paste("No occurrences inside basin for species:", species,
              "- setting year / preserved_specimen = NA")
      )
      
      species_rows <- spp %>%
        dplyr::filter(valid_name == species) %>%
        dplyr::mutate(
          year               = NA_integer_,
          preserved_specimen = NA_integer_
        )
      
      all_results <- dplyr::bind_rows(all_results, species_rows)
      next
    }
    
    # ------------------------------------------------------------
    # 3.1 Compute cumulative count of preserved specimens per year
    # ------------------------------------------------------------
    grouped_data <- occ_in_basin %>%
      sf::st_drop_geometry() %>%
      # keep only preserved specimen records
      dplyr::filter(is_preserved_specimen(basisOfRecord)) %>%
      # ensure we have a year column
      dplyr::filter(!is.na(year)) %>%
      dplyr::select(species, year, basisOfRecord) %>%
      dplyr::group_by(species, year) %>%
      dplyr::summarise(
        preserved_specimen = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::arrange(year) %>%
      dplyr::group_by(species) %>%
      dplyr::mutate(
        preserved_specimen = cumsum(preserved_specimen),
        species            = species
      ) %>%
      dplyr::ungroup()
    
    # If after filtering there are no preserved specimen records:
    # keep spp rows, add NA
    if (nrow(grouped_data) == 0) {
      cli_alert_warning(
        paste("No preserved specimens for species:", species,
              "- setting year / preserved_specimen = NA")
      )
      
      species_rows <- spp %>%
        dplyr::filter(valid_name == species) %>%
        dplyr::mutate(
          year               = NA_integer_,
          preserved_specimen = NA_integer_
        )
      
      all_results <- dplyr::bind_rows(all_results, species_rows)
      next
    }
    
    # ------------------------------------------------------------
    # 3.2 Join species-year data back into spp rows for this species
    #     - grouped_data has: species, year, preserved_specimen
    #     - spp has: valid_name + any other fields
    #     - join key: valid_name = species
    # ------------------------------------------------------------
    species_rows <- spp %>%
      dplyr::filter(valid_name == species) %>%
      dplyr::left_join(
        grouped_data,
        by = c("valid_name" = "species")
      )
    
    # Append to basin-wide result table
    all_results <- dplyr::bind_rows(all_results, species_rows)
    
    cli_alert_info(
      paste("Finished species", species, "in basin", basin_id_value)
    )
  }
  
  # ------------------------------------------------------------
  # 4. Write basin-level output: spp + year + preserved_specimen
  # ------------------------------------------------------------
  out_file <- paste0(path_output_dir, basin_id_value, ".csv")
  write.csv(all_results, file = out_file, row.names = FALSE)
  
  cli_alert_success(
    paste("All species processed for basin", basin_id_value, "→ saved to", out_file)
  )
  
  # Cleanup after basin
  gc()
}


# ------------------------------------------------------------
# Run for all basins
# ------------------------------------------------------------
purrr::walk(basin_ids, process_basin)

cli_progress_done()
gc()
