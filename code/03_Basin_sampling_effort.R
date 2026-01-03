library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(sf)
library(cli)

# Suppress all warnings
options(warn = -1)

# ------------------------------------------------------------
# File paths
# ------------------------------------------------------------
path_species_dir  <- "input/processed/basin_sp/"
path_occ_dir      <- "input/processed/occ/"
path_output_dir   <- "input/data_prep/basin_sampling_effort/"

# Ensure output directory exists
if (!dir.exists(path_output_dir)) {
  dir.create(path_output_dir, recursive = TRUE)
  cli_alert_success("Created directory: {path_output_dir}")
}

# ------------------------------------------------------------
# Load and fix basin geometry
# ------------------------------------------------------------
inland <- readRDS("input/raw/basin/basin_sf_v1.rds") %>%
  st_transform(4326) %>%
  suppressWarnings(st_make_valid()) %>%
  # Wrap around the dateline to fix global geometries
  st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE) %>%
  st_collection_extract("POLYGON", warn = FALSE) %>%
  st_as_sf()

# Get all basin IDs
basin_ids <- inland %>%
  st_drop_geometry() %>%
  pull(basin_id)


# Create progress bar
cli_progress_bar("Processing basins", total = length(basin_ids))

# ------------------------------------------------------------
# Function: process a single basin
# ------------------------------------------------------------
process_basin <- function(basin_id_value) {
  
  cli_h2(paste("Processing basin:", basin_id_value))
  
  # ------------------------------------------------------------
  # 1. Read species list for this basin
  # ------------------------------------------------------------
  spp_file <- paste0(path_species_dir, basin_id_value, ".csv")
  if (!file.exists(spp_file)) {
    cli_alert_warning(paste("File not found:", spp_file, "- skipping"))
    return(NULL)
  }
  
  # Keep all columns from spp so that we can "add" fields on top
  spp <- read.csv(spp_file)[, c(1, 3)]
  
  # Basic validation
  if (!is.data.frame(spp)) {
    cli_alert_warning("Invalid species table - skipping")
    return(NULL)
  }
  if (!"valid_name" %in% colnames(spp)) {
    cli_alert_warning("'valid_name' column missing - skipping")
    return(NULL)
  }
  
  # ------------------------------------------------------------
  # 2. Prepare a container for spp + sampling fields
  #    (same structure as spp, will append species-by-species)
  # ------------------------------------------------------------
  all_results <- spp[0, , drop = FALSE]  # empty with same columns as spp
  
  # ------------------------------------------------------------
  # 3. Loop through all unique species in this basin
  # ------------------------------------------------------------
  for (species in unique(spp$valid_name)) {
    
    # Basin geometry (same basin for all species)
    basin <- inland %>%
      dplyr::filter(basin_id %in% basin_id_value)
    
    # Path to species occurrence file
    distribution_file <- paste0(
      path_occ_dir,
      gsub(" ", "_", species),
      ".rds"
    )
    
    # If occurrence file does not exist: keep spp rows, add NA sampling fields
    if (!file.exists(distribution_file)) {
      cli_alert_warning(
        paste("Occurrence file not found:", distribution_file, "- setting sampling fields to NA")
      )
      
      species_rows <- spp %>%
        dplyr::filter(valid_name == species) %>%
        dplyr::mutate(
          sampling_year   = NA_integer_,
          sampling_effort = NA_real_
        )
      
      all_results <- dplyr::bind_rows(all_results, species_rows)
      next
    }
    
    cli_alert_info(
      paste("Reading occurrence data for species:", species)
    )
    
    # Read occurrence data and convert to sf points
    distribution_sf <- readRDS(distribution_file) %>%
      sf::st_as_sf(
        coords = c("decimalLongitude", "decimalLatitude"),
        crs    = 4326
      )
    
    # --------------------------------------------------------
    # 3.1 Reproject basin to equal-area projection (Eckert IV)
    # --------------------------------------------------------
    basin_proj <- sf::st_transform(
      basin,
      crs = "+proj=eck4 +lon_0=0 +lat_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    )
    
    # --------------------------------------------------------
    # 3.2 Create 50 km fishnet grid clipped to basin boundary
    # --------------------------------------------------------
    fishnet_clipped_sf <- sf::st_make_grid(
      basin_proj,
      cellsize = c(50000, 50000),
      what     = "polygons"
    ) %>%
      suppressWarnings(sf::st_make_valid()) %>%
      sf::st_intersection(basin_proj) %>%
      sf::st_as_sf() %>%
      dplyr::mutate(grid_id = dplyr::row_number())
    
    # --------------------------------------------------------
    # 3.3 Project species points to the same CRS as fishnet
    # --------------------------------------------------------
    points_sf_in_basin <- sf::st_transform(
      distribution_sf,
      crs = sf::st_crs(fishnet_clipped_sf)
    )
    
    # --------------------------------------------------------
    # 3.4 Intersect points with grid cells and compute
    #     yearly sampling effort:
    #     - count records per grid cell per year
    #     - take mean(count) across grid cells per year
    # --------------------------------------------------------
    grouped_data <- sf::st_intersection(
      fishnet_clipped_sf,
      points_sf_in_basin
    ) %>%
      sf::st_drop_geometry() %>%
      dplyr::filter(!is.na(year)) %>%
      dplyr::group_by(year, grid_id) %>%
      dplyr::summarise(
        count = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::select(year, count) %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(
        sampling_effort = mean(count, na.rm = TRUE),
        .groups = "drop"
      )
    
    # If for some reason no valid yearly records → set NA
    if (nrow(grouped_data) == 0) {
      cli_alert_warning(
        paste("No valid yearly sampling effort for species:", species, "- setting fields to NA")
      )
      
      species_rows <- spp %>%
        dplyr::filter(valid_name == species) %>%
        dplyr::mutate(
          sampling_year   = NA_integer_,
          sampling_effort = NA_real_
        )
      
      all_results <- dplyr::bind_rows(all_results, species_rows)
      next
    }
    
    # --------------------------------------------------------
    # 3.5 Attach basin_id, valid_name and rename year → sampling_year
    #     This is the "robust" sampling table:
    #     basin_id, valid_name, sampling_year, sampling_effort
    # --------------------------------------------------------
    grouped_data <- grouped_data %>%
      dplyr::mutate(
        basin_id   = basin_id_value,
        valid_name = species
      ) %>%
      dplyr::rename(sampling_year = year) %>%
      dplyr::select(
        basin_id,
        valid_name,
        sampling_year,
        sampling_effort
      )
    
    # --------------------------------------------------------
    # 3.6 Join grouped_data back into spp rows of this species,
    #     so that all original spp columns are preserved
    #     and we simply "add" sampling_year / sampling_effort
    # --------------------------------------------------------
    species_rows <- spp %>%
      dplyr::filter(valid_name == species) %>%
      dplyr::left_join(grouped_data, by = c("basin_id", "valid_name"))
    
    all_results <- dplyr::bind_rows(all_results, species_rows)
    
    cli_alert_info(
      paste("Done:", species, "in basin", basin_id_value)
    )
  }
  
  # ------------------------------------------------------------
  # 4. Save results for this basin (spp + sampling fields)
  # ------------------------------------------------------------
  output_file <- paste0(path_output_dir, basin_id_value, ".csv")
  write.csv(all_results, file = output_file, row.names = FALSE)
  
  cli_alert_success(
    paste("Finished basin", basin_id_value, "→ Saved to", output_file)
  )
  
  gc()
}



# ------------------------------------------------------------
# Run all basins
# ------------------------------------------------------------
purrr::walk(basin_ids, process_basin)

cli_progress_done()
gc()
