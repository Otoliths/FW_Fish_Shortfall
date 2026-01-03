# ============================================================
# Country-level Sequencing Effort Calculation
# ------------------------------------------------------------
# For each Country, this script compiles NCBI sequencing records
# (specimen × gene types) of fish species, calculates annual
# and cumulative sequencing effort, and exports Country-level tables.
# ============================================================

library(furrr)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(sf)
library(arrow)
library(cli)

# ------------------------------------------------------------
# Global Settings
# ------------------------------------------------------------
options(warn = -1)  # Suppress warnings

# Directory configuration (centralized for easy modification)
dir_species       <- "input/processed/country_sp/"
dir_output        <- "input/data_prep/country_sequencing_effort/"


# Ensure output directory exists
if (!dir.exists(dir_output)) {
  dir.create(dir_output, recursive = TRUE)
  cli_alert_success("Created directory: {dir_output}")
}
# ------------------------------------------------------------
# Load and prepare country polygons
# ------------------------------------------------------------
country <- readRDS("input/processed/country_fix.rds") %>%
  st_transform(4326) %>%
  suppressWarnings(st_make_valid()) %>%
  # Wrap around the dateline to avoid invalid geometries across 180°
  st_wrap_dateline(
    options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
    quiet   = TRUE
  ) %>%
  st_collection_extract("POLYGON", warn = FALSE) %>%
  st_as_sf()

# Extract country IDs once
country_ids <- country %>%
  st_drop_geometry() %>%
  pull(iso3)

seq_annotation <- read_parquet("input/raw/seq_annotation_v1.parquet")


# ------------------------------------------------------------
# Helper function: process sequencing effort for a single country
# ------------------------------------------------------------
process_country <- function(country_id_value) {
  
  # ------------------------------------------------------------
  # Display the country being processed
  # ------------------------------------------------------------
  cli_h2(paste("Processing country:", country_id_value))
  
  # ------------------------------------------------------------
  # Read the country-level species list file
  # ------------------------------------------------------------
  spp_file <- paste0(dir_species, country_id_value, ".csv")
  
  # If file does not exist, skip this country
  if (!file.exists(spp_file)) {
    cli_alert_warning(paste("File not found:", spp_file, "- skipping"))
    return(NULL)
  }
  
  # Read species list and keep only required columns
  spp <- read.csv(spp_file)
  spp <- spp[, c(1, 2)]
  
  # Validate species table format
  if (!is.data.frame(spp)) {
    cli_alert_warning("Invalid species file - skipping country")
    return(NULL)
  }
  if (!"valid_name" %in% colnames(spp)) {
    cli_alert_warning("'valid_name' column missing - skipping country")
    return(NULL)
  }
  
  # ------------------------------------------------------------
  # Store final results for all species in this country
  # (One row per species × sequencing record)
  # ------------------------------------------------------------
  all_results <- tibble()
  
  # ------------------------------------------------------------
  # Loop through each species in the species list
  # ------------------------------------------------------------
  for (species in unique(spp$valid_name)) {
    
    # Filter GenBank annotation table for current species
    gb_file <- seq_annotation %>% 
      dplyr::filter(valid_name == species)
    
    # If no metadata exists for this species
    if (nrow(gb_file) == 0) {
      cli_alert_warning(
        paste("No sequence metadata found for species:", species,
              "- setting sequencing fields to NA")
      )
      
      # Append one row for this species with NA values
      merged_result <- spp %>%
        dplyr::filter(valid_name == species) %>%
        dplyr::mutate(
          sequencing_year          = NA_integer_,
          sequencing_effort        = NA_real_,
          sequencing_effort_cumsum = NA_real_
        )
      
      all_results <- dplyr::bind_rows(all_results, merged_result)
      next
    }
    
    cli_alert_info(paste("Reading NCBI metadata for species:", species))
    
    # Select required columns and remove duplicates
    gb <- gb_file %>%
      dplyr::select(specimen_id, valid_name, iso3, year, gene) %>%
      dplyr::distinct() %>%
      dplyr::filter(iso3 == country_id_value)
    
    # If the species has no sequences recorded for this country
    if (nrow(gb) == 0) {
      cli_alert_warning(
        paste("No sequence records for species in this country:", species,
              "- setting fields to NA")
      )
      
      merged_result <- spp %>%
        dplyr::filter(valid_name == species) %>%
        dplyr::mutate(
          sequencing_year          = NA_integer_,
          sequencing_effort        = NA_real_,
          sequencing_effort_cumsum = NA_real_
        )
      
      all_results <- dplyr::bind_rows(all_results, merged_result)
      next
    }
    
    # ------------------------------------------------------------
    # Compute sequencing effort:
    # 1. Count unique gene types per specimen per year
    # 2. Summarize total gene-type counts per year
    # 3. Compute cumulative sequencing effort
    # ------------------------------------------------------------
    sequencing_effort <- gb %>%
      dplyr::filter(!is.na(year), !is.na(specimen_id), !is.na(gene)) %>%
      dplyr::group_by(year, specimen_id) %>%
      dplyr::summarise(gene_type_count = dplyr::n_distinct(gene), .groups = "drop") %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(sequencing_effort = sum(gene_type_count), .groups = "drop") %>%
      dplyr::arrange(year) %>%
      dplyr::mutate(
        valid_name = species,
        sequencing_effort_cumsum = cumsum(sequencing_effort)
      ) %>%
      dplyr::rename(sequencing_year = year) %>%
      dplyr::mutate(sequencing_year = as.integer(sequencing_year))
    
    # ------------------------------------------------------------
    # Merge sequencing results into the spp entry for this species
    # Only join the rows for the current species
    # (so previous species results are not overwritten)
    # ------------------------------------------------------------
    merged_result <- spp %>%
      dplyr::filter(valid_name == species) %>%
      dplyr::left_join(sequencing_effort, by = "valid_name")
    
    # Ensure sequencing_year is integer if present
    if ("sequencing_year" %in% names(merged_result)) {
      merged_result <- merged_result %>%
        dplyr::mutate(sequencing_year = as.integer(sequencing_year))
    }
    
    # Append species result to global result table
    all_results <- dplyr::bind_rows(all_results, merged_result)
    
    cli_alert_info(
      paste("Completed species:", species, "for country:", country_id_value)
    )
  }
  
  # ------------------------------------------------------------
  # Save results for this country
  # ------------------------------------------------------------
  output_file <- paste0(dir_output, country_id_value, ".csv")
  write.csv(all_results, file = output_file, row.names = FALSE)
  
  cli_alert_success(
    paste("Finished country", country_id_value, "→ Saved to", output_file)
  )
  
  # Cleanup memory
  gc()
}

# ------------------------------------------------------------
# Execute processing for all country
# ------------------------------------------------------------
purrr::walk(country_ids, process_country)

cli_progress_done()
gc()
