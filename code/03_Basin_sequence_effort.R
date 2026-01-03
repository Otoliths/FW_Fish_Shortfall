# Sequencing effort was calculated at the basin-species-year level, with unique individuals 
# (identified by specimen_id) as the core unit to avoid overcounting biases from multiple 
# gene-region submissions per individual. Each accession number was treated as a sequencing event, 
# representing one gene-region submission linked to a single individual via specimen_id. For each species, 
# we first quantified molecular coverage per individual per year by counting the number of distinct gene-regions 
# sequenced for each specimen_id within a given year. We then aggregated these values across all individuals to 
# derive the annual sequencing effort—this metric inherently integrates two key dimensions of genetic sampling: 
# (1) temporal intensity (number of unique individuals sequenced annually) and 
# (2) molecular coverage (number of distinct gene-regions per individual). Finally, 
# annual sequencing effort values were cumulatively summed across years to generate a time series of 
# sequencing effort for each species within each basin. This quantification reflects comprehensive 
# sequencing investment in freshwater fish taxonomy, eliminating biases caused by multiple accessions derived from the same individual.


# ============================================================
# Basin-level Sequencing Effort Calculation
# ------------------------------------------------------------
# For each basin, this script compiles NCBI sequencing records
# (specimen × gene types) of fish species, calculates annual
# and cumulative sequencing effort, and exports basin-level tables.
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
dir_species       <- "input/processed/basin_sp/"
dir_output        <- "input/data_prep/basin_sequencing_effort/"


# Ensure output directory exists
if (!dir.exists(dir_output)) {
  dir.create(dir_output, recursive = TRUE)
  cli_alert_success("Created directory: {dir_output}")
}
# ------------------------------------------------------------
# Load and prepare basin information
# ------------------------------------------------------------
inland <- readRDS("input/raw/basin/basin_sf_v1.rds")

seq_annotation <- read_parquet("input/raw/seq_annotation_v1.parquet")

# Extract all basin IDs
basin_ids <- inland %>%
  st_drop_geometry() %>%
  pull(basin_id)

# ------------------------------------------------------------
# Helper function: process sequencing effort for a single basin
# ------------------------------------------------------------
process_basin <- function(basin_id_value) {
  
  # ------------------------------------------------------------
  # Display the basin being processed
  # ------------------------------------------------------------
  cli_h2(paste("Processing basin:", basin_id_value))
  
  # ------------------------------------------------------------
  # Read the basin-level species list file
  # ------------------------------------------------------------
  spp_file <- paste0(dir_species, basin_id_value, ".csv")
  
  # If file does not exist, skip this basin
  if (!file.exists(spp_file)) {
    cli_alert_warning(paste("File not found:", spp_file, "- skipping"))
    return(NULL)
  }
  
  # Read species list and keep only required columns
  # (assuming column 1 is basin_id, column 3 is valid_name)
  spp <- read.csv(spp_file)
  spp <- spp[, c(1, 3)]
  
  # Validate species table format
  if (!is.data.frame(spp)) {
    cli_alert_warning("Invalid species file - skipping basin")
    return(NULL)
  }
  if (!"valid_name" %in% colnames(spp)) {
    cli_alert_warning("'valid_name' column missing - skipping basin")
    return(NULL)
  }
  
  # ------------------------------------------------------------
  # Store final results for all species in this basin
  # (one row per spp entry × sequencing year, if available)
  # ------------------------------------------------------------
  all_results <- tibble::tibble()
  
  # ------------------------------------------------------------
  # Loop through each species in the basin species list
  # ------------------------------------------------------------
  for (species in unique(spp$valid_name)) {
    
    # Filter GenBank annotation table for the current species
    gb_file <- seq_annotation %>%
      dplyr::filter(valid_name == species)
    
    # If no metadata exists for this species at all
    if (nrow(gb_file) == 0) {
      cli_alert_warning(
        paste("No sequence metadata found for species:", species,
              "- setting sequencing fields to NA")
      )
      
      # Append spp rows for this species with NA in the three new fields
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
    
    # Select required columns and remove duplicates, then filter by basin
    gb <- gb_file %>%
      dplyr::select(specimen_id, valid_name, basin_id, year, gene) %>%
      dplyr::distinct() %>%
      dplyr::filter(basin_id == basin_id_value)
    
    # If the species has no sequence records in this basin
    if (nrow(gb) == 0) {
      cli_alert_warning(
        paste("No sequence records for species in this basin:", species,
              "- setting sequencing fields to NA")
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
    # Compute sequencing effort for this species in this basin:
    # 1. Count unique gene types per specimen per year
    # 2. Summarize total gene-type counts per year
    # 3. Compute cumulative sequencing effort across years
    # ------------------------------------------------------------
    sequencing_effort <- gb %>%
      dplyr::filter(!is.na(year), !is.na(specimen_id), !is.na(gene)) %>%
      dplyr::group_by(year, specimen_id) %>%
      dplyr::summarise(
        gene_type_count = dplyr::n_distinct(gene),
        .groups = "drop"
      ) %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(
        sequencing_effort = sum(gene_type_count),
        .groups = "drop"
      ) %>%
      dplyr::arrange(year) %>%
      dplyr::mutate(
        valid_name = species,
        sequencing_effort_cumsum = cumsum(sequencing_effort)
      ) %>%
      dplyr::rename(sequencing_year = year) %>%
      dplyr::mutate(sequencing_year = as.integer(sequencing_year))
    
    # ------------------------------------------------------------
    # Merge sequencing results into the spp entry for this species
    # Only filter spp for the current species before joining,
    # so that previous species results are not overwritten.
    # ------------------------------------------------------------
    merged_result <- spp %>%
      dplyr::filter(valid_name == species) %>%
      dplyr::left_join(sequencing_effort, by = "valid_name")
    
    # Ensure sequencing_year is integer if present
    if ("sequencing_year" %in% names(merged_result)) {
      merged_result <- merged_result %>%
        dplyr::mutate(sequencing_year = as.integer(sequencing_year))
    }
    
    # Append current species result to the global result table
    all_results <- dplyr::bind_rows(all_results, merged_result)
    
    cli_alert_info(
      paste("Completed species:", species, "in basin:", basin_id_value)
    )
  }
  
  # ------------------------------------------------------------
  # Save results for this basin
  # ------------------------------------------------------------
  output_file <- paste0(dir_output, basin_id_value, ".csv")
  write.csv(all_results, file = output_file, row.names = FALSE)
  
  cli_alert_success(
    paste("Finished basin", basin_id_value, "→ Saved to", output_file)
  )
  
  # Cleanup memory
  gc()
}


# ------------------------------------------------------------
# Execute processing for all basins
# ------------------------------------------------------------
purrr::walk(basin_ids, process_basin)

cli_progress_done()
gc()
