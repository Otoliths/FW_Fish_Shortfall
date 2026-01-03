library(furrr)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(sf)
library(cli)
library(arrow)
library(magrittr)
# Suppress all warnings
options(warn = -1)

inland <- readRDS("input/raw/basin/basin_sf_v1.rds")

# Parallel plan (use all but one core)
future::plan(future::multisession(workers = parallel::detectCores() - 1))

seq_annotation <- read_parquet("input/raw/seq_annotation_v1.parquet")

# Extract basin IDs once
basin_id_value <- inland %>%
  st_drop_geometry() %>%
  pull(basin_id)

# Progress bar
cli_progress_bar("Processing basin ID", total = length(basin_id_value))

# Process a single basin_id
process_basin <- function(basin_id_value) {
  
  # Log current basin
  cli_h2(paste("Processing basin:", basin_id_value))
  
  # Species table for this basin
  spp_file <- paste0("input/processed/basin_sp/", basin_id_value, ".csv")
  
  # Skip if file missing
  if (!file.exists(spp_file)) {
    cli_alert_warning(paste("File not found:", spp_file, "- Skipping"))
    return(NULL)
  }
  
  # Read species data
  spp <- read.csv(spp_file)
  
  
  # Validate data frame
  if (!is.data.frame(spp)) {
    cli_alert_warning("spp is not a data frame - Skipping")
    return(NULL)
  }
  
  # Ensure 'valid_name' exists
  if (!"valid_name" %in% colnames(spp)) {
    cli_alert_warning("'valid_name' column missing - Skipping")
    return(NULL)
  }
  
  # Iterate over species in this basin
  for (species in unique(spp$valid_name)) {
    
    # Path to per-species NCBI/sequence metadata (CSV named by species)
    gb_file <- seq_annotation %>% dplyr::filter(valid_name %in% species)
    
    # If missing, set min_year to NA
    if (is_empty(gb_file)) {
      cli_alert_warning(paste("No sequence metadata found for species:", species, "- setting year_sequence = NA"))
      min_year <- NA
    } else {
      # Read NCBI/sequence data
      cli_alert_info("Reading NCBI data for species: {.bold {species}}")
      
      gb <- gb_file %>% 
        dplyr::select(valid_name, basin_id, year) %>%
        na.omit() %>%
        distinct() %>%
        dplyr::filter(basin_id %in% basin_id_value)
      
      # Year of description for this species (from spp table)
      current_year <- spp %>% 
        filter(valid_name == species) %>% 
        pull(year_description) %>% 
        unique()
      
      # Earliest sequence year later than description year
      min_year <- if (any(!is.na(gb$year) & gb$year > current_year)) {
        min(gb$year[gb$year > current_year], na.rm = TRUE)
      } else {
        NA
      }
    }
    
    # Ensure 'year_sequence' exists; then update for this species
    if (!"year_sequence" %in% colnames(spp)) {
      spp$year_sequence <- NA  # create if missing
    }
    
    spp %<>%
      mutate(
        year_sequence = ifelse(valid_name == species, min_year, year_sequence),
        basin_id = basin_id_value
      )
    
    # Optional column subset:
    # spp <- spp[,c("basin_id","biogeographic_realm","Ecoregion","valid_name","year_description","year_geolocation","year_sequence")]
    
    cli_alert_info("Completed species: {.bold {species}}, basin: {.bold {basin_id_value}}")
  }
  
  # Save updated species table (after all species are processed)
  write.csv(spp, file = spp_file, row.names = FALSE)
  
  cli_alert_info(paste("All species processed for basin: {.bold {basin_id_value}}"))
  
  # Clean up after each basin
  gc()
}

# Run over all basins
purrr::walk(basin_id_value, process_basin)

# Finish progress
cli_progress_done()

# Final cleanup (parallel session memory is minimal here, no plan teardown needed)
gc()
