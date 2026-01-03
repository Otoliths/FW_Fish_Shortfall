library(furrr)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(sf)
library(cli)
library(magrittr)

# Suppress all warnings
options(warn = -1)

inland <- readRDS("input/raw/basin/basin_sf_v1.rds")

# Configure parallel processing (single basin here, but keep parallel setup)
future::plan(future::multisession(workers = parallel::detectCores() - 30))  # Use available cores minus 30

# Get basin IDs
basin_id_value <- inland %>%
  st_drop_geometry() %>%
  pull(basin_id)

###########################-----------------------------------------------------

# Initialize progress bar
cli_progress_bar("Processing basin ID", total = length(basin_id_value))

# Function to process a single basin_id
process_basin <- function(basin_id_value) {
  
  # Log current basin being processed
  cli_h2(paste("Processing basin:", basin_id_value))
  
  # Load species data for the basin
  spp_file <- paste0("input/processed/basin_sp/", basin_id_value, ".csv")
  
  # Skip if the file does not exist
  if (!file.exists(spp_file)) {
    cli_alert_warning(paste("File not found:", spp_file, "- Skipping"))
    return(NULL)
  }
  
  # Read species data
  spp <- read.csv(spp_file)
  
  # Ensure the object is a data frame
  if (!is.data.frame(spp)) {
    cli_alert_warning("spp is not a data frame - Skipping")
    return(NULL)
  }
  
  # Ensure 'valid_name' column exists
  if (!"valid_name" %in% colnames(spp)) {
    cli_alert_warning("Missing 'valid_name' column - Skipping")
    return(NULL)
  }
  
  # Iterate through all unique species names
  for (species in unique(spp$valid_name)) {
    
    # Get basin geometry
    basin <- inland %>% filter(basin_id %in% basin_id_value)
    
    # Build path to species occurrence file (replace spaces with underscores)
    distribution_file <- paste0("input/processed/occ/", gsub(" ", "_", species), ".rds")
    
    # Handle missing occurrence file
    if (!file.exists(distribution_file)) {
      cli_alert_warning(paste("Distribution file missing:", distribution_file, "- setting year_geolocation = NA"))
      min_year <- NA
    } else {
      # Read occurrence data
      cli_alert_info("Reading distribution data for species: {.bold {species}}")
      
      distribution_sf <- readRDS(distribution_file) %>%
        st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
      
      # Spatially intersect occurrence points with the basin polygon
      sf_points_in_basin <- st_intersection(distribution_sf, basin)
      
      # Extract year of description for the species
      current_year <- spp %>% 
        filter(valid_name == species) %>% 
        pull(year_description) %>% 
        unique()
      
      # Compute the earliest georeferenced record later than the description year
      min_year <- if (any(!is.na(sf_points_in_basin$year) & sf_points_in_basin$year > current_year)) {
        min(sf_points_in_basin$year[sf_points_in_basin$year > current_year], na.rm = TRUE)
      } else {
        NA
      }
    }
    
    # Ensure 'year_geolocation' column exists and update it
    if (!"year_geolocation" %in% colnames(spp)) {
      spp$year_geolocation <- NA
    }
    
    spp %<>%
      mutate(year_geolocation = ifelse(valid_name == species, min_year, year_geolocation))
    
    cli_alert_info("Completed species: {.bold {species}} in basin: {.bold {basin_id_value}}")
  }
  
  # Save updated species data
  write.csv(spp, file = spp_file, row.names = FALSE)
  
  cli_alert_info(paste("All species processed for basin: {.bold {basin_id_value}}"))
  
  # Perform garbage collection after each basin
  gc()
}

# Process all basin IDs
purrr::walk(basin_id_value, process_basin)

# Finish progress bar
cli_progress_done()

# Shut down parallel workers and free memory
future::plan(future::sequential)
gc()








