library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(sf)
library(ggplot2)
library(cli)

# Suppress all warnings
options(warn = -1)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
path_species_dir   <- "input/processed/basin_sp/"
path_occ_dir       <- "input/processed/occ/"
path_output_dir    <- "input/data_prep/basin_range_size/"

# Make sure output directory exists
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
  # Wrap around the dateline to avoid huge self-crossing polygons
  st_wrap_dateline(
    options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
    quiet = TRUE
  ) %>%
  st_collection_extract("POLYGON", warn = FALSE) %>%
  st_as_sf()

# Get all basin IDs
basin_ids <- inland %>%
  st_drop_geometry() %>%
  pull(basin_id)


# Progress bar for all basins
cli_progress_bar("Processing basins", total = length(basin_ids))

# ------------------------------------------------------------
# Function: process a single basin
# ------------------------------------------------------------
process_basin <- function(basin_id_value) {
  
  cli_h2(paste("Processing basin:", basin_id_value))
  
  # Path to species list for this basin
  spp_file <- paste0(path_species_dir, basin_id_value, ".csv")
  
  # Skip this basin if species list is missing
  if (!file.exists(spp_file)) {
    cli_alert_warning(paste("File not found:", spp_file, "- skipping"))
    return(NULL)
  }
  
  # Read species table and keep relevant columns
  spp <- read.csv(spp_file)[, c(1, 3)]
  
  # Basic validation
  if (!is.data.frame(spp)) {
    cli_alert_warning("Species table is not a data frame - skipping")
    return(NULL)
  }
  if (!"valid_name" %in% colnames(spp)) {
    cli_alert_warning("'valid_name' column missing - skipping")
    return(NULL)
  }
  
  # Ensure output column exists
  if (!"range_size" %in% colnames(spp)) {
    spp$range_size <- NA_real_
  }
  
  # Loop over all species in this basin
  for (species in unique(spp$valid_name)) {
    
    # Get the geometry for the current basin
    basin <- inland %>% filter(basin_id %in% basin_id_value)
    
    # Path to occurrence data for this species
    distribution_file <- paste0(
      path_occ_dir,
      gsub(" ", "_", species),
      ".rds"
    )
    
    # Skip species if no occurrence data
    if (!file.exists(distribution_file)) {
      cli_alert_warning(paste("Occurrence file not found:", distribution_file, "- skipping"))
      next
    }
    
    cli_alert_info("Reading occurrence data for species: {.bold {species}}")
    
    # Load occurrences and convert to sf
    distribution_sf <- readRDS(distribution_file) %>%
      st_as_sf(
        coords = c("decimalLongitude", "decimalLatitude"),
        crs = 4326
      )
    
    # Reproject basin to an equal-area projection (Eckert IV)
    basin_proj <- st_transform(
      basin,
      crs = "+proj=eck4 +lon_0=0 +lat_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    )
    
    # Build a 50 km grid over the basin and clip it to the basin boundary
    fishnet_clipped_sf <- st_make_grid(
      basin_proj,
      cellsize = c(50000, 50000),
      what = "polygons"
    ) %>%
      suppressWarnings(st_make_valid()) %>%
      st_intersection(basin_proj) %>%
      st_as_sf() %>%
      mutate(grid_id = row_number())
    
    # Reproject occurrences to match the grid
    points_sf_in_basin <- st_transform(distribution_sf, crs = st_crs(fishnet_clipped_sf))
    
    # Count points per grid cell
    fishnet_clipped_sf$n_colli <- lengths(
      st_intersects(fishnet_clipped_sf, points_sf_in_basin)
    )
    
    # Keep only grid cells that have at least one record
    fishnet_count <- fishnet_clipped_sf %>%
      filter(n_colli > 0)
    
    # Update range_size for this species:
    # use the number of occupied grid cells as a proxy for range size
    spp <- spp %>%
      mutate(
        range_size = ifelse(valid_name == species, nrow(fishnet_count), range_size)
      )
    
    
    cli_alert_info(
      "Finished species {.bold {species}} in basin {.bold {basin_id_value}}"
    )
  }
  
  # Write basin-level results
  out_file <- paste0(path_output_dir, basin_id_value, ".csv")
  write.csv(spp, file = out_file, row.names = FALSE)
  
  cli_alert_success(
    "All species processed for basin {.bold {basin_id_value}} → saved to {.path {out_file}}"
  )
  
  # Clean up memory for large loops
  gc()
}

# ------------------------------------------------------------
# Run for all basins
# ------------------------------------------------------------
purrr::walk(basin_ids, process_basin)

cli_progress_done()
gc()
