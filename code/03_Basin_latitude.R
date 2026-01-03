# For each species–basin combination we quantified its latitudinal position (“basin_latitude”) as the mean latitude of that species’ occurrence records within the hydrological basin.
# We first filtered all occurrence points to those falling inside the basin polygon and removed duplicate records with identical coordinates.
# To reduce oversampling bias from repeatedly monitored river reaches (e.g. fish passages, survey stations), we then applied 1 km spatial thinning: within each basin, no two retained records were allowed to be within 1 km of each other (distance calculated in projected metric space).
# The resulting thinned set of occurrence coordinates was used to compute the basin-level mean latitude.

# We projected all occurrence coordinates to the Eckert IV equal-area projection (EPSG:54012 equivalent; +proj=eck4) to ensure consistent metric distance calculations across latitudes.
# Spatial thinning was then applied at a 1 km radius (Euclidean distance in projected meters) to reduce oversampling bias while preserving global comparability.

library(dplyr)
library(purrr)
library(stringr)
library(sf)
library(cli)
sf::sf_use_s2(F) 
options(warn = -1)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
path_species_dir   <- "input/processed/basin_sp/"
path_occ_dir       <- "input/processed/occ/"
path_output_dir    <- "input/data_prep/basin_latitude/"

if (!dir.exists(path_output_dir)) {
  dir.create(path_output_dir, recursive = TRUE)
  cli_alert_success("Created directory: {path_output_dir}")
}

# ------------------------------------------------------------
# Basins (global)
# ------------------------------------------------------------
inland <- readRDS("input/raw/basin/basin_sf_v1.rds") %>%
  st_transform(4326) %>%
  suppressWarnings(st_make_valid()) %>%
  st_wrap_dateline(
    options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
    quiet = TRUE
  ) %>%
  st_collection_extract("POLYGON", warn = FALSE) %>%
  st_as_sf()

basin_ids <- inland %>%
  st_drop_geometry() %>%
  pull(basin_id)


cli_progress_bar("Processing basins", total = length(basin_ids))

# ------------------------------------------------------------
# Helper: 1 km spatial thinning (greedy)
# Input: sf POINT layer in EPSG:4326, with columns decimalLongitude/decimalLatitude
# Output: subset of rows, thinned so that no two kept points are within ~1 km
# ------------------------------------------------------------
thin_by_distance_1km <- function(points_sf_wgs84, radius_m = 1000) {
  # project to equal-area Eckert IV (used consistently in this project)
  pts_m <- st_transform(
    points_sf_wgs84,
    crs = "+proj=eck4 +lon_0=0 +lat_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  )
  
  coords <- st_coordinates(pts_m)
  
  keep_idx <- logical(nrow(pts_m))
  kept_coords <- matrix(NA_real_, 0, 2)
  
  for (i in seq_len(nrow(pts_m))) {
    this_xy <- coords[i, , drop = FALSE]
    if (nrow(kept_coords) == 0) {
      keep_idx[i] <- TRUE
      kept_coords <- this_xy
    } else {
      d2 <- (kept_coords[,1] - this_xy[,1])^2 + (kept_coords[,2] - this_xy[,2])^2
      dmin <- sqrt(min(d2))
      if (dmin > radius_m) {
        keep_idx[i] <- TRUE
        kept_coords <- rbind(kept_coords, this_xy)
      } else {
        keep_idx[i] <- FALSE
      }
    }
  }
  
  points_sf_wgs84[keep_idx, ]
}


# ------------------------------------------------------------
# Per-basin processing
# ------------------------------------------------------------
process_basin <- function(basin_id_value) {
  
  cli_h2(paste("Processing basin:", basin_id_value))
  
  spp_file <- paste0(path_species_dir, basin_id_value, ".csv")
  if (!file.exists(spp_file)) {
    cli_alert_warning(paste("File not found:", spp_file, "- skipping"))
    return(NULL)
  }
  
  # keep species name + valid_name
  spp <- read.csv(spp_file)[, c(1, 3)]
  
  if (!is.data.frame(spp)) {
    cli_alert_warning("Species table is not a data frame - skipping")
    return(NULL)
  }
  if (!"valid_name" %in% colnames(spp)) {
    cli_alert_warning("'valid_name' column missing - skipping")
    return(NULL)
  }
  
  # prepare output column
  spp$basin_latitude <- NA_real_
  
  # basin polygon in WGS84
  basin <- inland %>%
    filter(basin_id %in% basin_id_value) %>%
    st_make_valid()
  
  # loop species
  for (species in unique(spp$valid_name)) {
    
    occ_file <- paste0(
      path_occ_dir,
      gsub(" ", "_", species),
      ".rds"
    )
    
    if (!file.exists(occ_file)) {
      cli_alert_warning(paste("Occurrence file not found:", occ_file, " - skipping {species}"))
      next
    }
    
    cli_alert_info("Reading occurrence data for species: {.bold {species}}")
    
    occ_raw <- readRDS(occ_file)
    
    if (!all(c("decimalLongitude", "decimalLatitude") %in% colnames(occ_raw))) {
      cli_alert_warning("Occurrence data lacks coordinates - skipping {species}")
      next
    }
    
    # Make sf (WGS84)
    occ_sf <- st_as_sf(
      occ_raw,
      coords = c("decimalLongitude", "decimalLatitude"),
      crs = 4326,
      remove = FALSE
    )
    
    # Keep only points spatially within the basin polygon
    occ_in_basin <- st_join(occ_sf, basin, join = st_within, left = FALSE)
    if (nrow(occ_in_basin) == 0) {
      cli_alert_warning("No occurrence points in basin for {species}")
      next
    }
    
    # Step 1: strict deduplicate exact same coords
    occ_unique <- occ_in_basin %>%
      distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE)
    
    # Step 2: 1 km thinning to avoid oversampling at monitoring stations
    occ_thinned <- thin_by_distance_1km(occ_unique, radius_m = 1000)
    
    if (nrow(occ_thinned) == 0) {
      # fallback: if thinning somehow killed everything, use unique set
      occ_thinned <- occ_unique
    }
    
    # Compute mean latitude using the thinned points
    this_basin_latitude <- mean(occ_thinned$decimalLatitude, na.rm = TRUE)
    
    spp <- spp %>%
      mutate(
        basin_latitude = ifelse(
          valid_name == species,
          this_basin_latitude,
          basin_latitude
        )
      )
    
    cli_alert_info(
      "Finished {.bold {species}} in {.bold {basin_id_value}} | basin_latitude = {round(this_basin_latitude, 3)}"
    )
  }
  
  # save basin result
  out_file <- paste0(path_output_dir, basin_id_value, ".csv")
  write.csv(spp, file = out_file, row.names = FALSE)
  
  cli_alert_success(
    "Completed basin {.bold {basin_id_value}} → saved to {.path {out_file}}"
  )
  
  gc()
}

# ------------------------------------------------------------
# Run across basins
# ------------------------------------------------------------
purrr::walk(basin_ids, process_basin)

cli_progress_done()
gc()
