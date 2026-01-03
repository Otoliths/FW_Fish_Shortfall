library(dplyr)
library(purrr)
library(stringr)
library(sf)
library(cli)

options(warn = -1)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
path_species_dir   <- "input/processed/country_sp/"
path_occ_dir       <- "input/processed/occ/"
path_output_dir    <- "input/data_prep/country_latitude/"

if (!dir.exists(path_output_dir)) {
  dir.create(path_output_dir, recursive = TRUE)
  cli_alert_success("Created directory: {path_output_dir}")
}

# ------------------------------------------------------------
# country (global)
# ------------------------------------------------------------
country <- readRDS("input/processed/country_fix.rds") %>%
  st_transform(4326) %>%
  suppressWarnings(st_make_valid()) %>%
  st_wrap_dateline(
    options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
    quiet = TRUE
  ) %>%
  st_collection_extract("POLYGON", warn = FALSE) %>%
  st_as_sf()

country_ids <- country %>%
  st_drop_geometry() %>%
  pull(iso3)

cli_progress_bar("Processing country", total = length(country_ids))

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
# Per-country processing
# ------------------------------------------------------------
process_country <- function(country_id_value) {
  
  cli_h2(paste("Processing country:", country_id_value))
  
  spp_file <- paste0(path_species_dir, country_id_value, ".csv")
  if (!file.exists(spp_file)) {
    cli_alert_warning(paste("File not found:", spp_file, "- skipping"))
    return(NULL)
  }
  
  # keep species name + valid_name
  spp <- read.csv(spp_file)[, c(1,2)]
  
  if (!is.data.frame(spp)) {
    cli_alert_warning("Species table is not a data frame - skipping")
    return(NULL)
  }
  if (!"valid_name" %in% colnames(spp)) {
    cli_alert_warning("'valid_name' column missing - skipping")
    return(NULL)
  }
  
  # prepare output column
  spp$country_latitude <- NA_real_
  
  # country polygon in WGS84
  country <- country %>%
    filter(iso3 %in% country_id_value) %>%
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
    
    # Keep only points spatially within the country polygon
    occ_in_country <- st_join(occ_sf, country, join = st_within, left = FALSE)
    if (nrow(occ_in_country) == 0) {
      cli_alert_warning("No occurrence points in country for {species}")
      next
    }
    
    # Step 1: strict deduplicate exact same coords
    occ_unique <- occ_in_country %>%
      distinct(decimalLongitude, decimalLatitude, .keep_all = TRUE)
    
    # Step 2: 1 km thinning to avoid oversampling at monitoring stations
    occ_thinned <- thin_by_distance_1km(occ_unique, radius_m = 1000)
    
    if (nrow(occ_thinned) == 0) {
      # fallback: if thinning somehow killed everything, use unique set
      occ_thinned <- occ_unique
    }
    
    # Compute mean latitude using the thinned points
    this_country_latitude <- mean(occ_thinned$decimalLatitude, na.rm = TRUE)
    
    spp <- spp %>%
      mutate(
        country_latitude = ifelse(
          valid_name == species,
          this_country_latitude,
          country_latitude
        )
      )
    
    cli_alert_info(
      "Finished {.bold {species}} in {.bold {country_id_value}} | country_latitude = {round(this_country_latitude, 3)}"
    )
  }
  
  # save country result
  out_file <- paste0(path_output_dir, country_id_value, ".csv")
  write.csv(spp, file = out_file, row.names = FALSE)
  
  cli_alert_success(
    "Completed country {.bold {country_id_value}} → saved to {.path {out_file}}"
  )
  
  gc()
}

# ------------------------------------------------------------
# Run across country
# ------------------------------------------------------------
purrr::walk(country_ids, process_country)

cli_progress_done()
gc()



