# Amatulli, G. et al. A suite of global, cross-scale topographic variables
# for environmental and biodiversity modeling. Sci. Data 5, 180040 (2018).
# http://www.earthenv.org/topography

# ------------------------------- Libraries -----------------------------------
library(dplyr)
library(sf)
library(exactextractr)
library(raster)
library(purrr)
library(cli)
library(stringr)
library(readr)

# Global options
options(warn = -1)        # suppress warnings
options(sf_use_s2 = FALSE)

# ------------------------------- Config --------------------------------------
EQUAL_AREA_CRS <- "+proj=eck4 +lon_0=0 +lat_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
GRID_KM        <- 50
OCC_DIR        <- "input/processed/occ"
SPP_DIR        <- "input/processed/country_sp"
OUT_DIR        <- "input/data_prep/country_elevation"
ELEV_PATH      <- "input/raw/elevation/elevation_1KMmn_GMTEDmn.tif"

# Ensure output directory exists
if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
  cli_alert_success("Created directory: {OUT_DIR}")
} else {
  cli_alert_info("Output directory exists: {OUT_DIR}")
}

# ------------------------------- Data Load -----------------------------------
# Load country polygons once
country <- readRDS("input/processed/country_fix.rds") %>%
  st_transform(4326) %>%
  suppressWarnings(st_make_valid()) %>%
  # Normalize geometries that cross the date line, split big MULTIPOLYGONs cleanly
  st_wrap_dateline(
    options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
    quiet   = TRUE
  ) %>%
  st_collection_extract("POLYGON", warn = FALSE) %>%
  st_as_sf()

# Load elevation raster once (avoid repeated disk I/O in loop)
elev_rast <- raster::raster(ELEV_PATH)

# ------------------------------ Helpers --------------------------------------

# Convert occurrence data (either sf or lon/lat table) into sf with CRS
as_occ_sf <- function(x,
                      lon = "decimalLongitude",
                      lat = "decimalLatitude",
                      crs = 4326) {
  if (inherits(x, "sf")) {
    # If CRS missing, assign one
    if (is.na(sf::st_crs(x))) x <- sf::st_set_crs(x, crs)
    return(x)
  }
  stopifnot(all(c(lon, lat) %in% names(x)))
  
  x %>%
    filter(!is.na(.data[[lon]]), !is.na(.data[[lat]])) %>%
    st_as_sf(coords = c(lon, lat), crs = crs, remove = FALSE)
}

# Safe mean wrapper (NA-tolerant → returns NA if all NA)
safe_mean <- function(x) {
  if (all(is.na(x))) {
    NA_real_
  } else {
    mean(x, na.rm = TRUE)
  }
}

# Build a clipped equal-area grid (fishnet) for ONE country
# Returns: sf polygons with grid_id, or NULL on failure
build_country_grid <- function(country_sf, grid_km = GRID_KM) {
  # country_sf is assumed to be that country's polygon(s) still in lon/lat (EPSG:4326)
  
  # 1. Clean and dissolve country geometry, then project to equal-area
  country_eq <- country_sf %>%
    st_make_valid() %>%
    suppressWarnings(st_buffer(0)) %>%
    st_union() %>%                    # dissolve multi-island boundaries
    st_cast("MULTIPOLYGON") %>%
    st_transform(EQUAL_AREA_CRS) %>%
    st_make_valid() %>%
    suppressWarnings(st_buffer(0))
  
  # 2. Build regular grid in equal-area CRS
  cell_m <- rep(grid_km * 1000, 2)
  
  grid_raw <- st_make_grid(
    country_eq,
    cellsize = cell_m,
    what     = "polygons"
  )
  
  # 3. Clip the grid to the country's shape
  fishnet_sf <- tryCatch({
    st_as_sf(grid_raw) %>%
      st_make_valid() %>%
      suppressWarnings(st_buffer(0)) %>%
      st_intersection(country_eq) %>%      # <- this is where Fiji used to crash
      st_make_valid() %>%
      suppressWarnings(st_buffer(0)) %>%
      mutate(
        grid_id = dplyr::row_number(),
        .before = 1
      )
  }, error = function(e) {
    cli_alert_danger(
      "Grid intersection failed for this country: {conditionMessage(e)}"
    )
    return(NULL)
  })
  
  fishnet_sf
}

# --------------------------- Main worker --------------------------------------
# Process a single country_id (writes one CSV to OUT_DIR)
process_country <- function(country_id_value) {
  
  out_csv <- file.path(OUT_DIR, paste0(country_id_value, ".csv"))
  
  # Skip if already processed
  if (file.exists(out_csv)) {
    cli_alert_info("File already exists: {out_csv} - skipping {country_id_value}")
    return(invisible(out_csv))
  }
  
  cli_h2("Processing country: {country_id_value}")
  
  # ---------------- species table for this country ----------------
  spp_file <- file.path(SPP_DIR, paste0(country_id_value, ".csv"))
  if (!file.exists(spp_file)) {
    cli_alert_warning("Species table not found: {spp_file} - skipping {country_id_value}")
    return(NULL)
  }
  
  spp <- read.csv(spp_file, check.names = FALSE)
  
  if (!"valid_name" %in% names(spp)) {
    cli_alert_warning("'valid_name' column missing in {spp_file} - skipping {country_id_value}")
    return(NULL)
  }
  
  # keep first 2 columns (you did spp <- spp[,1:2])
  spp <- spp[, 1:2]
  # ensure elevation column exists
  if (!"elevation" %in% names(spp)) {
    spp$elevation <- NA_real_
  }
  
  # ---------------- extract geometry for THIS country ----------------
  country_this <- country %>%
    dplyr::filter(.data$iso3 == country_id_value)
  
  if (nrow(country_this) == 0) {
    cli_alert_warning(
      "No country geometry found for {country_id_value} - skipping"
    )
    return(NULL)
  }
  
  # ---------------- build the fishnet grid once ----------------
  fishnet_eq <- build_country_grid(country_this, grid_km = GRID_KM)
  
  if (is.null(fishnet_eq) || nrow(fishnet_eq) == 0) {
    cli_alert_warning(
      "Failed to build grid for {country_id_value} - skipping"
    )
    return(NULL)
  }
  
  # Reproject grid to raster CRS for extraction
  fishnet_rast_crs <- fishnet_eq %>%
    st_transform(crs = raster::crs(elev_rast)) %>%
    st_make_valid() %>%
    suppressWarnings(st_buffer(0))
  
  # ---------------- iterate over species ----------------
  unique_species <- unique(spp$valid_name)
  
  for (sp in unique_species) {
    
    occ_path <- file.path(
      OCC_DIR,
      paste0(str_replace_all(sp, " ", "_"), ".rds")
    )
    
    if (!file.exists(occ_path)) {
      cli_alert_warning(
        "Distribution file missing: {occ_path} - skipping species: {.bold {sp}}"
      )
      next
    }
    
    cli_alert_info("Reading distribution for species: {.bold {sp}}")
    
    occ_raw <- readRDS(occ_path)
    
    # coerce occurrence to sf
    occ_sf <- tryCatch(
      as_occ_sf(occ_raw),
      error = function(e) NULL
    )
    
    if (is.null(occ_sf) || nrow(occ_sf) == 0) {
      cli_alert_warning(
        "No valid occurrences for species {.bold {sp}} - skipping"
      )
      next
    }
    
    # transform occurrences to equal-area CRS used by fishnet_eq
    occ_eq <- st_transform(occ_sf, st_crs(fishnet_eq))
    
    # how many occurrence points fall in each grid cell?
    n_pts <- lengths(st_intersects(fishnet_eq, occ_eq))
    has_pts_idx <- which(n_pts > 0L)
    
    if (length(has_pts_idx) == 0L) {
      cli_alert_info(
        "No grid cells with points for {.bold {sp}} in country {country_id_value}"
      )
      # elevation stays NA for this species
      spp$elevation[spp$valid_name == sp] <- NA_real_
      next
    }
    
    # extract mean elevation only for occupied grid cells
    elev_vals <- tryCatch({
      exact_extract(
        elev_rast,
        fishnet_rast_crs[has_pts_idx, , drop = FALSE],
        "mean",
        progress = FALSE
      )
    }, error = function(e) {
      cli_alert_warning(
        "exact_extract() failed for species {.bold {sp}} in {country_id_value}: {conditionMessage(e)}"
      )
      list(NA_real_)
    })
    
    # assign species-level mean elevation (mean of grid means)
    spp$elevation[spp$valid_name == sp] <- safe_mean(unlist(elev_vals))
    
    cli_alert_info(
      "Done: {.bold {sp}} in {.bold {country_id_value}} | elevation (mean over occupied cells) = {round(spp$elevation[spp$valid_name==sp], 2)}"
    )
  }
  
  # ---------------- write result ----------------
  readr::write_csv(spp, out_csv)
  cli_alert_success("Saved country-level result: {out_csv}")
  
  invisible(out_csv)
}

# ------------------------------- Run ------------------------------------------
# list of country ISO3 codes from the country sf
country_ids <- country %>%
  st_drop_geometry() %>%
  pull(iso3)

# Run for all countries:
purrr::walk(country_ids, process_country)


