# Amatulli, G. et al. A suite of global, cross-scale topographic variables for environmental and biodiversity modeling. Sci. Data 5, 180040 (2018).
# http://www.earthenv.org/topography

# Core libs
library(dplyr)
library(sf)
library(exactextractr)
library(raster)
library(purrr)
library(cli)
library(stringr)

# ------------------------------- Config ---------------------------------------
EQUAL_AREA_CRS <- "+proj=eck4 +lon_0=0 +lat_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
GRID_KM        <- 50
OCC_DIR        <- "input/processed/occ"
SPP_DIR        <- "input/processed/basin_sp"
OUT_DIR        <- "input/data_prep/basin_elevation"
ELEV_PATH      <- "input/raw/elevation/elevation_1KMmn_GMTEDmn.tif"

# Ensure output directory exists
if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
  cli_alert_success("Created directory: {OUT_DIR}")
}

# Load basins once, fix geometry/topology, and normalize across the dateline
inland <- readRDS("input/raw/basin/basin_sf_v1.rds") %>%
  st_transform(4326) %>%
  suppressWarnings(st_make_valid()) %>%
  # Wrap near the antimeridian to avoid huge polygons and self-intersections
  st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE) %>%
  st_collection_extract("POLYGON", warn = FALSE) %>%
  st_as_sf()

# Load elevation raster once (avoid repeated disk I/O in the loop)
elev_rast <- raster::raster(ELEV_PATH)

# ------------------------------ Helpers ---------------------------------------

# Coerce occurrence data to sf (accepts either sf with geometry or plain data with lon/lat)
as_occ_sf <- function(x, lon = "decimalLongitude", lat = "decimalLatitude", crs = 4326) {
  if (inherits(x, "sf")) {
    # Ensure CRS is set
    if (is.na(sf::st_crs(x))) x <- sf::st_set_crs(x, crs)
    return(x)
  }
  stopifnot(all(c(lon, lat) %in% names(x)))
  x %>%
    filter(!is.na(.data[[lon]]), !is.na(.data[[lat]])) %>%
    st_as_sf(coords = c(lon, lat), crs = crs, remove = FALSE)
}

# Build a clipped equal-area fishnet for a single basin (computed once per basin)
build_basin_grid <- function(basin_sf, grid_km = GRID_KM) {
  basin_eq <- st_transform(basin_sf, EQUAL_AREA_CRS)
  cell_m   <- rep(grid_km * 1000, 2)
  
  st_make_grid(basin_eq, cellsize = cell_m, what = "polygons") %>%
    suppressWarnings(st_make_valid()) %>%
    st_intersection(basin_eq) %>%
    st_as_sf() %>%
    mutate(grid_id = dplyr::row_number(),
           .before = 1) # keep an explicit ID up front
}

# Safe mean wrapper (NA-tolerant)
safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

# --------------------------- Main worker --------------------------------------

# Process a single basin_id (writes one CSV to OUT_DIR)
process_basin <- function(basin_id_value) {
  out_csv <- file.path(OUT_DIR, paste0(basin_id_value, ".csv"))
  
  # Skip if already processed
  if (file.exists(out_csv)) {
    cli_alert_info("File already exists: {out_csv} - skipping this basin")
    return(invisible(out_csv))
  }
  
  cli_h2("Processing basin: {basin_id_value}")
  
  # Load species table for this basin and validate
  spp_file <- file.path(SPP_DIR, paste0(basin_id_value, ".csv"))
  if (!file.exists(spp_file)) {
    cli_alert_warning("Species table not found: {spp_file} - skipping basin")
    return(NULL)
  }
  
  spp <- read.csv(spp_file, check.names = FALSE)
  
  if (!"valid_name" %in% names(spp)) {
    cli_alert_warning("'valid_name' column missing in {spp_file} - skipping basin")
    return(NULL)
  }
  
  # Keep a minimal, explicit set of columns (avoid fragile positional slicing)
  keep_cols <- intersect(c("basin_id","valid_name", "something_else"), names(spp))
  if (length(keep_cols) == 0L) {
    # Fallback: keep first 3 columns if nothing matches (as in the original script)
    keep_cols <- names(spp)[seq_len(min(3L, ncol(spp)))]
  }
  spp <- tibble::as_tibble(spp[, keep_cols, drop = FALSE])
  
  # Make sure an elevation column exists for output
  if (!"elevation" %in% names(spp)) spp$elevation <- NA_real_
  
  # Extract this basin geometry once
  basin <- inland %>% filter(.data$basin_id %in% basin_id_value)
  if (nrow(basin) == 0) {
    cli_alert_warning("No basin geometry found for basin_id: {basin_id_value} - skipping")
    return(NULL)
  }
  
  # Build and cache the basin grid once
  fishnet_eq <- build_basin_grid(basin, grid_km = GRID_KM)
  
  # Prepare a version of the grid reprojected to the raster CRS (used for exact_extract)
  fishnet_rast_crs <- st_transform(fishnet_eq, crs = raster::crs(elev_rast))
  
  # Iterate over unique species
  for (sp in unique(spp$valid_name)) {
    occ_path <- file.path(OCC_DIR, paste0(str_replace_all(sp, " ", "_"), ".rds"))
    
    if (!file.exists(occ_path)) {
      cli_alert_warning("Distribution file missing: {occ_path} - skipping species: {.bold {sp}}")
      next
    }
    
    cli_alert_info("Reading distribution for species: {.bold {sp}}")
    occ_raw <- readRDS(occ_path)
    
    # Coerce to sf and transform to the fishnet CRS used for counting points
    occ_sf <- tryCatch(as_occ_sf(occ_raw), error = function(e) NULL)
    if (is.null(occ_sf) || nrow(occ_sf) == 0) {
      cli_alert_warning("No valid occurrences for species: {.bold {sp}} - skipping")
      next
    }
    occ_eq <- st_transform(occ_sf, st_crs(fishnet_eq))
    
    # Count occurrences per grid, keep only grids that have >=1 point
    n_pts <- lengths(st_intersects(fishnet_eq, occ_eq))
    has_pts_idx <- which(n_pts > 0L)
    
    if (length(has_pts_idx) == 0L) {
      cli_alert_info("No grid cells with points for species: {.bold {sp}} in basin {basin_id_value}")
      # Ensure elevation stays NA for this species entry
      spp$elevation[spp$valid_name == sp] <- NA_real_
      next
    }
    
    # Exact extract mean elevation over only the occupied grid cells
    elev_vals <- exact_extract(
      elev_rast,
      fishnet_rast_crs[has_pts_idx, , drop = FALSE],
      "mean",
      progress = FALSE
    )
    
    # Aggregate to a single basin-level elevation for this species
    spp$elevation[spp$valid_name == sp] <- safe_mean(unlist(elev_vals))
    cli_alert_info("Done: {.bold {sp}} in basin {.bold {basin_id_value}} | elevation (mean over occupied cells) = {round(spp$elevation[spp$valid_name==sp], 2)}")
  }
  
  # Write basin-level table
  readr::write_csv(spp, out_csv)
  cli_alert_success("Saved basin-level result: {out_csv}")
  invisible(out_csv)
}

# ------------------------------- Run ------------------------------------------
basin_ids <- inland %>% st_drop_geometry() %>% pull(basin_id)

# Use walk for side-effect processing; change to basin_ids to process all
purrr::walk(basin_ids, process_basin)


