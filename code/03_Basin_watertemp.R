# ------------------------------- Libraries ------------------------------------
library(dplyr)
library(sf)
library(exactextractr)
library(terra)        # use consistently for rasters
library(purrr)
library(tidyr)
library(stringr)
library(ggplot2)
library(cli)

# ------------------------------- Config ---------------------------------------
EQUAL_AREA_CRS <- "+proj=eck4 +lon_0=0 +lat_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
GRID_KM        <- 50
DIR_OCC        <- "input/processed/occ"
DIR_SPP        <- "input/processed/basin_sp"
DIR_OUT        <- "input/data_prep/basin_watertemp"
DIR_RAST       <- "input/processed/waterTemp_annual_mean_ensemble"

# Ensure output directory exists
if (!dir.exists(DIR_OUT)) {
  dir.create(DIR_OUT, recursive = TRUE)
  cli_alert_info("Created directory: {DIR_OUT}")
} else {
  cli_alert_success("Directory already exists: {DIR_OUT}")
}

# ------------------------------ Basins (once) ---------------------------------
# Read, normalize around the dateline, and fix topology up front
inland <- readRDS("input/raw/basin/basin_sf_v1.rds") %>%
  st_transform(4326) %>%
  suppressWarnings(st_make_valid()) %>%
  st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE) %>%
  #st_collection_extract("POLYGON", warn = FALSE) %>%
  st_as_sf()

basin_ids <- inland %>% st_drop_geometry() %>% pull(basin_id)

# ------------------------------- Rasters (once) --------------------------------
# List candidate ensemble rasters only once; keep order stable and limit to first 50 if desired
ras_files <- list.files(
  DIR_RAST,
  pattern = "^waterTemp_annual_mean_ensemble.*\\.(tif|tiff)$",
  full.names = TRUE
)

# If you truly want only the first 50, keep the next line; otherwise remove it.
ras_files <- ras_files[seq_len(min(50L, length(ras_files)))]

if (length(ras_files) == 0L) {
  stop("No discharge ensemble rasters found in: ", DIR_RAST)
}

# Read the first raster to capture the common CRS (assumes all ensemble rasters share CRS)
ras0 <- terra::rast(ras_files[1])
ras_crs <- terra::crs(ras0)

# ------------------------------- Helpers --------------------------------------
# Convert an occurrence object to sf; supports either an sf or a data.frame with lon/lat
as_occ_sf <- function(x, lon = "decimalLongitude", lat = "decimalLatitude", crs = 4326) {
  if (inherits(x, "sf")) {
    if (is.na(sf::st_crs(x))) x <- sf::st_set_crs(x, crs)
    return(x)
  }
  stopifnot(all(c(lon, lat) %in% names(x)))
  x %>%
    filter(!is.na(.data[[lon]]), !is.na(.data[[lat]])) %>%
    st_as_sf(coords = c(lon, lat), crs = crs, remove = FALSE)
}

# Build a clipped equal-area grid for a basin (computed once per basin)
build_basin_grid <- function(basin_sf, grid_km = GRID_KM) {
  basin_eq <- st_transform(basin_sf, EQUAL_AREA_CRS)
  cell_m   <- rep(grid_km * 1000, 2)
  st_make_grid(basin_eq, cellsize = cell_m, what = "polygons") %>%
    suppressWarnings(st_make_valid()) %>%
    st_intersection(basin_eq) %>%
    st_as_sf() %>%
    mutate(grid_id = dplyr::row_number(), .before = 1)
}

# Small NA-safe mean
safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

# --------------------------------------------------------------------
# Helper: robust year parser (drop-in replacement)
# --------------------------------------------------------------------
# Reads years from layer names when available; otherwise falls back to
# the filename (works for single-layer annual tiles). Never uses
# terra::sources() so it won't break across terra versions.
layer_years <- function(r, filepath = NULL) {
  nm <- names(r)
  yr <- suppressWarnings(
    as.integer(stringr::str_extract(nm, "(?<!\\d)(18|19|20)\\d{2}(?!\\d)"))
  )
  if (length(yr) == terra::nlyr(r) && all(!is.na(yr))) return(yr)
  
  if (!is.null(filepath)) {
    bn  <- basename(filepath)
    yrf <- suppressWarnings(
      as.integer(stringr::str_extract(bn, "(?<!\\d)(18|19|20)\\d{2}(?!\\d)"))
    )
    if (!is.na(yrf)) return(rep(yrf, terra::nlyr(r)))
  }
  
  rep(NA_integer_, terra::nlyr(r))
}

# --------------------------------------------------------------------
# Main worker (depends on: inland, ras_files, as_occ_sf, build_basin_grid, safe_mean)
# --------------------------------------------------------------------
process_basin <- function(basin_id_value) {
  DIR_SPP <- DIR_SPP
  DIR_OCC <- DIR_OCC
  DIR_OUT <- DIR_OUT
  
  out_csv <- file.path(DIR_OUT, paste0(basin_id_value, ".csv"))
  
  # Skip if already processed
  if (file.exists(out_csv)) {
    cli::cli_alert_info("File already exists: {out_csv} - skipping basin")
    return(invisible(out_csv))
  }
  
  cli::cli_h2("Processing basin: {basin_id_value}")
  
  # ---- species table ----
  spp_path <- file.path(DIR_SPP, paste0(basin_id_value, ".csv"))
  if (!file.exists(spp_path)) {
    cli::cli_alert_warning("Species table not found: {spp_path} - skipping basin")
    return(NULL)
  }
  
  spp <- read.csv(spp_path, check.names = FALSE)
  if (!is.data.frame(spp) || !"valid_name" %in% names(spp)) {
    cli::cli_alert_warning("Invalid species table (missing 'valid_name'): {spp_path} - skipping basin")
    return(NULL)
  }
  
  # Keep an explicit, minimal schema; fallback to first 3 columns if needed
  keep_cols <- intersect(c("basin_id","valid_name", "anything_else"), names(spp))
  if (length(keep_cols) == 0L) keep_cols <- names(spp)[seq_len(min(3L, ncol(spp)))]
  spp <- tibble::as_tibble(spp[, keep_cols, drop = FALSE])
  
  # ---- basin geometry & grid (once per basin) ----
  basin <- inland %>% dplyr::filter(.data$basin_id %in% basin_id_value)
  if (nrow(basin) == 0L) {
    cli::cli_alert_warning("No geometry for basin_id: {basin_id_value} - skipping")
    return(NULL)
  }
  
  grid_eq <- build_basin_grid(basin, grid_km = GRID_KM)  # equal-area grid for joins/counts
  
  results <- list()
  
  # ---- iterate species (ALWAYS emit rows, even if NA) ----
  for (sp in unique(spp$valid_name)) {
    occ_path <- file.path(DIR_OCC, paste0(stringr::str_replace_all(sp, " ", "_"), ".rds"))
    
    # Helper: build a full NA table covering all rasters/layers for this species
    build_na_rows <- function() {
      purrr::map(ras_files, function(fp) {
        r <- terra::rast(fp)
        tibble::tibble(
          valid_name = sp,
          year       = layer_years(r, filepath = fp),
          watertemp  = rep(NA_real_, terra::nlyr(r))
        )
      }) %>% dplyr::bind_rows()
    }
    
    if (!file.exists(occ_path)) {
      cli::cli_alert_warning("Occurrence file missing: {occ_path} - emitting NA rows for {.bold {sp}}")
      sp_rows <- build_na_rows()
    } else {
      cli::cli_alert_info("Reading occurrences for species {.bold {sp}} in basin {.bold {basin_id_value}}")
      occ_raw <- readRDS(occ_path)
      occ_sf  <- tryCatch(as_occ_sf(occ_raw), error = function(e) NULL)
      
      if (is.null(occ_sf) || nrow(occ_sf) == 0L) {
        cli::cli_alert_warning("No valid occurrences for {.bold {sp}} - emitting NA rows")
        sp_rows <- build_na_rows()
      } else {
        # Count occurrences per cell on the equal-area grid
        occ_eq <- sf::st_transform(occ_sf, sf::st_crs(grid_eq))
        n_pts  <- lengths(sf::st_intersects(grid_eq, occ_eq))
        keep_i <- which(n_pts > 0L)
        
        if (length(keep_i) == 0L) {
          cli::cli_alert_info("No grid cells with points for {.bold {sp}} - emitting NA rows")
          sp_rows <- build_na_rows()
        } else {
          # Extract means for each raster file; reproject only the occupied cells per file
          sp_rows <- purrr::map(ras_files, function(fp) {
            r <- terra::rast(fp)
            grid_ras_crs_this <- sf::st_transform(grid_eq[keep_i, , drop = FALSE], terra::crs(r))
            
            vals <- exactextractr::exact_extract(
              r,
              grid_ras_crs_this,
              "mean",
              progress = FALSE
            )
            
            if (length(vals) == 0L) {
              lyr_means <- rep(NA_real_, terra::nlyr(r))
            } else {
              m <- do.call(rbind, lapply(vals, as.numeric))
              lyr_means <- if (is.null(dim(m))) safe_mean(m) else colMeans(m, na.rm = TRUE)
            }
            
            tibble::tibble(
              valid_name = sp,
              year       = layer_years(r, filepath = fp),
              watertemp  = as.numeric(lyr_means)
            )
          }) %>% dplyr::bind_rows()
        }
      }
    }
    
    # Merge species metadata (keeps your additional columns) and append
    merged <- spp %>%
      dplyr::filter(.data$valid_name == sp) %>%
      dplyr::right_join(sp_rows, by = "valid_name")
    
    results <- append(results, list(merged))
    cli::cli_alert_info("Finished {.bold {sp}} in basin {.bold {basin_id_value}} ({nrow(merged)} rows)")
  }
  
  # ---- write output (always write, even if only NA rows) ----
  out_tbl <- dplyr::bind_rows(results)
  if (!dir.exists(DIR_OUT)) dir.create(DIR_OUT, recursive = TRUE)
  readr::write_csv(out_tbl, out_csv)
  cli::cli_alert_success("Saved basin-level discharge table: {out_csv}")
  
  invisible(out_csv)
  gc()
}



# ------------------------------- Run ------------------------------------------
pb_id <- cli::cli_progress_bar("Processing basins", total = length(basin_ids))

for (i in seq_along(basin_ids)) {
  id <- basin_ids[i]
  
  # Run your worker
  try(process_basin(id), silent = TRUE)
  
  # Increment the same progress bar using its ID
  cli::cli_progress_update(id = pb_id, inc = 1)
}

# Close the bar explicitly
cli::cli_progress_done(id = pb_id)


