# ========================== Libraries =========================================
library(dplyr)        # data manipulation
library(sf)           # vector data
library(terra)        # raster handling / extraction
library(purrr)        # iteration / mapping
library(tidyr)        # tidy helpers
library(stringr)      # string cleanup
library(cli)          # CLI messages
library(exactextractr) # polygon-weighted raster summaries

# ========================== Config ============================================
EQUAL_AREA_CRS <- "+proj=eck4 +lon_0=0 +lat_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

GRID_KM        <- 50  # 50 km grid cells (equal-area projection)
DIR_SPP        <- "input/processed/country_sp"          # country -> species table
DIR_OCC        <- "input/processed/occ"               # per-species occurrence .rds
DIR_WORLDPOP   <- "input/raw/worldpop"                # population density rasters
DIR_OUT        <- "input/data_prep/country_population_density"

# ensure output dir exists
if (!dir.exists(DIR_OUT)) {
  dir.create(DIR_OUT, recursive = TRUE)
  cli::cli_alert_info("Created output directory: {DIR_OUT}")
} else {
  cli::cli_alert_success("Output directory already exists: {DIR_OUT}")
}

# ========================== country (once) =====================================
# Read country polygons once, normalize around the dateline, and fix topology
country <- readRDS("input/processed/country_fix.rds") %>%
  st_transform(4326) %>%
  suppressWarnings(st_make_valid()) %>%
  st_wrap_dateline(
    options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
    quiet   = TRUE
  ) %>%
  st_as_sf()

country_ids <- country %>%
  st_drop_geometry() %>%
  pull(iso3)

# ========================== Population rasters (once) =========================
# Load all population rasters up front into memory and keep filename metadata
worldpop_files <- list.files(DIR_WORLDPOP, full.names = TRUE)

if (length(worldpop_files) == 0L) {
  stop("No population rasters found in: ", DIR_WORLDPOP)
}

cli::cli_h2("Preloading population rasters into memory...")

worldpop_rasters_all <- lapply(worldpop_files, terra::rast)

cli::cli_alert_success(
  "Loaded {length(worldpop_rasters_all)} population rasters in memory"
)

# Split into "current" vs "past" groups, following original index logic
current_idx <- 1:33
past_idx    <- 34:57

current_files   <- worldpop_files[current_idx]
past_files      <- worldpop_files[past_idx]

# ========================== Helpers ===========================================

# NA-safe mean of a numeric vector
safe_mean <- function(x) {
  if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

# Extract year(s) from raster layer names or fallback to filename.
# (Same logic style as in your discharge code)
layer_years <- function(r, filepath = NULL) {
  nm <- names(r)
  yr <- suppressWarnings(
    as.integer(stringr::str_extract(nm, "(?<!\\d)(17|18|19|20)\\d{2}(?!\\d)"))
  )
  if (length(yr) == terra::nlyr(r) && all(!is.na(yr))) return(yr)
  
  if (!is.null(filepath)) {
    bn  <- basename(filepath)
    yrf <- suppressWarnings(
      as.integer(stringr::str_extract(bn, "(?<!\\d)(17|18|19|20)\\d{2}(?!\\d)"))
    )
    if (!is.na(yrf)) return(rep(yrf, terra::nlyr(r)))
  }
  
  rep(NA_integer_, terra::nlyr(r))
}

# convert raw occurrence data (either sf or data.frame with lon/lat) to sf in WGS84
as_occ_sf <- function(x,
                      lon = "decimalLongitude",
                      lat = "decimalLatitude",
                      crs = 4326) {
  if (inherits(x, "sf")) {
    if (is.na(sf::st_crs(x))) {
      x <- sf::st_set_crs(x, crs)
    }
    return(x)
  }
  
  req_cols <- c(lon, lat)
  if (!all(req_cols %in% names(x))) return(NULL)
  
  x %>%
    dplyr::filter(!is.na(.data[[lon]]), !is.na(.data[[lat]])) %>%
    sf::st_as_sf(coords = c(lon, lat), crs = crs, remove = FALSE)
}


# build an equal-area grid for ONE country, clipped to country extent
build_country_grid <- function(country_sf, grid_km = GRID_KM) {
  country_eq <- sf::st_transform(country_sf, EQUAL_AREA_CRS)
  cell_m   <- rep(grid_km * 1000, 2)
  country_eq <- country_eq %>%
    st_make_valid() %>%
    st_buffer(0) %>%                       
    st_collection_extract("POLYGON") %>%   
    st_set_precision(100) %>%              
    st_make_valid()
  country_eq <- st_union(country_eq) %>% st_as_sf()
  sf::st_make_grid(
    country_eq,
    cellsize = cell_m,
    what     = "polygons"
  ) %>%
    suppressWarnings(sf::st_make_valid()) %>%
    sf::st_intersection(country_eq) %>%
    sf::st_as_sf() %>%
    dplyr::mutate(grid_id = dplyr::row_number(), .before = 1)
}

# Summarize population density for ONE species in ONE country across occupied cells,
# using the discharge-style exactextractr workflow.
#
# Args:
#   ras_files   : character vector of raster file paths (e.g. current_files or past_files)
#   grid_eq     : sf grid (equal-area cells) for this country
#   keep_i      : integer indices of occupied cells in that grid
#   sp          : species name (string)
#
# Returns:
#   tibble(valid_name, year, population_density)
#
summarize_population_for_species <- function(ras_files, grid_eq, keep_i, sp) {
  # helper to emit NA rows if no occupied cells
  build_na_rows <- function() {
    purrr::map(ras_files, function(fp) {
      r <- terra::rast(fp)
      tibble::tibble(
        valid_name         = sp,
        year               = layer_years(r, filepath = fp),
        population_density = rep(NA_real_, terra::nlyr(r))
      )
    }) %>% dplyr::bind_rows()
  }
  
  if (length(keep_i) == 0L) {
    # species has no occupied cells in this country
    return(build_na_rows())
  }
  
  occ_cells_eq <- grid_eq[keep_i, , drop = FALSE]
  
  sp_rows <- purrr::map(ras_files, function(fp) {
    r <- terra::rast(fp)
    
    # Reproject occupied cells to the raster CRS
    cells_in_r_crs <- sf::st_transform(occ_cells_eq, terra::crs(r))
    
    # exact_extract returns a list: one entry per polygon
    vals_list <- exactextractr::exact_extract(
      r,
      cells_in_r_crs,
      "mean",
      progress = FALSE
    )
    
    # Convert list -> matrix (rows = polygons, cols = raster layers)
    if (length(vals_list) == 0L) {
      layer_means <- rep(NA_real_, terra::nlyr(r))
    } else {
      m <- do.call(rbind, lapply(vals_list, as.numeric))
      # single-cell edge case
      layer_means <- if (is.null(dim(m))) {
        safe_mean(m)
      } else {
        colMeans(m, na.rm = TRUE)
      }
    }
    
    tibble::tibble(
      valid_name         = sp,
      year               = layer_years(r, filepath = fp),
      population_density = as.numeric(layer_means)
    )
  }) %>%
    dplyr::bind_rows()
  
  sp_rows
}

# ========================== Worker ============================================
# For one country:
#   1. read species list
#   2. build equal-area grid
#   3. for each species:
#        - read occurrences
#        - find occupied grid cells
#        - summarize population density (current + past rasters)
#   4. write CSV
process_country <- function(country_id_value) {
  out_csv <- file.path(DIR_OUT, paste0(country_id_value, ".csv"))
  
  # optional skip-on-existing behavior
  if (file.exists(out_csv)) {
    cli::cli_alert_info("File already exists: {out_csv} - skipping country")
    return(invisible(out_csv))
  }
  
  cli::cli_h2("Processing country: {country_id_value}")
  
  # ---- species table ----
  spp_path <- file.path(DIR_SPP, paste0(country_id_value, ".csv"))
  if (!file.exists(spp_path)) {
    cli::cli_alert_warning("Species table not found: {spp_path} - skipping country")
    return(NULL)
  }
  
  spp_raw <- read.csv(spp_path, stringsAsFactors = FALSE)
  if (!is.data.frame(spp_raw) || nrow(spp_raw) == 0L) {
    cli::cli_alert_warning("Species table empty/invalid for country {country_id_value} - skipping")
    return(NULL)
  }
  if (!"valid_name" %in% names(spp_raw)) {
    cli::cli_alert_warning("'valid_name' missing in {spp_path} - skipping country {country_id_value}")
    return(NULL)
  }
  
  spp_tbl <- spp_raw[, 1:2, drop = FALSE]
  
  # ---- country geometry ----
  country_sf <- country %>% dplyr::filter(.data$iso3 %in% country_id_value)
  if (nrow(country_sf) == 0L) {
    cli::cli_alert_warning("No geometry for country {country_id_value} - skipping")
    return(NULL)
  }
  
  # ---- build equal-area grid for this country (robust) ----
  grid_eq <- tryCatch(
    {
      build_country_grid(country_sf, grid_km = GRID_KM)
    },
    error = function(e) {
      cli::cli_alert_warning(
        "build_country_grid() failed for {.bold {country_id_value}}: {conditionMessage(e)}"
      )
      NULL
    }
  )
  
  # Helper: build NA table for one species if grid is unusable
  make_na_row_for_species <- function(sp_name) {
    tibble::tibble(
      iso3           = country_id_value,
      valid_name         = sp_name,
      year               = NA_character_,
      population_density = NA_real_
    )
  }
  
  # If we totally failed to build a grid (grid_eq NULL or 0 rows),
  # we still write an output for every species but with NA.
  if (is.null(grid_eq) || nrow(grid_eq) == 0L) {
    cli::cli_alert_warning(
      "No usable grid for {.bold {country_id_value}} (NULL or 0 rows). Writing NA densities."
    )
    
    out_tbl <- purrr::map_dfr(
      unique(spp_tbl$valid_name),
      make_na_row_for_species
    )
    
    # final cleanup: year format and column order
    out_tbl <- out_tbl %>%
      dplyr::mutate(
        year = stringr::str_extract(.data$year, "\\d{4}(?:AD)?"),
        year = as.character(year)
      ) %>%
      dplyr::select(
        iso3,
        valid_name,
        year,
        population_density
      )
    
    # write anyway
    write.csv(out_tbl, file = out_csv, row.names = FALSE)
    cli::cli_alert_success("Saved country-level population table: {out_csv}")
    
    return(invisible(out_csv))
  }
  
  # ---- iterate species ----
  results <- list()
  
  for (sp in unique(spp_tbl$valid_name)) {
    cli::cli_alert_info(
      "Processing species {.bold {sp}} in country {.bold {country_id_value}}"
    )
    
    # read occurrence
    occ_path <- file.path(
      DIR_OCC,
      paste0(stringr::str_replace_all(sp, " ", "_"), ".rds")
    )
    
    occ_sf <- if (file.exists(occ_path)) {
      occ_raw <- readRDS(occ_path)
      as_occ_sf(occ_raw)
    } else {
      NULL
    }
    
    # occupied grid cells
    if (is.null(occ_sf) || nrow(occ_sf) == 0L) {
      cli::cli_alert_warning(
        "No valid occurrences for {.bold {sp}} in country {.bold {country_id_value}}"
      )
      keep_i <- integer(0)
    } else {
      occ_eq <- sf::st_transform(occ_sf, sf::st_crs(grid_eq))
      n_pts  <- lengths(sf::st_intersects(grid_eq, occ_eq))
      keep_i <- which(n_pts > 0L)
    }
    
    # summarize pop density for current & past rasters
    dens_current <- summarize_population_for_species(
      ras_files = current_files,
      grid_eq   = grid_eq,
      keep_i    = keep_i,
      sp        = sp
    )
    dens_past <- summarize_population_for_species(
      ras_files = past_files,
      grid_eq   = grid_eq,
      keep_i    = keep_i,
      sp        = sp
    )
    
    dens_all <- dplyr::bind_rows(dens_current, dens_past)
    
    if (nrow(dens_all) == 0L) {
      merged_sp <- make_na_row_for_species(sp)
    } else {
      merged_sp <- dens_all %>%
        dplyr::mutate(
          iso3 = country_id_value
        ) %>%
        dplyr::select(
          iso3,
          valid_name,
          year,
          population_density
        )
    }
    
    results <- append(results, list(merged_sp))
    
    cli::cli_alert_info(
      "Finished {.bold {sp}} in country {.bold {country_id_value}} ({nrow(merged_sp)} rows)"
    )
  }
  
  # ---- combine, clean year, enforce columns ----
  out_tbl <- dplyr::bind_rows(results) %>%
    dplyr::mutate(
      year = stringr::str_extract(.data$year, "\\d{4}(?:AD)?"),
      year = as.character(year)
    ) %>%
    dplyr::select(
      iso3,
      valid_name,
      year,
      population_density
    )
  
  # safety fallback in case results was length 0
  if (nrow(out_tbl) == 0L) {
    out_tbl <- tibble::tibble(
      iso3           = country_id_value,
      valid_name         = NA_character_,
      year               = NA_character_,
      population_density = NA_real_
    )
  }
  
  # ---- write ----
  write.csv(out_tbl, file = out_csv, row.names = FALSE)
  cli::cli_alert_success("Saved country-level population table: {out_csv}")
  
  invisible(out_csv)
}



# ========================== Run ===============================================
pb_id <- cli::cli_progress_bar("Processing country", total = length(country_ids))

for (i in seq_along(country_ids)) { #
  id_now <- country_ids[i]
  
  # guard the loop so one bad basin doesn't kill the run
  try(process_country(id_now), silent = TRUE)
  
  cli::cli_progress_update(id = pb_id, inc = 1)
}

cli::cli_progress_done(id = pb_id)