library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(sf)
library(cli)

# Suppress all warnings
options(warn = -1)
options(sf_use_s2 = FALSE)
# ------------------------------------------------------------
# File paths
# ------------------------------------------------------------
path_species_dir  <- "input/processed/country_sp/"
path_occ_dir      <- "input/processed/occ/"
path_output_dir   <- "input/data_prep/country_sampling_effort/"

# Ensure output directory exists
if (!dir.exists(path_output_dir)) {
  dir.create(path_output_dir, recursive = TRUE)
  cli_alert_success("Created directory: {path_output_dir}")
}

# ------------------------------------------------------------
# Load and fix country geometry
# ------------------------------------------------------------
country <- readRDS("input/processed/country_fix.rds") %>%
  st_transform(4326) %>%
  suppressWarnings(st_make_valid()) %>%
  # Wrap around the dateline to fix global geometries
  st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"), quiet = TRUE) %>%
  st_collection_extract("POLYGON", warn = FALSE) %>%
  st_as_sf()

# Get all country IDs
country_ids <- country %>%
  st_drop_geometry() %>%
  pull(iso3)

# Create progress bar
cli_progress_bar("Processing country", total = length(country_ids))

# ------------------------------------------------------------
# Function: process a single country
# ------------------------------------------------------------
process_country <- function(country_id_value) {
  cli_h2(paste("Processing country:", country_id_value))
  
  # ---------------------------------------------------------------------------
  # helper: safely intersect grid with country polygon
  # returns an sf with at least columns (grid_id, geometry) or NULL
  # ---------------------------------------------------------------------------
  safe_intersection <- function(grid_sf, country_poly) {
    
    # First attempt: real geometric clip (nice coastal trimming)
    out_try <- tryCatch({
      st_intersection(grid_sf, country_poly) %>%
        st_as_sf()
    }, error = function(e) e)
    
    if (!inherits(out_try, "error")) {
      # Success: make sure geometry is valid
      out_ok <- suppressWarnings(st_make_valid(out_try))
      return(out_ok)
    }
    
    # If we get here, st_intersection failed, so we fall back:
    cli_alert_warning(paste(
      "st_intersection() failed for", country_id_value, "→",
      conditionMessage(out_try),
      "Falling back to spatial filter (st_intersects)."
    ))
    
    # 1. Find which grid cells even touch the country polygon
    suppressWarnings({
      hits <- st_intersects(grid_sf, country_poly)
    })
    keep_idx <- which(lengths(hits) > 0)
    
    if (length(keep_idx) == 0) {
      return(NULL)
    }
    
    # 2. Subset those cells BUT keep sf-ness and CRS explicitly
    #    grid_sf[keep_idx, , drop = FALSE] *should* stay sf,
    #    but we've seen edge cases where it degrades,
    #    so we'll rebuild explicitly just to be safe.
    subset_cells <- grid_sf[keep_idx, , drop = FALSE]
    
    # enforce sf class + geometry column + CRS
    subset_cells <- st_as_sf(
      subset_cells,
      crs = st_crs(grid_sf)
    )
    
    # final validity cleanup
    subset_cells <- suppressWarnings(st_make_valid(subset_cells))
    subset_cells
  }
  
  # ---------------------------------------------------------------------------
  # 1. Load species list for this country
  # ---------------------------------------------------------------------------
  spp_file <- paste0(path_species_dir, country_id_value, ".csv")
  if (!file.exists(spp_file)) {
    cli_alert_warning(paste("File not found:", spp_file, "- skipping"))
    return(invisible(NULL))
  }
  
  spp <- read.csv(spp_file)[, c(1, 2)]
  
  if (!is.data.frame(spp)) {
    cli_alert_warning("Invalid species table - skipping")
    return(invisible(NULL))
  }
  if (!"valid_name" %in% colnames(spp)) {
    cli_alert_warning("'valid_name' column missing - skipping")
    return(invisible(NULL))
  }
  
  # ---------------------------------------------------------------------------
  # 2. Get country polygon and clean aggressively
  # ---------------------------------------------------------------------------
  country_sf <- country %>%
    dplyr::filter(iso3 %in% country_id_value)
  
  if (nrow(country_sf) == 0) {
    cli_alert_warning(paste("No geometry for country:", country_id_value, "- skipping"))
    return(invisible(NULL))
  }
  
  # basic clean in lon/lat
  country_sf <- country_sf %>%
    st_transform(4326) %>%
    suppressWarnings(st_make_valid()) %>%
    suppressWarnings(st_buffer(0))
  
  # project to equal-area CRS (Eckert IV)
  ECK4 <- "+proj=eck4 +lon_0=0 +lat_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  country_proj <- country_sf %>%
    st_transform(ECK4) %>%
    suppressWarnings(st_make_valid()) %>%
    suppressWarnings(st_buffer(0))
  
  # dissolve multiparts/islands
  country_proj_union <- suppressWarnings(st_union(country_proj))
  
  # ---------------------------------------------------------------------------
  # 3. Build 50 km grid in projected space and clip it safely
  # ---------------------------------------------------------------------------
  grid_raw <- st_make_grid(
    country_proj_union,
    cellsize = c(50000, 50000),
    what = "polygons"
  )
  
  grid_sf <- st_as_sf(grid_raw) %>%
    dplyr::mutate(grid_id = dplyr::row_number()) %>%
    suppressWarnings(st_make_valid(.))
  
  fishnet_clipped_sf <- safe_intersection(grid_sf, country_proj_union)
  
  # If still NULL / empty, we cannot spatially bin effort in this country.
  # We still output NA rows for each species so downstream code is consistent.
  if (is.null(fishnet_clipped_sf) || nrow(fishnet_clipped_sf) == 0) {
    cli_alert_warning(paste(
      "No usable fishnet for country", country_id_value,
      "→ exporting NA rows only."
    ))
    
    na_rows <- tibble::tibble(
      iso3             = country_id_value,
      valid_name       = unique(spp$valid_name),
      sampling_year    = NA_integer_,
      sampling_effort  = NA_real_
    )
    
    output_file <- paste0(path_output_dir, country_id_value, ".csv")
    write.csv(na_rows, file = output_file, row.names = FALSE)
    
    cli_alert_success(
      "Finished country {.bold {country_id_value}} (geometry fallback only) → saved {.path {output_file}}"
    )
    
    gc()
    return(invisible(na_rows))
  }
  
  # NOTE: we no longer do select(grid_id, geometry) here.
  # fishnet_clipped_sf should ALREADY have grid_id + geometry and be sf.
  # So we trust it as-is.
  
  # ---------------------------------------------------------------------------
  # 4. Loop over each species and compute sampling effort
  # ---------------------------------------------------------------------------
  all_results <- dplyr::tibble()
  
  for (species in unique(spp$valid_name)) {
    
    distribution_file <- paste0(
      path_occ_dir,
      gsub(" ", "_", species),
      ".rds"
    )
    
    if (!file.exists(distribution_file)) {
      cli_alert_warning(paste("Occurrence file not found:", distribution_file, "- filling NA"))
      empty_row <- tibble::tibble(
        iso3             = country_id_value,
        valid_name       = species,
        sampling_year    = NA_integer_,
        sampling_effort  = NA_real_
      )
      all_results <- dplyr::bind_rows(all_results, empty_row)
      next
    }
    
    cli_alert_info("Reading occurrence data for species: {.bold {species}}")
    
    distribution_raw <- readRDS(distribution_file)
    
    # Check required fields
    if (!all(c("decimalLongitude", "decimalLatitude", "year") %in% names(distribution_raw))) {
      cli_alert_warning(paste("Missing coords/year in", distribution_file, "- filling NA"))
      empty_row <- tibble::tibble(
        iso3             = country_id_value,
        valid_name       = species,
        sampling_year    = NA_integer_,
        sampling_effort  = NA_real_
      )
      all_results <- dplyr::bind_rows(all_results, empty_row)
      next
    }
    
    # Build sf points
    distribution_sf <- distribution_raw %>%
      dplyr::filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
      st_as_sf(
        coords  = c("decimalLongitude", "decimalLatitude"),
        crs     = 4326,
        remove  = FALSE
      )
    
    if (nrow(distribution_sf) == 0) {
      cli_alert_warning(paste("No valid points for", species, "- filling NA"))
      empty_row <- tibble::tibble(
        iso3             = country_id_value,
        valid_name       = species,
        sampling_year    = NA_integer_,
        sampling_effort  = NA_real_
      )
      all_results <- dplyr::bind_rows(all_results, empty_row)
      next
    }
    
    # Reproject to match fishnet CRS
    points_proj <- st_transform(distribution_sf, st_crs(fishnet_clipped_sf))
    
    # Try fine intersection (points in cells). Fallback: st_join.
    points_in_grid <- tryCatch({
      suppressWarnings(st_intersection(fishnet_clipped_sf, points_proj))
    }, error = function(e) {
      cli_alert_warning(paste(
        "st_intersection() failed for species", species,
        "in", country_id_value, "→", conditionMessage(e),
        "Falling back to st_join()."
      ))
      suppressWarnings(st_join(fishnet_clipped_sf, points_proj, left = FALSE))
    })
    
    if (is.null(points_in_grid) || nrow(points_in_grid) == 0) {
      cli_alert_warning(paste("No grid overlap for", species, "in", country_id_value, "- filling NA"))
      empty_row <- tibble::tibble(
        iso3             = country_id_value,
        valid_name       = species,
        sampling_year    = NA_integer_,
        sampling_effort  = NA_real_
      )
      all_results <- dplyr::bind_rows(all_results, empty_row)
      next
    }
    
    # Aggregate sampling effort by year
    grouped_data <- points_in_grid %>%
      st_drop_geometry() %>%
      dplyr::filter(!is.na(year)) %>%
      dplyr::group_by(year, grid_id) %>%
      dplyr::summarise(
        count = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::select(year, count) %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(
        sampling_effort = mean(count, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(valid_name = species)
    
    if (nrow(grouped_data) > 0) {
      grouped_data <- grouped_data %>%
        dplyr::rename(sampling_year = year) %>%
        dplyr::mutate(iso3 = country_id_value) %>%
        dplyr::select(iso3, valid_name, sampling_year, sampling_effort)
      
      all_results <- dplyr::bind_rows(all_results, grouped_data)
      cli_alert_info("Done: {.bold {species}} in {.bold {country_id_value}}")
    } else {
      cli_alert_warning(paste("No usable year info for", species, "- filling NA"))
      empty_row <- tibble::tibble(
        iso3             = country_id_value,
        valid_name       = species,
        sampling_year    = NA_integer_,
        sampling_effort  = NA_real_
      )
      all_results <- dplyr::bind_rows(all_results, empty_row)
    }
  } # end species loop
  
  # ---------------------------------------------------------------------------
  # 5. Save output for this country
  # ---------------------------------------------------------------------------
  output_file <- paste0(path_output_dir, country_id_value, ".csv")
  write.csv(all_results, file = output_file, row.names = FALSE)
  
  cli_alert_success(
    "Finished country {.bold {country_id_value}} → saved {.path {output_file}}"
  )
  
  gc()
  invisible(all_results)
}



# ------------------------------------------------------------
# Run all country
# ------------------------------------------------------------
purrr::walk(country_ids, process_country)

cli_progress_done()
gc()
