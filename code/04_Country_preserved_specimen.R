library(furrr)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(sf)
library(ggplot2)
library(cli)

# Suppress all warnings
options(warn = -1)
options(sf_use_s2 = FALSE)
# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
path_species_dir   <- "input/processed/country_sp/"
path_occ_dir       <- "input/processed/occ/"
path_output_dir    <- "input/data_prep/country_preserved_specimen/"

# Ensure output directory exists
if (!dir.exists(path_output_dir)) {
  dir.create(path_output_dir, recursive = TRUE)
  cli_alert_success("Created directory: {path_output_dir}")
}

# ------------------------------------------------------------
# Load and prepare country polygons
# ------------------------------------------------------------
country <- readRDS("input/processed/country_fix.rds") %>%
  st_transform(4326) %>%
  suppressWarnings(st_make_valid()) %>%
  # Wrap around the dateline to avoid invalid geometries across 180°
  st_wrap_dateline(
    options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
    quiet   = TRUE
  ) %>%
  st_collection_extract("POLYGON", warn = FALSE) %>%
  st_as_sf()

# Extract country IDs once
country_ids <- country %>%
  st_drop_geometry() %>%
  pull(iso3)

# Create a global progress bar
cli_progress_bar("Processing country", total = length(country_ids))

# ------------------------------------------------------------
# Helper: does a basisOfRecord entry represent a preserved specimen?
# We accept multiple variants (case differences and formatting differences)
# ------------------------------------------------------------
is_preserved_specimen <- function(x) {
  # normalize to lowercase and remove spaces/underscores
  norm <- x %>%
    tolower() %>%
    gsub("_", "", ., fixed = TRUE) %>%
    gsub(" ", "", ., fixed = TRUE)
  
  norm %in% c(
    "PRESERVED_SPECIMEN",   # standard Darwin Core
    "preservedspecimen",   # underscore form
    "PreservedSpecimen"    # (kept multiple on purpose, harmless)
  )
}

# ------------------------------------------------------------
# Core worker: process a single country
# ------------------------------------------------------------
process_country <- function(country_id_value) {
  
  out_file <- paste0(path_output_dir, country_id_value, ".csv")
  
  # Skip if already processed
  if (file.exists(out_file)) {
    cli::cli_alert_info("File already exists: {out_file} - skipping country")
    return(invisible(out_file))
  }
  
  cli_h2(paste("Processing country:", country_id_value))
  
  # 1. species list ---------------------------------------------------------
  spp_file <- paste0(path_species_dir, country_id_value, ".csv")
  if (!file.exists(spp_file)) {
    cli_alert_warning(paste("Species file not found:", spp_file, "- skipping country"))
    return(invisible(NULL))
  }
  
  spp <- read.csv(spp_file, stringsAsFactors = FALSE)[, c(1, 2)]
  
  if (!is.data.frame(spp)) {
    cli_alert_warning("Species table is not a data frame - skipping this country")
    return(invisible(NULL))
  }
  if (!"valid_name" %in% colnames(spp)) {
    cli_alert_warning("'valid_name' column missing - skipping this country")
    return(invisible(NULL))
  }
  
  all_results <- dplyr::tibble()
  
  # 2. country polygon (once per country) -----------------------------------
  country_sf <- country %>%
    dplyr::filter(.data$iso3 == country_id_value)
  
  if (nrow(country_sf) == 0) {
    cli_alert_warning(paste("No matching country geometry for", country_id_value, "- skipping country"))
    return(invisible(NULL))
  }
  
  # make sure it's valid and in WGS84 lon/lat
  country_sf_4326 <- suppressWarnings(
    country_sf %>%
      sf::st_make_valid() %>%
      sf::st_transform(4326) %>%
      # buffer(0) to heal self-intersections that ruin st_within()
      sf::st_buffer(0)
  )
  
  # precompute bounding box of country in lon/lat for cheap filtering
  bb <- sf::st_bbox(country_sf_4326)
  minx <- bb["xmin"]; maxx <- bb["xmax"]
  miny <- bb["ymin"]; maxy <- bb["ymax"]
  
  # 3. iterate species ------------------------------------------------------
  for (sp in unique(spp$valid_name)) {
    cli_alert_info("Processing species {.bold {sp}} in {.bold {country_id_value}}")
    
    occ_file <- paste0(path_occ_dir, gsub(" ", "_", sp), ".rds")
    
    if (!file.exists(occ_file)) {
      cli_alert_warning(
        paste("Occurrence file not found:", occ_file, "- filling NA for", sp)
      )
      na_row <- tibble::tibble(
        iso3                = country_id_value,
        valid_name          = sp,
        year                = NA_integer_,
        preserved_specimen  = NA_integer_
      )
      all_results <- dplyr::bind_rows(all_results, na_row)
      next
    }
    
    occ_raw <- readRDS(occ_file)
    required_cols <- c("decimalLongitude", "decimalLatitude", "year", "basisOfRecord")
    
    if (!all(required_cols %in% names(occ_raw))) {
      cli_alert_warning(
        paste("Missing required columns in", occ_file, "- filling NA for", sp)
      )
      na_row <- tibble::tibble(
        iso3                = country_id_value,
        valid_name          = sp,
        year                = NA_integer_,
        preserved_specimen  = NA_integer_
      )
      all_results <- dplyr::bind_rows(all_results, na_row)
      next
    }
    
    # 3a. fast bbox pre-filter BEFORE sf conversion (huge memory saver)
    occ_filt <- occ_raw %>%
      dplyr::filter(
        !is.na(decimalLongitude),
        !is.na(decimalLatitude),
        decimalLongitude >= minx,
        decimalLongitude <= maxx,
        decimalLatitude  >= miny,
        decimalLatitude  <= maxy
      )
    
    if (!is.data.frame(occ_filt) || nrow(occ_filt) == 0L) {
      cli_alert_warning(
        paste("No records for", sp, "within bounding box of", country_id_value, "- filling NA")
      )
      na_row <- tibble::tibble(
        iso3                = country_id_value,
        valid_name          = sp,
        year                = NA_integer_,
        preserved_specimen  = NA_integer_
      )
      all_results <- dplyr::bind_rows(all_results, na_row)
      next
    }
    
    # 3b. convert bbox-filtered records to sf points in 4326
    occ_sf <- tryCatch(
      {
        sf::st_as_sf(
          occ_filt,
          coords = c("decimalLongitude", "decimalLatitude"),
          crs = 4326,
          remove = FALSE
        )
      },
      error = function(e) NULL
    )
    
    if (is.null(occ_sf) || nrow(occ_sf) == 0L) {
      cli_alert_warning(
        paste("st_as_sf failed or 0 rows for", sp, "after bbox filter - filling NA")
      )
      na_row <- tibble::tibble(
        iso3                = country_id_value,
        valid_name          = sp,
        year                = NA_integer_,
        preserved_specimen  = NA_integer_
      )
      all_results <- dplyr::bind_rows(all_results, na_row)
      next
    }
    
    # 3c. spatial filter points truly inside polygon using st_within()
    inside_idx <- tryCatch(
      {
        sf::st_within(occ_sf, country_sf_4326, sparse = TRUE)
      },
      error = function(e) NULL
    )
    
    if (is.null(inside_idx)) {
      cli_alert_warning(
        paste("st_within failed for", sp, "in", country_id_value, "- filling NA")
      )
      na_row <- tibble::tibble(
        iso3                = country_id_value,
        valid_name          = sp,
        year                = NA_integer_,
        preserved_specimen  = NA_integer_
      )
      all_results <- dplyr::bind_rows(all_results, na_row)
      next
    }
    
    keep_rows <- lengths(inside_idx) > 0
    occ_in_country <- occ_sf[keep_rows, , drop = FALSE]
    
    if (nrow(occ_in_country) == 0L) {
      cli_alert_warning(
        paste("No in-country points for", sp, "in", country_id_value, "- filling NA")
      )
      na_row <- tibble::tibble(
        iso3                = country_id_value,
        valid_name          = sp,
        year                = NA_integer_,
        preserved_specimen  = NA_integer_
      )
      all_results <- dplyr::bind_rows(all_results, na_row)
      next
    }
    
    # 3d. Keep preserved specimens and summarize by year
    occ_preserved <- occ_in_country %>%
      sf::st_drop_geometry() %>%
      dplyr::filter(!is.na(year)) %>%
      dplyr::filter(is_preserved_specimen(basisOfRecord))
    
    if (nrow(occ_preserved) == 0L) {
      cli_alert_warning(
        paste("No preserved specimens for", sp, "in", country_id_value, "- filling NA")
      )
      na_row <- tibble::tibble(
        iso3                = country_id_value,
        valid_name          = sp,
        year                = NA_integer_,
        preserved_specimen  = NA_integer_
      )
      all_results <- dplyr::bind_rows(all_results, na_row)
      next
    }
    
    summary_tbl <- occ_preserved %>%
      dplyr::group_by(year) %>%
      dplyr::summarise(
        preserved_specimen = dplyr::n(),
        .groups = "drop"
      ) %>%
      dplyr::arrange(year) %>%
      dplyr::mutate(
        valid_name = sp,
        iso3 = country_id_value
      ) %>%
      dplyr::select(
        iso3,
        valid_name,
        year,
        preserved_specimen
      )
    
    all_results <- dplyr::bind_rows(all_results, summary_tbl)
    cli_alert_info("Finished {.bold {sp}} in {.bold {country_id_value}}")
  }
  
  # 4. Write output ---------------------------------------------------------
  # even if all_results is empty tibble(), we still write an empty csv
  write.csv(all_results, file = out_file, row.names = FALSE)
  
  cli_alert_success(
    "All species processed for {.bold {country_id_value}} → saved to {.path {out_file}}"
  )
  
  gc()
  invisible(all_results)
}



# ------------------------------------------------------------
# Run for all country
# ------------------------------------------------------------
purrr::walk(country_ids, process_country)

cli_progress_done()
gc()
