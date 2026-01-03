library(furrr)
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
# Paths
# ------------------------------------------------------------
path_species_dir   <- "input/processed/country_sp/"          # per-country species lists
path_occ_dir       <- "input/processed/occ/"               # per-species occurrence .rds
path_output_dir    <- "input/data_prep/country_rarity/"      # output

# Make sure the output directory exists
if (!dir.exists(path_output_dir)) {
  dir.create(path_output_dir, recursive = TRUE)
  cli_alert_success("Created directory: {path_output_dir}")
}

# ------------------------------------------------------------
# Load and prepare country polygons
# IMPORTANT: country MUST already have area_km2 (country area in km²)
# ------------------------------------------------------------
country <- readRDS("input/processed/country_fix.rds") %>%
  st_transform(4326) %>%
  suppressWarnings(st_make_valid()) %>%
  # Wrap across the antimeridian to avoid invalid polygons
  st_wrap_dateline(
    options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
    quiet   = TRUE
  ) %>%
  #st_collection_extract("POLYGON", warn = FALSE) %>%
  st_as_sf()

# NOTE: if 'country' does NOT yet have area_km2, you need to precompute it once
# using an equal-area projection (Eckert IV) and left_join back into 'country'
# before running the rest of this script.

# Get all country IDs
country_ids <- country %>%
  st_drop_geometry() %>%
  pull(iso3)

# ------------------------------------------------------------
# Function: process a single country
# ------------------------------------------------------------

process_country <- function(country_id_value) { 
  out_file <- paste0(path_output_dir, country_id_value, ".csv")
  
  # Skip if already processed
  if (file.exists(out_file)) {
    cli::cli_alert_info("File already exists: {out_file} - skipping country")
    return(invisible(out_file))
  }
  cli::cli_h2(paste("Processing country:", country_id_value))
  # --- helper crs / constants ---------------------------------------------
  ECKERT4_CRS <- "+proj=eck4 +lon_0=0 +lat_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  CELL_M      <- 50000    # 50 km
  CELL_KM2    <- 2500     # 50 km * 50 km
  MIN_AREA_M2 <- 1        # drop garbage slivers
  
  # ------------------------------------------------------------------------
  # 1. species list for this country
  # ------------------------------------------------------------------------
  spp_file <- paste0(path_species_dir, country_id_value, ".csv")
  if (!file.exists(spp_file)) {
    cli_alert_warning(paste("Species file not found:", spp_file, "- skipping country"))
    return(NULL)
  }
  
  spp <- read.csv(spp_file, check.names = FALSE)[, c(1,2)]
  if (!"valid_name" %in% colnames(spp)) {
    cli_alert_warning(paste("No 'valid_name' col in", spp_file, "- skipping country"))
    return(NULL)
  }
  
  # ------------------------------------------------------------------------
  # 2. country polygon + area
  # ------------------------------------------------------------------------
  country_row <- country %>%
    dplyr::filter(iso3 == country_id_value)
  
  if (nrow(country_row) == 0) {
    cli_alert_warning(paste("country polygon not found for", country_id_value, "- skipping"))
    return(NULL)
  }
  
  country_area_km2 <- readRDS("input/data_prep/country_area.rds") %>%
    dplyr::filter(iso3 == country_id_value) %>%
    dplyr::pull(area_km2) %>%
    as.numeric()
  
  if (length(country_area_km2) == 0 || is.na(country_area_km2)) {
    cli_alert_danger(
      paste("No area_km2 found for", country_id_value,
            "in country_area.rds. Aborting for this country.")
    )
    return(NULL)
  }
  
  # ------------------------------------------------------------------------
  # 2a. Clean geometry in lon/lat BEFORE projecting.
  #     We do this step-by-step with guards because st_wrap_dateline()
  #     can return NULL or empty geometry for multi-island countries.
  # ------------------------------------------------------------------------
  
  # 2a-1. force to lon/lat
  country_lonlat <- tryCatch(
    { suppressWarnings(st_transform(country_row, 4326)) },
    error = function(e) NULL
  )
  if (is.null(country_lonlat) || nrow(country_lonlat) == 0) {
    cli_alert_danger(paste("st_transform(4326) failed for", country_id_value, "- skipping"))
    return(NULL)
  }
  
  # 2a-2. wrap at the antimeridian to split polygons that cross 180°
  country_wrapped <- tryCatch(
    {
      suppressWarnings(
        st_wrap_dateline(
          country_lonlat,
          options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
          quiet   = TRUE
        )
      )
    },
    error = function(e) NULL
  )
  # st_wrap_dateline() can return NULL; if so, fall back to original lon/lat
  if (is.null(country_wrapped)) {
    country_wrapped <- country_lonlat
  }
  
  # 2a-3. keep only polygonal parts (drop GEOMETRYCOLLECTION stuff)
  country_polys <- tryCatch(
    {
      suppressWarnings(st_collection_extract(country_wrapped, "POLYGON", warn = FALSE))
    },
    error = function(e) NULL
  )
  # If extraction failed or produced nothing, again fall back
  if (is.null(country_polys)) {
    country_polys <- country_wrapped
  }
  if (is.null(country_polys) || nrow(country_polys) == 0) {
    cli_alert_danger(paste("No polygonal geometry after wrap for", country_id_value, "- skipping"))
    return(NULL)
  }
  
  # 2a-4. heal invalid rings / self-intersections
  country_valid <- tryCatch(
    {
      country_polys %>%
        suppressWarnings(st_make_valid()) %>%
        # buffer(0) is a common trick to fix bow-tie polygons
        suppressWarnings(st_buffer(0))
    },
    error = function(e) NULL
  )
  if (is.null(country_valid) || nrow(country_valid) == 0) {
    cli_alert_danger(paste("Geometry repair failed for", country_id_value, "- skipping"))
    return(NULL)
  }
  
  # ------------------------------------------------------------------------
  # 2b. project to equal-area CRS and drop microscopic slivers
  # ------------------------------------------------------------------------
  country_proj <- tryCatch(
    {
      suppressWarnings(st_transform(country_valid, ECKERT4_CRS)) %>%
        suppressWarnings(st_make_valid()) %>%
        suppressWarnings(st_buffer(0))
    },
    error = function(e) NULL
  )
  if (is.null(country_proj) || nrow(country_proj) == 0) {
    cli_alert_danger(paste("Projection to equal-area CRS failed for", country_id_value, "- skipping"))
    return(NULL)
  }
  
  # remove polygons with near-zero area (can appear after wrapping at 180°)
  poly_area_m2 <- suppressWarnings(sf::st_area(country_proj))
  keep_idx <- as.numeric(poly_area_m2) > MIN_AREA_M2
  country_proj <- country_proj[keep_idx, , drop = FALSE]
  
  if (nrow(country_proj) == 0) {
    cli_alert_danger(
      paste("All polygons degenerate (<", MIN_AREA_M2, "m²) for", country_id_value, "- skipping")
    )
    return(NULL)
  }
  
  # ------------------------------------------------------------------------
  # 3. build ~50 km grid and clip to the country polygon
  # ------------------------------------------------------------------------
  grid_raw <- suppressWarnings(
    st_make_grid(
      country_proj,
      cellsize = c(CELL_M, CELL_M),
      what     = "polygons"
    )
  )
  
  grid_sf <- st_as_sf(grid_raw) %>%
    suppressWarnings(st_make_valid()) %>%
    suppressWarnings(st_buffer(0))
  
  country_proj2 <- country_proj %>%
    st_make_valid() %>%
    st_buffer(0) %>%                       
    st_collection_extract("POLYGON") %>%   
    st_set_precision(100) %>%              
    st_make_valid()
  country_proj2 <- st_union(country_proj2) %>% st_as_sf()
  fishnet_clipped_sf <- tryCatch(
    {
      suppressWarnings(st_intersection(grid_sf, country_proj2)) %>%
        st_as_sf()
    },
    error = function(e) {
      cli_alert_warning(paste(
        "Grid/country intersection failed for", country_id_value,
        ":", conditionMessage(e),
        "→ using empty grid."
      ))
      empty_sf <- grid_sf[0, , drop = FALSE]
      empty_sf$grid_id <- integer(0)
      empty_sf
    }
  )
  
  if (nrow(fishnet_clipped_sf) > 0) {
    fishnet_clipped_sf <- fishnet_clipped_sf %>%
      dplyr::mutate(grid_id = dplyr::row_number())
  } else {
    fishnet_clipped_sf$grid_id <- integer(0)
  }
  
  # ------------------------------------------------------------------------
  # 4. iterate over species
  # ------------------------------------------------------------------------
  all_results <- data.frame(
    iso3             = character(),
    valid_name       = character(),
    sampling_year    = numeric(),
    relative_rarity  = numeric(),
    stringsAsFactors = FALSE
  )
  
  add_na_row_for_sp <- function(sp_name) {
    data.frame(
      iso3             = country_id_value,
      valid_name       = sp_name,
      sampling_year    = NA_real_,
      relative_rarity  = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  
  for (sp in unique(spp$valid_name)) {
    
    occ_file <- paste0(path_occ_dir, gsub(" ", "_", sp), ".rds")
    
    if (!file.exists(occ_file)) {
      cli_alert_warning(
        paste("Occurrence file not found:", occ_file, "- skipping species", sp)
      )
      all_results <- dplyr::bind_rows(all_results, add_na_row_for_sp(sp))
      next
    }
    
    occ_raw <- readRDS(occ_file)
    required_cols <- c("decimalLongitude", "decimalLatitude", "year")
    if (!all(required_cols %in% names(occ_raw))) {
      cli_alert_warning(
        paste("Missing required columns in", occ_file, "for", sp,
              "- producing NA output row.")
      )
      all_results <- dplyr::bind_rows(all_results, add_na_row_for_sp(sp))
      next
    }
    
    occ_sf <- st_as_sf(
      occ_raw,
      coords = c("decimalLongitude", "decimalLatitude"),
      crs = 4326
    )
    occ_proj <- suppressWarnings(st_transform(occ_sf, st_crs(country_proj)))
    
    if (nrow(fishnet_clipped_sf) == 0) {
      all_results <- dplyr::bind_rows(all_results, add_na_row_for_sp(sp))
      next
    }
    
    occ_by_cell_year <- tryCatch(
      {
        suppressWarnings(st_intersection(fishnet_clipped_sf, occ_proj))
      },
      error = function(e) {
        cli_alert_warning(paste(
          "Point/grid intersection failed for", sp, "in", country_id_value,
          ":", conditionMessage(e),
          "→ NA row."
        ))
        NULL
      }
    )
    
    if (is.null(occ_by_cell_year) || nrow(occ_by_cell_year) == 0) {
      all_results <- dplyr::bind_rows(all_results, add_na_row_for_sp(sp))
      next
    }
    
    occ_by_cell_year <- occ_by_cell_year %>%
      st_drop_geometry() %>%
      dplyr::filter(!is.na(year)) %>%
      dplyr::distinct(year, grid_id) %>%
      dplyr::arrange(year, grid_id)
    
    if (nrow(occ_by_cell_year) == 0) {
      all_results <- dplyr::bind_rows(all_results, add_na_row_for_sp(sp))
      next
    }
    
    yr_list <- sort(unique(occ_by_cell_year$year))
    
    cum_df <- purrr::map_dfr(
      yr_list,
      function(y) {
        n_unique_cells_up_to_y <- occ_by_cell_year %>%
          dplyr::filter(year <= y) %>%
          dplyr::pull(grid_id) %>%
          unique() %>%
          length()
        data.frame(
          sampling_year  = y,
          cum_range_size = n_unique_cells_up_to_y,
          stringsAsFactors = FALSE
        )
      }
    )
    
    cum_df <- cum_df %>%
      group_by(sampling_year) %>%
      dplyr::mutate(
        iso3       = country_id_value,
        valid_name = sp,
        occupancy_fraction = ifelse(
          country_area_km2 > 0,
          pmin((cum_range_size * CELL_KM2/ country_area_km2), 1),
          NA_real_
        ),
        relative_rarity = ifelse(
          !is.na(occupancy_fraction) & occupancy_fraction > 0,
          log10((1 / occupancy_fraction) + 1),
          NA_real_
        )
      ) %>%
      dplyr::select(
        iso3,
        valid_name,
        sampling_year,
        relative_rarity
      )
    
    if (nrow(cum_df) == 0) {
      cum_df <- add_na_row_for_sp(sp)
    }
    
    all_results <- dplyr::bind_rows(all_results, cum_df)
  }
  
  # ------------------------------------------------------------------------
  # 5. fallback if nothing was produced
  # ------------------------------------------------------------------------
  if (nrow(all_results) == 0) {
    all_results <- data.frame(
      iso3            = country_id_value,
      valid_name      = NA_character_,
      sampling_year   = NA_real_,
      relative_rarity = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  
  all_results <- all_results %>%
    dplyr::select(iso3, valid_name, sampling_year, relative_rarity)
  
  # ------------------------------------------------------------------------
  # 6. write
  # ------------------------------------------------------------------------
  
  write.csv(
    all_results,
    file = out_file,
    row.names = FALSE
  )
  
  cli::cli_alert_success(paste(
    "Finished country", country_id_value, "→ saved", out_file
  ))
  
  invisible(all_results)
}


# ------------------------------------------------------------
# Run all country
# ------------------------------------------------------------
purrr::walk(country_ids, process_country)

gc()
cli_progress_done()

