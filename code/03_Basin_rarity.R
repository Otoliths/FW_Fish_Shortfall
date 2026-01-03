library(furrr)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(sf)
library(cli)

# Suppress all warnings
options(warn = -1)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
path_species_dir   <- "input/processed/basin_sp/"          # per-basin species lists
path_occ_dir       <- "input/processed/occ/"               # per-species occurrence .rds
path_output_dir    <- "input/data_prep/basin_rarity/"      # output

# Make sure the output directory exists
if (!dir.exists(path_output_dir)) {
  dir.create(path_output_dir, recursive = TRUE)
  cli_alert_success("Created directory: {path_output_dir}")
}

# ------------------------------------------------------------
# Load and prepare basin polygons
# IMPORTANT: inland MUST already have area_km2 (basin area in km²)
# ------------------------------------------------------------
inland <- readRDS("input/raw/basin/basin_sf_v1.rds") %>%
  st_transform(4326) %>%
  suppressWarnings(st_make_valid()) %>%
  # Wrap across the antimeridian to avoid invalid polygons
  st_wrap_dateline(
    options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
    quiet   = TRUE
  ) %>%
  st_collection_extract("POLYGON", warn = FALSE) %>%
  st_as_sf()

# NOTE: if 'inland' does NOT yet have area_km2, you need to precompute it once
# using an equal-area projection (Eckert IV) and left_join back into 'inland'
# before running the rest of this script.

# Get all basin IDs
basin_ids <- inland %>%
  st_drop_geometry() %>%
  pull(basin_id)


# ------------------------------------------------------------
# Function: process a single basin
# ------------------------------------------------------------
process_basin <- function(basin_id_value) { 
  
  cli::cli_h2(paste("Processing basin:", basin_id_value))
  
  out_file <- paste0(path_output_dir, basin_id_value, ".csv")
  
  # Skip if already processed
  if (file.exists(out_file)) {
    cli::cli_alert_info("File already exists: {out_file} - skipping country")
    return(invisible(out_file))
  }
  
  # 0. species list for this basin -------------------------------------------
  spp_file <- paste0(path_species_dir, basin_id_value, ".csv")
  if (!file.exists(spp_file)) {
    cli_alert_warning(paste("Species file not found:", spp_file, "- skipping basin"))
    return(NULL)
  }
  
  spp <- read.csv(spp_file, check.names = FALSE)[, c(1, 3)]
  if (!"valid_name" %in% colnames(spp)) {
    cli_alert_warning(paste("No 'valid_name' col in", spp_file, "- skipping basin"))
    return(NULL)
  }
  
  # 1. basin polygon and area -----------------------------------------------
  basin_row <- inland %>% dplyr::filter(basin_id == basin_id_value)
  if (nrow(basin_row) == 0) {
    cli_alert_warning(paste("Basin polygon not found for", basin_id_value, "- skipping"))
    return(NULL)
  }
  
  basin_area_km2 <- readRDS("input/data_prep/watershed_area.rds") %>%
    dplyr::filter(basin_id == basin_id_value) %>%
    dplyr::pull(area_km2) %>%
    as.numeric()
  if (length(basin_area_km2) == 0 || is.na(basin_area_km2)) {
    cli_alert_warning(paste("No area_km2 for", basin_id_value, "- skipping"))
    return(NULL)
  }
  
  # 2. build equal-area grid (~50 km) ---------------------------------------
  basin_proj <- st_transform(
    basin_row,
    crs = "+proj=eck4 +lon_0=0 +lat_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
  )
  
  fishnet_clipped_sf <- st_make_grid(
    basin_proj,
    cellsize = c(50000, 50000),
    what = "polygons"
  ) %>%
    suppressWarnings(st_make_valid()) %>%
    st_intersection(basin_proj) %>%
    st_as_sf() %>%
    dplyr::mutate(grid_id = dplyr::row_number())
  
  # container for all species rows ------------------------------------------
  all_results <- data.frame(
    basin_id        = character(),
    valid_name      = character(),
    sampling_year   = numeric(),
    relative_rarity = numeric(),
    stringsAsFactors = FALSE
  )
  
  # 3. loop over species -----------------------------------------------------
  for (species in unique(spp$valid_name)) {
    
    occ_file <- paste0(path_occ_dir, gsub(" ", "_", species), ".rds")
    
    if (!file.exists(occ_file)) {
      cli_alert_warning(paste("Occurrence file not found:", occ_file, "- skipping species"))
      
      # add placeholder row for this species so it's not lost
      all_results <- dplyr::bind_rows(
        all_results,
        data.frame(
          basin_id        = basin_id_value,
          valid_name      = species,
          sampling_year   = NA_real_,
          relative_rarity = NA_real_,
          stringsAsFactors = FALSE
        )
      )
      next
    }
    
    occ_raw <- readRDS(occ_file)
    if (!all(c("decimalLongitude", "decimalLatitude", "year") %in% names(occ_raw))) {
      cli_alert_warning(paste("Missing required columns in", occ_file, "for", species))
      
      all_results <- dplyr::bind_rows(
        all_results,
        data.frame(
          basin_id        = basin_id_value,
          valid_name      = species,
          sampling_year   = NA_real_,
          relative_rarity = NA_real_,
          stringsAsFactors = FALSE
        )
      )
      next
    }
    
    # Occurrences as sf
    occ_sf <- st_as_sf(
      occ_raw,
      coords = c("decimalLongitude", "decimalLatitude"),
      crs = 4326
    )
    occ_proj <- st_transform(occ_sf, st_crs(fishnet_clipped_sf))
    
    # If basin produced no usable grid cells, still add placeholder
    if (nrow(fishnet_clipped_sf) == 0) {
      all_results <- dplyr::bind_rows(
        all_results,
        data.frame(
          basin_id        = basin_id_value,
          valid_name      = species,
          sampling_year   = NA_real_,
          relative_rarity = NA_real_,
          stringsAsFactors = FALSE
        )
      )
      next
    }
    
    # Intersect occurrences with grid
    occ_by_cell_year <- suppressWarnings(
      st_intersection(fishnet_clipped_sf, occ_proj)
    )
    
    if (nrow(occ_by_cell_year) == 0) {
      # no overlaps at all in this basin for this species
      all_results <- dplyr::bind_rows(
        all_results,
        data.frame(
          basin_id        = basin_id_value,
          valid_name      = species,
          sampling_year   = NA_real_,
          relative_rarity = NA_real_,
          stringsAsFactors = FALSE
        )
      )
      next
    }
    
    occ_by_cell_year <- occ_by_cell_year %>%
      st_drop_geometry() %>%
      dplyr::filter(!is.na(year)) %>%
      dplyr::distinct(year, grid_id) %>%
      dplyr::arrange(year, grid_id)
    
    if (nrow(occ_by_cell_year) == 0) {
      # all 'year' were NA, so same fallback
      all_results <- dplyr::bind_rows(
        all_results,
        data.frame(
          basin_id        = basin_id_value,
          valid_name      = species,
          sampling_year   = NA_real_,
          relative_rarity = NA_real_,
          stringsAsFactors = FALSE
        )
      )
      next
    }
    # -------------------------------------------------------------------------
    # Calculate annual relative rarity (RR) for each species
    #
    # Cumulative knowledge approach:
    #   - For each year, consider ALL grid cells where the species has EVER
    #     been observed up to and including that year.
    #   - Let cum_range_size = number of unique occupied grid cells known
    #     by that year.
    #   - Convert that into an area fraction of the basin, then invert.
    #
    # Ecological meaning:
    #   • High RR  → species known from very few places (narrow/rare)
    #   • Low RR   → species known broadly (widespread)
    # -------------------------------------------------------------------------
    # Build cumulative occupied cell counts by year -------------------------
    cell_area_km2 <- 2500  # 50 km x 50 km grid cell
    
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
    
    # Compute occupancy_fraction and rarity ---------------------------------
    cum_df <- cum_df %>%
      group_by(sampling_year) %>%
      dplyr::mutate(
        basin_id   = basin_id_value,
        valid_name = species,
        occupancy_fraction = ifelse(
          basin_area_km2 > 0,
          pmin((cum_range_size * cell_area_km2) / basin_area_km2, 1),
          NA_real_
        ),
        # robust rarity transform: log10(1/occ + 1)
        relative_rarity = ifelse(
          !is.na(occupancy_fraction) & occupancy_fraction > 0,
          log10((1 / occupancy_fraction) + 1),
          NA_real_
        )
      ) %>%
      dplyr::select(
        basin_id,
        valid_name,
        sampling_year,
        relative_rarity
      )
    
    # append to all_results
    all_results <- dplyr::bind_rows(all_results, cum_df)
  } # end species loop
  
  # after loop: if still empty (extreme edge case), create 1 fallback row ----
  if (nrow(all_results) == 0) {
    all_results <- data.frame(
      basin_id        = basin_id_value,
      valid_name      = NA_character_,
      sampling_year   = NA_real_,
      relative_rarity = NA_real_,
      stringsAsFactors = FALSE
    )
  }
  
  # final select (columns are guaranteed to exist now)
  all_results <- all_results %>%
    dplyr::select(basin_id, valid_name, sampling_year, relative_rarity)
  
  # write CSV ---------------------------------------------------------------
  
  write.csv(all_results, out_file, row.names = FALSE)
  cli::cli_alert_success(paste("Finished basin", basin_id_value, "→ saved", out_file))
  
  invisible(all_results)
}

# ------------------------------------------------------------
# Run all basins
# ------------------------------------------------------------
purrr::walk(basin_ids, process_basin)

gc()
cli_progress_done()

