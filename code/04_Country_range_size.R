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
path_output_dir    <- "input/data_prep/country_range_size/"

# Make sure output directory exists
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
  # Wrap around the dateline to avoid huge self-crossing polygons
  st_wrap_dateline(
    options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
    quiet = TRUE
  ) %>%
  st_collection_extract("POLYGON", warn = FALSE) %>%
  st_as_sf()

# Get all country IDs
country_ids <- country %>%
  st_drop_geometry() %>%
  pull(iso3)

# Progress bar for all basins
cli_progress_bar("Processing basins", total = length(country_ids))


# ------------------------------------------------------------
# Function: process a single country
# ------------------------------------------------------------
process_country <- function(country_id_value) {
  
  cli_h2(paste("Processing country:", country_id_value))
  out_file <- paste0(path_output_dir, country_id_value, ".csv")
  
  # # Skip if already processed
  # if (file.exists(out_file)) {
  #   cli::cli_alert_info("File already exists: {out_file} - skipping country")
  #   return(invisible(out_file))
  # }
  
  # -------------------------------------------------
  # 1. Species table
  # -------------------------------------------------
  spp_file <- paste0(path_species_dir, country_id_value, ".csv")
  if (!file.exists(spp_file)) {
    cli_alert_warning(paste("File not found:", spp_file, "- skipping"))
    return(NULL)
  }
  
  spp <- read.csv(spp_file, stringsAsFactors = FALSE)[, c(1,2)]
  
  if (!is.data.frame(spp)) {
    cli_alert_warning("Species table is not a data frame - skipping")
    return(NULL)
  }
  if (!"valid_name" %in% colnames(spp)) {
    cli_alert_warning("'valid_name' column missing - skipping")
    return(NULL)
  }
  
  # make sure we have a mutable column to fill
  if (!"range_size" %in% colnames(spp)) {
    spp$range_size <- NA_real_
  }
  
  # -------------------------------------------------
  # 2. Country polygon (get once)
  # -------------------------------------------------
  country_row <- country %>% dplyr::filter(.data$iso3 %in% country_id_value)
  if (nrow(country_row) == 0L) {
    cli_alert_warning("No polygon for country {country_id_value} - skipping")
    return(NULL)
  }
  
  # Heal geometry in lon/lat first to avoid GEOS explosions
  country_lonlat <- suppressWarnings(sf::st_transform(country_row, 4326)) %>%
    suppressWarnings(sf::st_make_valid()) %>%
    # buffer(0) fixes self-intersections
    suppressWarnings(sf::st_buffer(0))
  
  # Project once to equal-area (Eckert IV)
  country_proj <- suppressWarnings(
    sf::st_transform(
      country_lonlat,
      crs = "+proj=eck4 +lon_0=0 +lat_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    )
  ) %>%
    suppressWarnings(sf::st_make_valid()) %>%
    suppressWarnings(sf::st_buffer(0))
  
  # Drop degenerate slivers (area 0 or near 0) that sometimes appear
  poly_area <- suppressWarnings(sf::st_area(country_proj))
  keep_idx  <- as.numeric(poly_area) > 1  # m² threshold
  if (!any(keep_idx)) {
    cli_alert_warning("Country geometry degenerated for {country_id_value} - skipping")
    return(NULL)
  }
  country_proj <- country_proj[keep_idx, , drop = FALSE]
  
  # -------------------------------------------------
  # 3. Build fishnet grid ONCE for this country
  # -------------------------------------------------
  country_proj2 <- country_proj %>%
    st_make_valid() %>%
    st_buffer(0) %>%                       
    st_collection_extract("POLYGON") %>%   
    st_set_precision(100) %>%              
    st_make_valid()
  country_proj2 <- st_union(country_proj2) %>% st_as_sf()
  fishnet_clipped_sf <- tryCatch(
    {
      # big memory step; do it once only
      grid_raw <- sf::st_make_grid(
        country_proj,
        cellsize = c(50000, 50000),  # 50 km
        what     = "polygons"
      )
      
      grid_sf <- sf::st_as_sf(grid_raw) %>%
        suppressWarnings(sf::st_make_valid()) %>%
        suppressWarnings(sf::st_buffer(0))
      
      clipped <- suppressWarnings(
        sf::st_intersection(grid_sf, country_proj2)
      ) %>%
        sf::st_as_sf() %>%
        dplyr::mutate(grid_id = dplyr::row_number())
      
      clipped
    },
    error = function(e) {
      cli_alert_warning(paste(
        "Failed to build fishnet for", country_id_value, ":", conditionMessage(e),
        "→ Will write NA range_size for all species."
      ))
      NULL
    }
  )
  
  # If we could not build a grid at all (NULL or 0 rows), just write NA rows
  if (is.null(fishnet_clipped_sf) || nrow(fishnet_clipped_sf) == 0L) {
    cli_alert_warning(
      "No usable fishnet for {country_id_value}; writing NA range_size for all species."
    )
    write.csv(spp, file = out_file, row.names = FALSE)
    cli_alert_success(
      "All species processed for country {.bold {country_id_value}} → saved to {.path {out_file}}"
    )
    gc()
    return(invisible(out_file))
  }
  
  # -------------------------------------------------
  # 4. Loop species and update range_size
  # -------------------------------------------------
  for (species in unique(spp$valid_name)) {
    
    cli_alert_info("Reading occurrence data for species: {.bold {species}}")
    
    distribution_file <- paste0(
      path_occ_dir,
      gsub(" ", "_", species),
      ".rds"
    )
    
    if (!file.exists(distribution_file)) {
      cli_alert_warning(
        "Occurrence file not found: {distribution_file} - skipping {.bold {species}}"
      )
      next
    }
    
    # Load occurrences and convert to sf, with safety
    distribution_sf <- tryCatch(
      {
        readRDS(distribution_file) %>%
          sf::st_as_sf(
            coords = c("decimalLongitude", "decimalLatitude"),
            crs = 4326
          )
      },
      error = function(e) {
        cli_alert_warning(
          "Failed to read/convert occurrences for {.bold {species}}: {conditionMessage(e)}"
        )
        NULL
      }
    )
    
    if (is.null(distribution_sf) || nrow(distribution_sf) == 0L) {
      # leave this species' range_size as NA
      cli_alert_info(
        "No valid occurrences for {.bold {species}} in {.bold {country_id_value}}"
      )
      next
    }
    
    # Reproject occurrences to match the grid CRS
    points_sf_in_country <- tryCatch(
      {
        suppressWarnings(sf::st_transform(distribution_sf, sf::st_crs(fishnet_clipped_sf)))
      },
      error = function(e) {
        cli_alert_warning(
          "st_transform failed for {.bold {species}} in {.bold {country_id_value}}: {conditionMessage(e)}"
        )
        NULL
      }
    )
    
    if (is.null(points_sf_in_country) || nrow(points_sf_in_country) == 0L) {
      next
    }
    
    # Count how many points fall in each grid cell
    n_colli_vec <- tryCatch(
      {
        lengths(sf::st_intersects(fishnet_clipped_sf, points_sf_in_country))
      },
      error = function(e) {
        cli_alert_warning(
          "st_intersects failed for {.bold {species}} in {.bold {country_id_value}}: {conditionMessage(e)}"
        )
        NULL
      }
    )
    
    if (is.null(n_colli_vec)) {
      next
    }
    
    # Keep only cells with at least one record
    occupied_n <- sum(n_colli_vec > 0, na.rm = TRUE)
    
    # Update range_size for this species
    spp$range_size[spp$valid_name == species] <- occupied_n
    
    cli_alert_info(
      "Finished {.bold {species}} in country {.bold {country_id_value}} (occupied cells: {occupied_n})"
    )
  }
  
  # -------------------------------------------------
  # 5. Write country table
  # -------------------------------------------------
  write.csv(spp, file = out_file, row.names = FALSE)
  
  cli_alert_success(
    "All species processed for country {.bold {country_id_value}} → saved to {.path {out_file}}"
  )
  
  gc()
  invisible(out_file)
}


# ------------------------------------------------------------
# Run for all country
# ------------------------------------------------------------
purrr::walk(country_ids, process_country)

cli_progress_done()
gc()
