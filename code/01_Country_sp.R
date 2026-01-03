# =============================================================================
# Load required packages
# =============================================================================
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(sf)
library(cli)
library(rfishbase) # version 5.0.1
# Disable s2 spherical geometry (avoids issues with polygon intersections)
sf::sf_use_s2(FALSE)

# =============================================================================
# 0. Global settings and data input
# =============================================================================
OCC_DIR <- "input/processed/occ"

country_raw <- readRDS("input/raw/country.rds")
df_taxa     <- openxlsx::read.xlsx("input/raw/cas_freshwater_v1.xlsx") # species table
biog_list   <- read.csv("input/raw/biogeographic_list.csv")
iso_table   <- readxl::read_excel("input/raw/country_iso.xls")
inland_raw  <- readRDS("input/raw/basin/basin_sf_v1.rds")

# =============================================================================
# 1. Prepare and clean the country shapefile
# =============================================================================
prep_country_layer <- function(country_raw) {
  country_raw %>%
    st_transform(4326) %>%                        # ensure consistent global CRS (WGS84)
    suppressWarnings(st_make_valid()) %>%         # repair invalid geometries
    suppressWarnings(st_buffer(0)) %>%            # close small cracks or gaps
    select(iso3, geometry)                     # keep only essential columns
}

country_clean <- prep_country_layer(country_raw)

# =============================================================================
# 2. Extract country ISO3 codes for a single species
# =============================================================================
occurrence_countries <- function(species, occ_dir, country_sf) {
  occ_file <- file.path(
    occ_dir,
    paste0(gsub(" ", "_", species), ".rds")
  )
  
  # Skip species with no occurrence file
  if (!file.exists(occ_file)) {
    cli::cli_alert_warning("File not found for {.bold {species}}: {occ_file}")
    return(tibble(species = species, iso3 = NA_character_))
  }
  
  cli::cli_alert_info("Reading occurrence data for {.bold {species}}")
  occ_raw <- readRDS(occ_file)
  
  # ---- Convert to sf point geometry ----
  occ_sf <- tryCatch({
    if (inherits(occ_raw, "sf")) {
      # ensure CRS exists
      if (is.na(st_crs(occ_raw))) occ_raw <- st_set_crs(occ_raw, 4326)
      occ_raw
    } else {
      # assume GBIF-style longitude/latitude columns
      stopifnot(all(c("decimalLongitude", "decimalLatitude") %in% names(occ_raw)))
      occ_raw %>%
        filter(!is.na(decimalLongitude), !is.na(decimalLatitude)) %>%
        st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
                 crs = 4326, remove = FALSE)
    }
  }, error = function(e) {
    cli::cli_alert_warning("Failed to convert {.bold {species}} ({e$message})")
    return(NULL)
  })
  
  # Return empty if no valid points
  if (is.null(occ_sf) || nrow(occ_sf) == 0L) {
    cli::cli_alert_warning("No valid points for {.bold {species}}")
    return(tibble(species = species, iso3 = NA_character_))
  }
  
  occ_sf <- st_transform(occ_sf, 4326)
  
  # ---- Spatial join: assign country to each occurrence ----
  sf::sf_use_s2(FALSE)
  joined <- suppressWarnings(
    st_join(occ_sf, country_sf, join = st_within, left = TRUE)
  )
  
  # ---- Aggregate ISO3 per species ----
  out_tab <- joined %>%
    st_drop_geometry() %>%
    filter(!is.na(iso3)) %>%
    summarise(iso3 = paste0(sort(unique(iso3)), collapse = ";"), .groups = "drop")
  
  # Return NA if no match found
  if (nrow(out_tab) == 0L) {
    return(tibble(species = species, iso3 = NA_character_))
  }
  
  out_tab %>% mutate(species = species, .before = 1)
}

# =============================================================================
# 3. Extract country ISO3 codes for all species
# =============================================================================
all_iso_stats <- map_df(df_taxa$valid_name, function(sp) {
  occurrence_countries(species = sp, occ_dir = OCC_DIR, country_sf = country_clean)
})

# Optionally save and reload to avoid reprocessing
# saveRDS(all_iso_stats, "input/raw/all_iso_stats.rds")
# all_iso_stats <- readRDS("input/raw/all_iso_stats.rds")

# =============================================================================
# 4. Identify species missing ISO3 codes
# =============================================================================
remain_species <- all_iso_stats %>%
  mutate(iso3 = as.character(iso3)) %>%
  filter(is.na(iso3) | str_trim(iso3) == "") %>%
  pull(species)

# remain_species_fishbase <- country(remain_species)
# Optionally save and reload to avoid reprocessing
# rrr <- remain_species_fishbase %>% left_join(iso_table[,1:2],by = "country")
# openxlsx::write.xlsx(rrr,"input/raw/remain_species_fishbase.xlsx")

remain_species_fishbase <- openxlsx::read.xlsx("input/raw/remain_species_fishbase.xlsx")
length(unique(remain_species_fishbase$Species)) #2960
length(unique(remain_species_fishbase$country))
remain_species_fishbase <- remain_species_fishbase %>%
  select(species,country) %>%
  distinct()


rrrr <- setdiff(remain_species,remain_species_fishbase$species)

# =============================================================================
# 5. Infer ISO3 codes for missing species using basin-country mapping
# =============================================================================
# (a) Extract basins for missing species
df_missing <- df_taxa %>%
  filter(valid_name %in% rrrr) %>%
  mutate(basin = str_split(basin, ";"),
         basin = lapply(basin, unique),
         basin = lapply(basin, sort)) %>%
  unnest(basin) %>%
  select(valid_name, year_description, basin)

# (b) Link basins to basin_id
df_missing <- df_missing %>%
  left_join(biog_list, by = "basin") %>%
  select(basin_id, valid_name, year_description)

# (c) Map basin_id → ISO3 using country lookup
basin_iso <- inland_raw %>%
  st_drop_geometry() %>%
  select(basin_id, country) %>%
  mutate(country = str_split(country, ";"),
         country = lapply(country, unique),
         country = lapply(country, sort)) %>%
  unnest(country) %>%
  left_join(iso_table, by = "country") %>%
  group_by(basin_id) %>%
  summarise(iso3 = paste0(unique(iso3), collapse = ";"), .groups = "drop")

# (d) Aggregate inferred ISO3 by species
df_remain <- df_missing %>%
  left_join(basin_iso, by = "basin_id") %>%
  group_by(valid_name) %>%
  summarise(iso3 = paste0(unique(iso3), collapse = ";"), .groups = "drop") %>%
  mutate(iso3 = str_split(iso3, ";"),
         iso3 = lapply(iso3, unique),
         iso3 = lapply(iso3, sort)) %>%
  unnest(iso3) %>%
  group_by(valid_name) %>%
  summarise(iso3 = paste0(unique(iso3), collapse = ";"), .groups = "drop")

# =============================================================================
# 6. Merge inferred ISO3 values into the main table
# =============================================================================
remain_species_fishbase <- openxlsx::read.xlsx("input/raw/remain_species_fishbase.xlsx") %>%
  select(valid_name = species,iso3 = iso3) %>%
  distinct() %>%
  group_by(valid_name) %>%
  summarise(iso3 = paste0(unique(iso3), collapse = ";"), .groups = "drop")

df_remain_all <- rbind(df_remain,remain_species_fishbase)

all_iso <- all_iso_stats %>%
  left_join(df_remain_all, by = c("species" = "valid_name")) %>%
  mutate(iso3.x = as.character(iso3.x),
         iso3.y = as.character(iso3.y),
         iso3   = ifelse(is.na(iso3.x) | str_trim(iso3.x) == "", iso3.y, iso3.x)) %>%
  select(species, iso3)

# Check remaining missing ISO3 values
sum(is.na(all_iso$iso3))

# rr <- country_raw %>% filter(iso3 == "GRL")
# idx <- st_intersects(inland_raw, rr, sparse = FALSE)[, 1]
# basin_in_country <- inland_raw[idx, ]

all_iso$iso3 <- ifelse(all_iso$species == "Salvelinus alpinus",paste0(all_iso$iso3,";","GRL"),all_iso$iso3)
all_iso$iso3 <- ifelse(all_iso$species == "Salmo salar",paste0(all_iso$iso3,";","GRL"),all_iso$iso3)

# =============================================================================
# 7. Count the number of unique species per country
# =============================================================================
country_species <- all_iso %>%
  mutate(iso3 = str_split(iso3, ";")) %>%
  unnest(iso3) %>%
  mutate(iso3 = str_trim(iso3)) %>%
  group_by(iso3) %>%
  summarise(n_species = n_distinct(species), .groups = "drop") %>%
  arrange(desc(n_species))

head(country_species)

# =============================================================================
# 8. Attach country-level metadata and save results
# =============================================================================
country_fix <- country_raw %>%
  filter(iso3 %in% country_species$iso3) %>%
  left_join(country_species, by = "iso3") %>%
  select(iso3, country, region, income_group, n_species)

saveRDS(country_fix, "input/processed/country_fix.rds")

# =============================================================================
# 9. Generate species-country table and export one CSV per country
# =============================================================================
country_sp <- all_iso %>%
  left_join(df_taxa[, c("valid_name", "year_description")],
            by = c("species" = "valid_name")) %>%
  transmute(valid_name = species,
            iso3 = iso3,
            year_description = year_description)

saveRDS(country_sp, "input/processed/country_sp.rds")

# Create output directory if needed
out_dir <- "input/processed/country_sp"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Expand ISO3 codes into one row per country
final_long <- country_sp %>%
  mutate(iso3 = str_split(iso3, ";")) %>%
  unnest(iso3) %>%
  mutate(iso3 = str_trim(iso3))

# Quick validation
anyNA(final_long$iso3)
length(unique(final_long$iso3))       # expected: ~219
length(unique(final_long$valid_name)) # expected: ~18821


# Export one CSV per country
final_long %>%
  group_by(iso3) %>%
  group_walk(~ {
    this_iso3 <- .y$iso3[[1]]
    out_tbl   <- bind_cols(.y, .x)
    file_path <- file.path(out_dir, paste0(this_iso3, ".csv"))
    write.csv(out_tbl, file_path, row.names = FALSE)
    message("Saved: ", file_path)
  })


################################################################################
rm(list = ls())
gc()  # trigger garbage collection
library(furrr)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(sf)
library(cli)
library(magrittr)
sf::sf_use_s2(FALSE)
# Suppress all warnings
options(warn = -1)
country <- readRDS("input/processed/country_fix.rds")

# Configure parallel processing (single country here, but keep parallel setup)
#future::plan(future::multisession(workers = parallel::detectCores() - 30))  # Use available cores minus 30

# Get country IDs
iso3_id_value <- country %>%
  st_drop_geometry() %>%
  pull(iso3)

###########################-----------------------------------------------------

# Initialize progress bar
cli_progress_bar("Processing iso3 ID", total = length(iso3_id_value))

# --------------------------------------------------------------------
# Fast function to process one country (ISO3 code)
# --------------------------------------------------------------------
process_country <- function(iso3_id_value, country_sf) {
  
  cli_h2(paste("Processing country:", iso3_id_value))
  
  # ------------------------------------------------------------
  # 1. Subset the country polygon only once
  # ------------------------------------------------------------
  country_sel <- country_sf %>%
    dplyr::filter(iso3 %in% iso3_id_value)
  
  if (nrow(country_sel) == 0) {
    cli_alert_warning(paste("No geometry found for:", iso3_id_value, "- Skipping"))
    return(NULL)
  }
  
  # ------------------------------------------------------------
  # 2. Read country-level species list
  # ------------------------------------------------------------
  spp_file <- paste0("input/processed/country_sp/", iso3_id_value, ".csv")
  
  if (!file.exists(spp_file)) {
    cli_alert_warning(paste("Species file not found:", spp_file))
    return(NULL)
  }
  
  spp <- read.csv(spp_file)
  
  # Basic validation
  if (!"valid_name" %in% names(spp)) {
    cli_alert_warning("Column 'valid_name' missing - Skipping")
    return(NULL)
  }
  if (!"year_description" %in% names(spp)) {
    cli_alert_warning("Column 'year_description' missing - Skipping")
    return(NULL)
  }
  
  species_vec <- unique(spp$valid_name)
  
  # ------------------------------------------------------------
  # 3. Pre-allocate vector for geolocation years
  #    (avoid mutate inside loops)
  # ------------------------------------------------------------
  year_geo_vec <- rep(NA_integer_, length(species_vec))
  names(year_geo_vec) <- species_vec
  
  # Disable S2 for faster and more robust within() operations
  sf_use_s2(FALSE)
  
  # ------------------------------------------------------------
  # 4. Loop through species and compute earliest geolocation year
  # ------------------------------------------------------------
  for (i in seq_along(species_vec)) {
    species <- species_vec[i]
    
    # Build path to species occurrence file
    occ_file <- paste0(
      "input/processed/occ/",
      gsub(" ", "_", species),
      ".rds"
    )
    
    # Occurrence file missing
    if (!file.exists(occ_file)) {
      cli_alert_warning(paste("Missing occurrence file:", occ_file))
      year_geo_vec[i] <- NA_integer_
      next
    }
    
    cli_alert_info(paste("Reading occurrences for species:", species))
    
    distribution_sf <- readRDS(occ_file)
    
    # Ensure geometry exists
    if (!inherits(distribution_sf, "sf")) {
      distribution_sf <- st_as_sf(
        distribution_sf,
        coords = c("decimalLongitude", "decimalLatitude"),
        crs = 4326
      )
    }
    
    # Fast point-in-polygon filtering (avoid st_intersection)
    # pts_in_country <- st_join(
    #   distribution_sf,
    #   country_sel,
    #   join = st_within,
    #   left = FALSE
    # )
    
    idx <- sf::st_intersects(distribution_sf, country_sel, sparse = TRUE)
    keep <- lengths(idx) > 0
    pts_in_country <- distribution_sf[keep, , drop = FALSE]
    
    # Extract species description year
    current_year <- spp %>%
      dplyr::filter(valid_name == species) %>%
      dplyr::pull(year_description) %>%
      unique()
    
    if (length(current_year) == 0 || all(is.na(current_year))) {
      year_geo_vec[i] <- NA_integer_
      next
    }
    current_year <- current_year[1]
    
    # Extract geolocation years > description year
    if (!"year" %in% names(pts_in_country)) {
      cli_alert_warning(paste("Column 'year' missing in occurrences for", species))
      year_geo_vec[i] <- NA_integer_
      next
    }
    
    yrs_ok <- pts_in_country$year
    yrs_ok <- yrs_ok[!is.na(yrs_ok) & yrs_ok > current_year]
    
    # Determine earliest post-description geolocation year
    if (length(yrs_ok) == 0) {
      year_geo_vec[i] <- NA_integer_
    } else {
      year_geo_vec[i] <- min(yrs_ok)
    }
    
    cli_alert_info(paste("Finished species:", species))
  }
  
  # ------------------------------------------------------------
  # 5. Join results back to species table in one operation
  # ------------------------------------------------------------
  year_geo_df <- tibble::tibble(
    valid_name = species_vec,
    year_geolocation = as.integer(year_geo_vec)
  )
  
  spp_new <- spp %>%
    dplyr::left_join(year_geo_df, by = "valid_name")
  
  # Save updated file
  write.csv(spp_new, spp_file, row.names = FALSE)
  
  cli_alert_success(paste("Completed country:", iso3_id_value))
  
  gc()
  invisible(spp_new)
}

# Process all country IDs
purrr::walk(iso3_ids, process_country, country_sf = country)

# Finish progress bar
cli_progress_done()

# Shut down parallel workers and free memory
#future::plan(future::sequential)
gc()

################################################################################
library(furrr)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(sf)
library(cli)
library(arrow)
library(magrittr)
# Suppress all warnings
options(warn = -1)

country <- readRDS("input/processed/country_fix.rds")

# Configure parallel processing (single country here, but keep parallel setup)

# Get country IDs
iso3_id_value <- country %>%
  st_drop_geometry() %>%
  pull(iso3)

seq_annotation <- read_parquet("input/raw/seq_annotation_v1.parquet")


# Progress bar
cli_progress_bar("Processing country ID", total = length(iso3_id_value))

# Process a single country_id
process_country <- function(iso3_id_value) {
  
  # Log current country
  cli_h2(paste("Processing country:", iso3_id_value))
  
  # Species table for this country
  spp_file <- paste0("input/processed/country_sp/", iso3_id_value, ".csv")
  
  # Skip if file missing
  if (!file.exists(spp_file)) {
    cli_alert_warning(paste("File not found:", spp_file, "- Skipping"))
    return(NULL)
  }
  
  # Read species data
  spp <- read.csv(spp_file)
  
  
  # Validate data frame
  if (!is.data.frame(spp)) {
    cli_alert_warning("spp is not a data frame - Skipping")
    return(NULL)
  }
  
  # Ensure 'valid_name' exists
  if (!"valid_name" %in% colnames(spp)) {
    cli_alert_warning("'valid_name' column missing - Skipping")
    return(NULL)
  }
  
  # Iterate over species in this country
  for (species in unique(spp$valid_name)) {
    
    # Path to per-species NCBI/sequence metadata (CSV named by species)
    gb_file <- seq_annotation %>% dplyr::filter(valid_name == species)
    
    # If missing, set min_year to NA
    if (is_empty(gb_file)) {
      cli_alert_warning(paste("No sequence metadata found for species:", species, "- setting year_sequence = NA"))
      min_year <- NA
    } else {
      # Read NCBI/sequence data
      cli_alert_info("Reading NCBI data for species: {.bold {species}}")
      
      gb <- gb_file %>% 
        dplyr::select(valid_name, iso3, year) %>%
        na.omit() %>%
        distinct() %>%
        dplyr::filter(iso3 == iso3_id_value)
      
      # Year of description for this species (from spp table)
      current_year <- spp %>% 
        filter(valid_name == species) %>% 
        pull(year_description) %>% 
        unique()
      
      # Earliest sequence year later than description year
      min_year <- if (any(!is.na(gb$year) & gb$year > current_year)) {
        min(gb$year[gb$year > current_year], na.rm = TRUE)
      } else {
        NA
      }
    }
    
    # Ensure 'year_sequence' exists; then update for this species
    if (!"year_sequence" %in% colnames(spp)) {
      spp$year_sequence <- NA  # create if missing
    }
    
    spp %<>%
      mutate(
        year_sequence = ifelse(valid_name == species, min_year, year_sequence),
        iso3 = iso3_id_value
      )
    
    
    cli_alert_info("Completed species: {.bold {species}}, country: {.bold {iso3_id_value}}")
  }
  
  # Save updated species table (after all species are processed)
  write.csv(spp, file = spp_file, row.names = FALSE)
  
  cli_alert_info(paste("All species processed for country: {.bold {iso3_id_value}}"))
  
  # Clean up after each country
  gc()
}

# Run over all country
purrr::walk(iso3_id_value, process_country) 

# Finish progress
cli_progress_done()

# Final cleanup (parallel session memory is minimal here, no plan teardown needed)
gc()


