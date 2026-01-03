# ==============================================================
# Function: clean_species_occ
# Purpose : Clean and spatially filter occurrence records for a single species.
#            If the cleaned file already exists, skip processing.
# ==============================================================

clean_species_occ <- function(
    sp,
    data,             # Data frame linking species to basins (columns: valid_name, basin)
    polygon,          # sf polygon layer of basins (must contain 'basin')
    occ,              # Occurrence records with: species, decimalLongitude, decimalLatitude, year, basisOfRecord
    remove_fossil = TRUE,  # Remove fossil specimen records
    drop_na_year = TRUE,   # Drop records missing year
    out_dir = "input/processed/occ"  # Directory to save cleaned results
) {
  # ---- 0. Prepare output path and check if species already processed ----
  if (!dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  safe_sp <- sp %>%
    stringi::stri_trans_general("Latin-ASCII") %>%  # Remove non-ASCII characters
    gsub("[^A-Za-z0-9_\\-]+", "_", .)               # Replace symbols with underscores
  out_file <- file.path(out_dir, paste0(safe_sp, ".rds"))
  
  # If file already exists, skip processing
  if (file.exists(out_file)) {
    message("✅ Already processed: ", sp)
    return(invisible(NULL))
  }
  
  # ---- 1. Get basin range for this species ----
  basins <- data %>%
    filter(valid_name == sp) %>%
    distinct(basin) %>%
    pull(basin) %>%
    na.omit()
  if (length(basins) == 0) {
    message("⚠️ No basin found for species: ", sp)
    return(invisible(NULL))
  }
  
  sp_range <- polygon %>%
    filter(basin %in% basins) %>%
    st_make_valid() %>%
    st_union() %>%
    st_as_sf()
  
  # Ensure the CRS is WGS84
  if (is.na(st_crs(sp_range)) || st_crs(sp_range) != 4326)
    sp_range <- st_transform(sp_range, 4326)
  
  # ---- 2. Filter occurrence data for this species ----
  x <- occ %>%
    filter(species == sp) %>%
    filter(!is.na(decimalLongitude), !is.na(decimalLatitude))
  
  if (remove_fossil && "basisOfRecord" %in% names(x))
    x <- x %>% filter(basisOfRecord != "FOSSIL_SPECIMEN")
  if (drop_na_year && "year" %in% names(x))
    x <- x %>% filter(!is.na(year))
  if (nrow(x) == 0) {
    message("⚠️ No valid occurrences after filtering for: ", sp)
    return(invisible(NULL))
  }
  
  # ---- 3. Clean coordinates ----
  clean <- x %>%
    cc_val(lon = "decimalLongitude", lat = "decimalLatitude") %>%
    cc_equ(lon = "decimalLongitude", lat = "decimalLatitude") %>%
    cc_zero(lon = "decimalLongitude", lat = "decimalLatitude") %>%
    cc_dupl(lon = "decimalLongitude", lat = "decimalLatitude", additions = "year") %>%
    cc_outl(method = "distance",tdi = 1600)
  # outlier population from main range (~1,600 kilometers away)
  # https://www.natureserve.org/sites/default/files/eo_rank_specifications-generic_guidelines_and_decision_key_05.08.2020.pdf
  if (nrow(clean) == 0) {
    message("⚠️ All coordinates removed after cleaning for: ", sp)
    return(invisible(NULL))
  }
  
  # ---- 4. Keep only records within the basin polygon ----
  pts <- st_as_sf(clean,
                  coords = c("decimalLongitude", "decimalLatitude"),
                  crs = 4326, remove = FALSE)
  inside <- lengths(st_intersects(pts, sp_range)) > 0
  clean_sf <- pts[inside, , drop = FALSE] %>% st_drop_geometry()
  if (nrow(clean_sf) == 0) {
    message("⚠️ No points inside polygon for: ", sp)
    return(invisible(NULL))
  }
  
  # ---- 5. Save cleaned data ----
  saveRDS(clean_sf, out_file)
  message("💾 Saved cleaned file for: ", sp)
  
  invisible(clean_sf)
}
