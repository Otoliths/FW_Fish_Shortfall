# ------------------------------------------------------------------------------
# Project: Calculate accessibility (travel time to cities >50K population) for inland basins
# Core Reference: 
# Weiss, D. J., Nelson, A., Gibson, H. S., Temperley, W., Peedell, S., Lieber, A., … & Gething, P. W. (2018). 
# A global map of travel time to cities to assess inequalities in accessibility in 2015. Nature, 553(7688), 333-336.
# Data Source: https://figshare.com/articles/dataset/Travel_time_to_cities_and_ports_in_the_year_2015/7638134/4
#
# Nelson et al. 2019 Data Details (Travel Time Layers):
# File pattern: travel_time_to_cities_x.tif (x = 1-12; each = urban population-defined layer)
# Pixel value: Estimated travel time (minutes) to nearest urban area (2015)
# Key Population Ranges (for reference):
# x | Pop Min (≥) | Pop Max (<)    | Note
# 1 | 5,000,000   | 50,000,000     | -
# 2 | 1,000,000   | 5,000,000      | -
# 3 | 500,000			| 1,000,000      | -
# 4 | 200,000			| 500,000        | -
# 5 | 100,000			| 200,000        | -
# 6 | 50,000			| 100,000        | -
# 7 | 20,000			| 50,000         | -
# 8 | 10,000			| 20,000         | -
# 9 | 5,000				| 10,000         | -
# 10| 20,000			| 110,000,000    | -
# 11| 50,000      | 50,000,000     | ✅ Target layer (cities >50K; used in this workflow)
# 12| 5,000       | 110,000,000    | -
# ------------------------------------------------------------------------------

# Load required packages:
# - mapme.biodiversity: Spatial resource/indicator calculation workflow
# - sf: Vector spatial data handling (polygons, intersections)
# - tidyr/dplyr: Data reshaping and manipulation
# - terra: Raster processing backend
library(mapme.biodiversity)
library(sf)
library(tidyr)
library(dplyr)
library(terra)
library(ggplot2)

# Source custom script to fix bugs in mapme.biodiversity's get_nelson_et_al()/get_resources()
source("code/functions/zzz.R")

# Disable S2 geometry (avoids precision issues in global spatial operations)
options(sf_use_s2 = FALSE)

# --------------------------
# 1. Setup output directories
# --------------------------
# Main output dir for mapme.biodiversity; create if missing
outdir <- file.path(getwd(), "input/raw/mapme-data")
dir.create(outdir, showWarnings = FALSE)

# Subdir for Nelson et al. travel time raster (organize data)
nelson_dir <- file.path(outdir, "nelson_et_al")
dir.create(nelson_dir, showWarnings = FALSE)

# Configure mapme.biodiversity: Set output dir + disable verbose logging
mapme_options(
  outdir = outdir,
  verbose = FALSE
)

# --------------------------
# 2. Download & validate target travel time raster (x=11: Pop >50K)
# --------------------------
# Define local path for the downloaded raster (matches Pop >50K layer)
local_raster_path <- file.path(nelson_dir, "traveltime-50k_50mio.tif")

# Download from permanent Figshare URL (target: x=11 layer, Pop 50K-50M)
if (!file.exists(local_raster_path)) {
  dir.create(nelson_dir, recursive = TRUE, showWarnings = FALSE)  # Create folder quietly if needed
  download.file(
    url = "https://figshare.com/ndownloader/files/14189849",
    destfile = local_raster_path,
    mode = "wb"  # Critical for TIFF (binary file)
  )
  message("✅ Raster downloaded successfully.")
} else {
  message("✅ Raster file already exists (skip download).")
}

# Validate file integrity via MD5 checksum (ensure no corruption during download)
local_md5 <- tools::md5sum(local_raster_path)
target_md5 <- "befbda61c99903bd831147414fba7378"  # Expected checksum for x=11 layer

if (local_md5 == target_md5) {
  message("✅ File verification successful: MD5 checksum matches (x=11 layer, Pop >50K)")
} else {
  warning("❌ File verification failed: MD5 checksum mismatch\n",
          "Local checksum: ", local_md5, "\n",
          "Expected checksum: ", target_md5, "\n",
          "Re-download the x=11 layer from Figshare.")
}

# --------------------------
# 3. Prepare spatial input data (inland basins + country boundaries)
# --------------------------
# download global country boundaries
# https://geodata.ucdavis.edu/gadm/gadm4.1/gadm_410-levels.zip (The current version is 4.1)

st_layers("input/raw/gadm_410-levels/gadm_410-levels.gpkg")
# Driver: GPKG 
# Available layers:
#   layer_name geometry_type features fields crs_name
# 1      ADM_0 Multi Polygon      263      2   WGS 84
# 2      ADM_1 Multi Polygon     3662     11   WGS 84
# 3      ADM_2 Multi Polygon    47217     13   WGS 84
# 4      ADM_3 Multi Polygon   144193     16   WGS 84
# 5      ADM_4 Multi Polygon   153410     14   WGS 84
# 6      ADM_5 Multi Polygon    51427     15   WGS 84

country <- read_sf(dsn = "input/raw/gadm_410-levels/gadm_410-levels.gpkg",
                  layer = "ADM_0"  # Replace with the country layer name from st_layers()
                  )

# Load inland basin vector data (shapefile); drop 8th column (unused attribute)
inland <- st_read("input/raw/basin/InlandRealms_20250323.shp") %>%
  dplyr::select(-8)  # Remove unused 8th column

# Spatial intersection: Assign countries to each inland basin
# Keep key attributes (basin ID/name, GID_0/COUNTRY) + remove duplicate rows
basin_country <- st_intersection(inland, country) %>%
  dplyr::select(New_basin, basin_id, GID_0, COUNTRY) %>%
  distinct_all()

length(unique(basin_country$basin_id))
#3364

# --------------------------
# 4. Calculate accessibility statistics (travel time to cities >50K)
# --------------------------
# Workflow:
# 1. Fetch Nelson et al. travel time data (50K-50M population cities)
# 2. Compute zonal stats (mean/median/sd/min/max) with exactextract engine
# 3. Reshape results to long format (1 row per basin-statistic pair)
# Record start time of the accessibility calculation
start_time <- Sys.time()

result <- basin_country %>%
  get_resources(get_nelson_et_al_fix(ranges = "50k_50mio")) %>%
  calc_indicators(
    calc_traveltime(
      engine = "exactextract",
      stats = c("mean", "median", "sd", "min", "max")
    )
  ) %>%
  portfolio_long()

# Record end time and calculate duration (convert to hours)
end_time <- Sys.time()
duration_hours <- as.numeric(difftime(end_time, start_time, units = "secs")) / 3600

# Print concise time summary (formatted for readability)
cat(sprintf(
  "Calculation Complete:\nStart: %s\nEnd:   %s\nDuration: %.4f hours\n",
  format(start_time, "%Y-%m-%d %H:%M:%S"),
  format(end_time, "%Y-%m-%d %H:%M:%S"),
  duration_hours
))

# Calculation Complete:
# Start: 2025-10-05 22:51:54
# End:   2025-10-06 08:56:22
# Duration: 10.0745 hours

length(unique(result$basin_id))

write.csv(
  result,
  file.path("input/processed/cost_field_sampling", "basin_result.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------------------------
# Guidance on Selecting Accessibility Statistics (Aligned with Weiss et al. 2018 & Nelson et al. 2019)
# ------------------------------------------------------------------------------
# Context: The raster data (x=11 layer) provides pixel-level travel time to the NEAREST city with ≥50K population (2015).
# When aggregating these pixels to basin polygons, choose statistics based on your research goal:
#
# 1. Minimum (min) Travel Time
#    - Meaning: The shortest travel time from ANY point in the basin to a ≥50K city (e.g., a basin edge adjacent to a city).
#    - Limitation: Overly optimistic—does not represent the basin’s overall accessibility (ignores hard-to-reach areas).
#    - Use case: Only if your goal is to identify "best-case" accessibility (e.g., "Is there any part of the basin near a city?").
#
# 2. Maximum (max) Travel Time
#    - Meaning: The longest travel time from ANY point in the basin to a ≥50K city (e.g., a remote corner of a large basin like the Amazon).
#    - Limitation: Overly pessimistic—disproportionately influenced by extreme, sparsely populated areas.
#    - Use case: Only if your goal is to highlight "worst-case" accessibility (e.g., "What’s the hardest-to-reach area in the basin?").
#
# 3. Mean vs. Median Travel Time
#    - Meaning: Both represent "central tendency" (overall basin accessibility), but differ in robustness to outliers.
#    - Mean: Sensitive to extreme values (e.g., a single remote pixel can inflate the mean for small basins).
#    - Median: Insensitive to outliers—represents the "middle" travel time (half the basin has shorter times, half longer).
#    - Alignment with literature: Weiss et al. (2018) and Nelson et al. (2019) prioritize mean/median as default metrics for global accessibility analysis.
#
# Recommended Statistical Choice for Basin-Level Accessibility
# ------------------------------------------------------------------------------
# For most research goals (e.g., comparing overall accessibility across basins, assessing broad-scale inequalities):
# ✅ Primary Indicator: Median travel time to the nearest ≥50K city
#    - Rationale: Balances optimism/pessimism, avoids distortion from outliers, and aligns with standard practices in accessibility research.
#
# For supplementary context (to quantify spatial variability):
# ✅ Secondary Indicator: Standard deviation (sd) of travel time
#    - Rationale: Measures how much accessibility varies within the basin (e.g., a high sd means large differences between easy-to-reach and hard-to-reach areas).
#
# Example: To filter results to the recommended metrics for downstream analysis:

basin_clean <- result %>%
  # extract statistic from "variable" column (confirm prefix via head(result))
  # From your code, "variable" likely looks like "50k_50mio_traveltime_median" (adds "traveltime_")
  mutate(stat = stringr::str_remove(variable, "50k_50mio_traveltime_")) %>%
  # Filter to keep only median and sd (recommended stats)
  filter(stat %in% c("median", "sd")) %>%
  # Calculate traveltime_hours FIRST (before using it in select)
  mutate(
    traveltime_hours = value / 60,  # Convert minutes to hours
    TravelTime_Hours = round(traveltime_hours, 1)  # Round for readability (create final column here)
  ) %>%
  # Select columns (use the pre-calculated "TravelTime_Hours" directly)
  dplyr::select(
    Basin_Name = New_basin,
    basin_id = basin_id,
    Country = COUNTRY,
    Country_ISO = GID_0,
    Statistic = stat,
    TravelTime_Minutes = value,
    TravelTime_Hours,  # Now uses the rounded column from mutate()
    geometry
  )


# Check for missing values (e.g., NA in travel time or country)
cat("=== Data Quality Check ===\n")
missing_values <- basin_clean %>%
  st_drop_geometry() %>%
  summarize(
    Total_Rows = n(),
    Missing_Median_SD = sum(is.na(Statistic) | is.na(TravelTime_Minutes)),
    Missing_Country = sum(is.na(Country)),
    Unique_Basins = n_distinct(basin_id),
    Unique_Countries = n_distinct(Country)
  )
print(missing_values, row.names = FALSE)

# If there are missing values (e.g., NA in Statistic), filter them out
if (missing_values$Missing_Median_SD > 0) {
  basin_clean <- basin_clean %>%
    filter(!is.na(Statistic) & !is.na(TravelTime_Minutes))
  warning(glue::glue("Removed {missing_values$Missing_Median_SD} rows with missing median/sd data"))
}


write.csv(
  basin_clean %>% st_drop_geometry(),
  file.path("input/processed/cost_field_sampling", "basin_clean.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# Reshape data for non-spatial analysis (wide format: 1 row per basin)
# ------------------------------------------------------------------------------
# Calculate basin-level accessibility metrics (handling transboundary rivers with multiple countries)
basin_accessibility_wide <- basin_clean %>%
  st_drop_geometry() %>%  # Remove spatial info (keep tabular data only)
  # Aggregate by basin + statistic (average across countries in the same basin)
  group_by(basin_id, Statistic) %>%  # Group by basin and metric type (e.g., "median" or "sd")
  summarise(
    # Average travel times across all countries in the basin
    TravelTime_Minutes = mean(TravelTime_Minutes, na.rm = TRUE),
    TravelTime_Hours = mean(TravelTime_Hours, na.rm = TRUE),
    .groups = "drop"  # Ungroup after aggregation
  ) %>%
  # Reshape to wide format (columns for "median" and "sd")
  pivot_wider(
    names_from = Statistic,  # Spread "median"/"sd" into separate columns
    values_from = c(TravelTime_Minutes, TravelTime_Hours),
    names_glue = "{Statistic}_{.value}"  # Rename columns (e.g., "median_TravelTime_Hours")
  ) %>%
  # Classify accessibility and calculate variability
  mutate(
    # Categorize basins by median travel time
    Accessibility_Class = case_when(
      median_TravelTime_Hours <= 2 ~ "Very High (≤2h)",
      median_TravelTime_Hours <= 6 ~ "High (2-6h)",
      median_TravelTime_Hours <= 12 ~ "Medium (6-12h)",
      median_TravelTime_Hours <= 24 ~ "Low (12-24h)",
      TRUE ~ "Very Low (>24h)"
    ),
    # Coefficient of variation (relative variability: SD / Median)
    CV_TravelTime = case_when(
      is.na(median_TravelTime_Hours) | is.na(sd_TravelTime_Hours) ~ NA_real_,
      median_TravelTime_Hours == 0 ~ 0,  # Explicitly set CV to 0 when median is 0
      TRUE ~ round(sd_TravelTime_Hours / median_TravelTime_Hours, 2)
    )
  ) %>%
  # Sort by accessibility (most accessible first)
  arrange(median_TravelTime_Hours)

# Export wide-format data (for Excel/SPSS analysis)
write.csv(
  basin_accessibility_wide,
  file.path("input/processed/cost_field_sampling", "Basin_Accessibility_Wide_Format.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
#saveRDS(basin_accessibility_wide,"input/processed/basin_cost_field_sampling.rds")

# ------------------------------------------------------------------------------
# Spatial Visualization (Map of Basin Accessibility)
# ------------------------------------------------------------------------------
# Merge wide-format metrics back with geometry (for mapping)
basin_accessibility_spatial <- basin_accessibility_wide %>%
    left_join(unique(inland[,7]), by = "basin_id")


# Map of median travel time (choropleth map)
basin_map <- ggplot() +
  # Add base world map (low resolution for geographic context)
  geom_sf(
    data = ne_countries(scale = "medium", returnclass = "sf"), 
    fill = "gray90",    # Light gray for land background
    color = "white",    # White borders between countries
    size = 0.2          # Thin borders to avoid clutter
  ) +
  # Add basin layer (color-coded by accessibility class)
  geom_sf(
    data = basin_accessibility_spatial,
    aes(fill = Accessibility_Class, geometry = geometry),
    color = "white",    # White borders between basins
    size = 0.1          # Thin basin borders
  ) +
  # Custom color scale: FIXED order + intuitive color coding
  scale_fill_manual(
    # Critical: Explicitly set legend order (Very High → Very Low)
    limits = c("Very High (≤2h)", "High (2-6h)", "Medium (6-12h)", "Low (12-24h)", "Very Low (>24h)"),
    # Match colors to accessibility classes (green → red gradient)
    values = c(
      "Very High (≤2h)" = "#2E8B57",  # Dark green (most accessible)
      "High (2-6h)" = "#90EE90",      # Light green
      "Medium (6-12h)" = "#FFFFE0",   # Light yellow
      "Low (12-24h)" = "#FFB6C1",     # Light pink
      "Very Low (>24h)" = "#DC143C"   # Dark red (least accessible)
    ),
    # Legend title (line break for readability)
    name = "Accessibility Class\n(Median Travel Time)",
    # Optional: Add legend labels if you want to clarify further
    labels = c(
      "Very High (≤2 hours)",
      "High (2–6 hours)",
      "Medium (6–12 hours)",
      "Low (12–24 hours)",
      "Very Low (>24 hours)"
    )
  ) +
  scale_y_continuous(limits = c(-56,90))+
  theme_minimal() +
  theme(
    axis.text = element_blank(),        # Hide latitude/longitude labels (unnecessary for global map)
    axis.title = element_blank(),       # Hide axis titles
    legend.position = "bottom",         # Place legend at the bottom (avoids covering map)
    legend.key.width = unit(2, "cm"),   # Widen legend keys for easier color comparison
    legend.key.height = unit(0.8, "cm"),# Adjust legend key height
    legend.text = element_text(size = 9),# Smaller legend text to fit
    plot.title = element_text(          # Format main title
      size = 12, 
      face = "bold", 
      hjust = 0.5,
      margin = margin(b = 10)           # Add bottom margin to separate from subtitle
    ),
    plot.subtitle = element_text(       # Format subtitle
      size = 10, 
      hjust = 0.5, 
      color = "gray60",
      margin = margin(b = 15)           # Add bottom margin to separate from map
    ),
    plot.background = element_rect(     # White background (avoids transparent edges in reports)
      fill = "white", 
      color = NA
    )
  ) +
  # Map titles and captions (optional but informative)
  labs(
    title = "Accessibility of Inland Basins to Cities with ≥50K Population (2015)",
    subtitle = "Color-coded by Median Travel Time to the Nearest Major City",
    caption = "Data Source: Weiss et al. (2018) | Base Map: Natural Earth"  # Add data credit
  )

# Save the fixed map (high resolution for publications/presentations)
ggsave(
  filename = file.path("input/processed/cost_field_sampling", "Basin_Accessibility_Map.png"),
  plot = basin_map,
  width = 14,    # Wide enough for global map
  height = 8,    # Balanced height
  dpi = 300,     # High resolution (print quality)
  bg = "white",  # Ensure white background (critical for embedding in documents)
  device = "png" # Save as PNG (universally compatible)
)

# ------------------------------------------------------------------------------
# Country-Level Aggregation (Summarize by Country)
# ------------------------------------------------------------------------------
# Calculate country-level accessibility metrics (aggregates basin data to country scale)
country_accessibility <- basin_clean %>%
  st_drop_geometry() %>%  # Remove spatial geometry (keep tabular data)
  # Aggregate by Country + ISO + Statistic (average across basins in the same country)
  group_by(Country, Country_ISO, Statistic) %>%  # Group by country (ISO for standardization) + metric type (median/sd)
  summarise(
    # Average travel times across all basins within the country
    TravelTime_Minutes = mean(TravelTime_Minutes, na.rm = TRUE),
    TravelTime_Hours = mean(TravelTime_Hours, na.rm = TRUE),
    .groups = "drop"  # Clear grouping after aggregation
  ) %>%
  # Reshape to wide format (columns for "median" and "sd" metrics)
  pivot_wider(
    names_from = Statistic,  # Spread "median"/"sd" into separate columns
    values_from = c(TravelTime_Minutes, TravelTime_Hours),
    names_glue = "{Statistic}_{.value}"  # Rename columns (e.g., "median_TravelTime_Hours")
  ) %>%
  # Classify accessibility + handle CV_TravelTime NA (from zero median)
  mutate(
    # Categorize countries by median travel time (add "Immediate" for 0h)
    Accessibility_Class = case_when(
      median_TravelTime_Hours == 0 ~ "Immediate (0h)",  # Special class for instant access
      median_TravelTime_Hours <= 2 ~ "Very High (≤2h)",
      median_TravelTime_Hours <= 6 ~ "High (2-6h)",
      median_TravelTime_Hours <= 12 ~ "Medium (6-12h)",
      median_TravelTime_Hours <= 24 ~ "Low (12-24h)",
      TRUE ~ "Very Low (>24h)"
    ),
    # Calculate CV (relative variability) + fix NA from zero median/ missing values
    CV_TravelTime = case_when(
      is.na(median_TravelTime_Hours) | is.na(sd_TravelTime_Hours) ~ NA_real_,  # Missing inputs → NA
      median_TravelTime_Hours == 0 ~ 0,  # Zero median → no variability (CV = 0)
      TRUE ~ round(sd_TravelTime_Hours / median_TravelTime_Hours, 2)  # Normal CV calculation
    ),
    # Flag countries with zero median (for transparency in analysis)
    Zero_Median_Flag = median_TravelTime_Hours == 0
  ) %>%
  
  # Sort by accessibility (most accessible first)
  arrange(median_TravelTime_Hours)

# Print top 5 most accessible countries
cat("\n=== Top 5 Countries by Average Basin Accessibility ===\n")
print(head(country_accessibility, 5), row.names = FALSE)

# # Export country-level summary
write.csv(
  country_accessibility,
  file.path("input/processed/cost_field_sampling", "Country_Level_Accessibility_Summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ------------------------------------------------------------------------------
# Statistical Test (Optional: Compare Accessibility by Ecoregion)
# ------------------------------------------------------------------------------

# Merge Ecoregion from 'inland' (shapefile) to basin accessibility data
# Ensure join uses correct ID columns (basin_id in wide data ↔ basin_id in inland)
basin_with_ecoregion <- basin_clean %>%
  st_drop_geometry() %>%
  left_join(
    inland[, c("basin_id", "Ecoregion")],  # Extract only needed cols from inland
    by = c("basin_id" = "basin_id")        # Explicit join key (avoid ambiguity)
  ) %>%
  filter(!is.na(Ecoregion)) %>%           # Remove basins with missing ecoregion
  # Optional: Reorder Ecoregion factor to match your table order (for consistent plotting)
  mutate(
    Ecoregion = factor(
      Ecoregion,
      levels = c("Afrotropic", "Australasia", "Indomalayan", 
                "Nearctic", "Neotropic", "Oceania", "Palearctic")
    )
  )

# Verify Ecoregion distribution (match your provided table)
cat("=== Ecoregion Distribution (Basin Count) ===\n")
ecoregion_count <- table(basin_with_ecoregion$Ecoregion)
print(ecoregion_count)

# ------------------------------------------------------------------------------
# Ecoregion-Level Accessibility Summary (Key Metrics)
# ------------------------------------------------------------------------------
# Calculate core statistics per Ecoregion (accessibility class)
ecoregion_accessibility <- basin_with_ecoregion %>%
  # Step 1: Aggregate by Ecoregion + Statistic (average across basins in the same ecoregion)
  group_by(Ecoregion, Statistic) %>%  # Group by ecoregion + metric type (median/sd)
  summarise(
    # Average travel times across all basins within the ecoregion
    TravelTime_Minutes = mean(TravelTime_Minutes, na.rm = TRUE),
    TravelTime_Hours = mean(TravelTime_Hours, na.rm = TRUE),
    .groups = "drop"  # Clear grouping after aggregation
  ) %>%
  
  # Step 2: Reshape to wide format (columns for "median" and "sd")
  pivot_wider(
    names_from = Statistic,  # Spread "median"/"sd" into separate columns
    values_from = c(TravelTime_Minutes, TravelTime_Hours),
    names_glue = "{Statistic}_{.value}"  # Rename columns (e.g., "median_TravelTime_Hours")
  ) %>%
  
  # Step 3: Classify accessibility + handle CV_TravelTime (consistent with prior logic)
  mutate(
    # Categorize ecoregions by median travel time
    Accessibility_Class = case_when(
      median_TravelTime_Hours <= 2 ~ "Very High (≤2h)",
      median_TravelTime_Hours <= 6 ~ "High (2-6h)",
      median_TravelTime_Hours <= 12 ~ "Medium (6-12h)",
      median_TravelTime_Hours <= 24 ~ "Low (12-24h)",
      TRUE ~ "Very Low (>24h)"
    ),
    # Calculate CV with handling for zero median/missing values
    CV_TravelTime = case_when(
      is.na(median_TravelTime_Hours) | is.na(sd_TravelTime_Hours) ~ NA_real_,  # Missing data → NA
      median_TravelTime_Hours == 0 ~ 0,  # Zero median → no variability (CV = 0)
      TRUE ~ round(sd_TravelTime_Hours / median_TravelTime_Hours, 2)  # Standard CV calculation
    ),
    # Flag ecoregions with zero median travel time
    Zero_Median_Flag = median_TravelTime_Hours == 0
  ) %>%
  
  # Sort by accessibility (most accessible first)
  arrange(median_TravelTime_Hours)

# Export full Ecoregion summary to CSV (for reports/appendices)
write.csv(
  ecoregion_accessibility,
  file.path("input/processed/cost_field_sampling", "Ecoregion_Accessibility_Summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)


# ------------------------------------------------------------------------------
# Boxplot (Basin Median Time by Ecoregion)
# ------------------------------------------------------------------------------
# Shows distribution of basin-level median travel time within each Ecoregion
# Highlights median, quartiles, and outliers (extreme basins)
# Create faceted boxplots to compare travel time distributions by ecoregion AND statistic (e.g., median vs. sd)
ecoregion_boxplot_facet <- ggplot(
  basin_with_ecoregion,  # Use full dataset (no Statistic filter—keep both median/sd for faceting)
  aes(
    x = Ecoregion,        # X-axis: Ecoregions (predefined order from earlier code)
    y = TravelTime_Hours, # Y-axis: Travel time in hours
    fill = Ecoregion      # Fill boxplots with ecoregion-specific colors (for visual distinction)
  )
) +
  # Add boxplots: Adjust transparency (alpha) to avoid visual overload; minimize outliers
  geom_boxplot(
    alpha = 0.7,          # Soften fill color to prevent clashing with other elements
    outlier.size = 1,     # Small outlier points to reduce clutter
    outlier.alpha = 0.3   # Faint outliers so they don’t distract from main distribution
  ) +
  # Overlay mean value as a white dot (easier to compare central tendency vs. median in boxplots)
  stat_summary(
    fun = mean,           # Calculate mean for each ecoregion-statistic group
    geom = "point",       # Plot mean as a point
    shape = 21,           # Round point with border (to stand out against boxplots)
    size = 2.5,           # Point size (visible but not overwhelming)
    fill = "white",       # White fill for contrast
    color = "black"       # Black border to define the point clearly
  ) +
  # Facet by Statistic (split into side-by-side plots for median vs. sd)
  facet_wrap(
    vars(Statistic),      # Create one facet per unique Statistic (e.g., "median", "sd")
    nrow = 1,             # Arrange facets in a single row (horizontal comparison)
    scales = "free_x"     # Let x-axis adjust per facet (hides empty ecoregion slots if any)
  ) +
  # Color scale: Use Set3 (colorblind-friendly, distinct hues for 7 ecoregions)
  scale_fill_brewer(
    palette = "Set3",     # Predefined colorblind-safe palette (good for categorical data)
    name = "Ecoregion"    # Legend title (redundant later, but keeps scale definition clear)
  ) +
  # Y-axis: Format travel time (hours) and truncate to 95th percentile (reduce extreme outlier impact)
  scale_y_continuous(
    name = "Travel Time (Hours)",  # Y-axis label
    labels = scales::comma,        # Use comma formatting for large numbers (e.g., 1,000 vs. 1000)
    limits = c(0, quantile(        # Truncate y-axis to 0 and 95th percentile
      basin_with_ecoregion$TravelTime_Hours, 
      0.95, 
      na.rm = TRUE                 # Ignore NA values when calculating percentile
    )) 
  ) +
  # X-axis: Label ecoregions (no custom formatting beyond rotation)
  scale_x_discrete(name = "Ecoregion") +
  # Theme adjustments for readability (clean, publication-ready style)
  theme_minimal() +
  theme(
    axis.text.x = element_text(    # Format x-axis labels (ecoregion names)
      angle = 45,                  # Rotate labels 45° to avoid overlap
      hjust = 1,                   # Align label ends with tick marks
      size = 9                     # Smaller font to fit all names
    ),
    axis.text.y = element_text(size = 10),  # Slightly larger y-axis text for readability
    axis.title = element_text(              # Bold axis titles for emphasis
      size = 11, 
      face = "bold"
    ),
    legend.position = "none",               # Hide legend (redundant—x-axis already labels ecoregions)
    plot.title = element_text(              # Center and bold main title
      size = 12, 
      face = "bold", 
      hjust = 0.5
    ),
    plot.subtitle = element_text(           # Gray, centered subtitle (explains plot elements)
      size = 10, 
      hjust = 0.5, 
      color = "gray60"
    ),
    plot.background = element_rect(         # White background (avoids transparent edges in exports)
      fill = "white", 
      color = NA
    ),
    strip.text.x = element_text(            # Format facet titles (e.g., "median", "sd")
      size = 10, 
      face = "bold"                        # Bold to distinguish facet categories
    )
  ) +
  # Plot labels: Title, subtitle, and caption (explain truncation for transparency)
  labs(
    title = "Basin Accessibility Distribution by Ecoregion and Statistic",
    subtitle = "Box = Interquartile Range (IQR), White Dot = Mean, Outliers = Extreme Basins",
    caption = "Note: Y-axis truncated at 95th percentile for readability"
  )
# Save boxplot
ggsave(
  filename = file.path("input/processed/cost_field_sampling", "Ecoregion_Basin_TravelTime_Boxplot.png"),
  plot = ecoregion_boxplot_facet,
  width = 12, height = 8, dpi = 300, bg = "white"
)


# ------------------------------------------------------------------------------
# Stacked Bar Chart (Accessibility Class by Ecoregion)
# ------------------------------------------------------------------------------
# Shows percentage of basins in each accessibility class (per Ecoregion)
# First, reshape accessibility class data to long format
ecoregion_class_long <- basin_with_ecoregion %>%
  st_drop_geometry() %>%
  filter(Statistic == "median") %>%
  mutate(
    # Categorize ecoregions by median travel time
    Accessibility_Class = case_when(
      TravelTime_Hours <= 2 ~ "Very High (≤2h)",
      TravelTime_Hours <= 6 ~ "High (2-6h)",
      TravelTime_Hours <= 12 ~ "Medium (6-12h)",
      TravelTime_Hours <= 24 ~ "Low (12-24h)",
      TRUE ~ "Very Low (>24h)"
    )) %>%
  dplyr::select(Ecoregion, Accessibility_Class) %>%
  # Count basins per Ecoregion + Class
  count(Ecoregion, Accessibility_Class, name = "Basin_Count") %>%
  # Calculate percentage of total basins per Ecoregion
  group_by(Ecoregion) %>%
  mutate(Pct = round((Basin_Count / sum(Basin_Count)) * 100, 1)) %>%
  ungroup() %>%
  # Reorder Accessibility Class for logical stacking (Very High → Very Low)
  mutate(
    Accessibility_Class = factor(
      Accessibility_Class,
      levels = c("Very High (≤2h)", "High (2-6h)", "Medium (6-12h)", 
                "Low (12-24h)", "Very Low (>24h)")
    )
  )

# Create stacked bar chart
ecoregion_stacked_bar <- ggplot(ecoregion_class_long,
                                  aes(x = Ecoregion, 
                                      y = Pct,
                                      fill = Accessibility_Class)) +
  # Stacked bars with percentage labels (only show if ≥5% to avoid clutter)
  geom_col(position = "stack", alpha = 0.8) +
  geom_text(
    aes(label = ifelse(Pct >= 5, paste0(Pct, "%"), "")),  # Hide small labels
    position = position_stack(vjust = 0.5),
    size = 3, color = "black"
  ) +
  # Custom color scale (match earlier map for consistency)
  scale_fill_manual(
    values = c(
      "Very High (≤2h)" = "#2E8B57",
      "High (2-6h)" = "#90EE90",
      "Medium (6-12h)" = "#FFFFE0",
      "Low (12-24h)" = "#FFB6C1",
      "Very Low (>24h)" = "#DC143C"
    ),
    name = "Accessibility Class\n(Median Travel Time)",
    # Reorder legend to match stacking order
    limits = c("Very High (≤2h)", "High (2-6h)", "Medium (6-12h)", 
              "Low (12-24h)", "Very Low (>24h)")
  ) +
  # Y-axis: Percentage (0-100)
  scale_y_continuous(
    name = "Percentage of Basins (%)",
    labels = scales::percent_format(scale = 1)  # Show as "50%" instead of "0.5"
  ) +
  # X-axis: Rotate labels
  scale_x_discrete(name = "Ecoregion") +
  # Theme adjustments
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "right",
    legend.key.height = unit(1.2, "cm"),  # Taller legend keys for color comparison
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray60"),
    plot.background = element_rect(fill = "white", color = NA)
  ) +
  # Plot labels
  labs(
    title = "Accessibility Class Distribution by Ecoregion",
    subtitle = "Only labels for classes ≥5% shown to avoid clutter",
    caption = "Data: Weiss et al. (2018) | Ecoregion Classification: Inland Basins Shapefile"
  )

# Save stacked bar chart
ggsave(
  filename = file.path("input/processed/cost_field_sampling", "Ecoregion_Accessibility_Class_StackedBar.png"),
  plot = ecoregion_stacked_bar,
  width = 14, height = 9, dpi = 300, bg = "white"
)


# ------------------------------------------------------------------------------
# Statistical Test – Kruskal-Wallis H Test (Ecoregion Differences)
# ------------------------------------------------------------------------------
# Test if median travel time differs across Ecoregions (non-parametric, since data may not be normal)
# Null hypothesis: All Ecoregions have the same median travel time distribution
kruskal_test <- kruskal.test(
  formula = TravelTime_Hours ~ Ecoregion,
  data = basin_with_ecoregion %>% filter(Statistic == "median")
)

cat("\n=== Kruskal-Wallis Test Results (Ecoregion Accessibility Differences) ===\n")
print(kruskal_test)

# Kruskal-Wallis rank sum test
# 
# data:  TravelTime_Hours by Ecoregion
# Kruskal-Wallis chi-squared = 441.61, df = 6, p-value < 2.2e-16


# If significant (p < 0.05), run post-hoc tests (Dunn's test) to find which Ecoregions differ
if (kruskal_test$p.value < 0.05) {
  # Install if needed: install.packages("dunn.test")
  library(dunn.test)
  
  # Post-hoc Dunn's test (with Bonferroni correction for multiple comparisons)
  dunn_posthoc <- dunn.test(
    x = basin_with_ecoregion[which(basin_with_ecoregion$Statistic == "median"),]$TravelTime_Hours,
    g = basin_with_ecoregion[which(basin_with_ecoregion$Statistic == "median"),]$Ecoregion,
    method = "bonferroni",  # Controls family-wise error rate
    altp = TRUE  # Show adjusted p-values
  )
  
  # Export post-hoc results to text file
  sink(file.path(outdir, "Ecoregion_KruskalWallis_Posthoc.txt"))
  cat("Kruskal-Wallis H Test Results (Ecoregion Accessibility Comparison)\n")
  cat("Null Hypothesis: All Ecoregions have identical median travel time distributions\n")
  cat("Alternative Hypothesis: At least one Ecoregion differs\n")
  cat("-------------------------------------------------------------------------\n")
  cat("Kruskal-Wallis Test Statistic:", round(kruskal_test$statistic, 3), "\n")
  cat("Kruskal-Wallis p-value:", format(kruskal_test$p.value, scientific = TRUE), "\n\n")
  cat("Post-hoc Dunn's Test (Bonferroni Correction for Multiple Comparisons)\n")
  cat("-------------------------------------------------------------------------\n")
  print(dunn_posthoc)
  sink()  # Close text sink
  
  message("✅ Kruskal-Wallis + post-hoc results exported to: ", file.path(outdir, "Ecoregion_KruskalWallis_Posthoc.txt"))
} else {
  message("❌ Kruskal-Wallis test not significant (p ≥ 0.05) – no post-hoc tests needed")
}