library(sf)
library(dplyr)

# Define equal-area projection (Eckert IV)
crs_equal_area <- "+proj=eck4 +lon_0=0 +lat_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# 1. Read basin polygons in lon/lat (assumed WGS84) and fix geometry
inland_raw <- readRDS("input/raw/basin/basin_sf_v1.rds") %>%
  #suppressWarnings(st_make_valid()) %>%
  # Handle polygons that cross the antimeridian (wrap around dateline)
  st_wrap_dateline(
    options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
    quiet   = TRUE
  ) %>%
  # If multipolygons exist, extract polygon parts
  #st_collection_extract("POLYGON", warn = FALSE) %>%
  st_as_sf()

# 2. Project to equal-area CRS for correct area calculation
inland_equal_area <- st_transform(inland_raw, crs_equal_area)

# 3. Compute area per basin_id
#    st_area() returns units in m^2 because the CRS uses meters.
basin_area <- inland_equal_area %>%
  mutate(area_m2 = st_area(geometry),
         area_km2 = as.numeric(area_m2) / 1e6) %>%  # convert to km^2
  st_drop_geometry() %>%
  select(basin_id, area_m2, area_km2)

saveRDS(basin_area,"input/data_prep/watershed_area.rds")

basin_area %>%
  filter(area_km2 < 2500) %>%
  select(basin_id, area_km2)

