library(sf)
library(dplyr)

# Define equal-area projection (Eckert IV)
crs_equal_area <- "+proj=eck4 +lon_0=0 +lat_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# 1. Read basin polygons in lon/lat (assumed WGS84) and fix geometry
country_raw <- readRDS("input/raw/country.rds") %>%
  suppressWarnings(st_make_valid()) %>%
  # Handle polygons that cross the antimeridian (wrap around dateline)
  st_wrap_dateline(
    options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"),
    quiet   = TRUE
  ) %>%
  # If multipolygons exist, extract polygon parts
  st_collection_extract("POLYGON", warn = FALSE) %>%
  st_as_sf()

# 2. Project to equal-area CRS for correct area calculation
country_equal_area <- st_transform(country_raw, crs_equal_area)

# 3. Compute area per iso3
#    st_area() returns units in m^2 because the CRS uses meters.
country_area <- country_equal_area %>%
  dplyr::filter(!is.na(.data$iso3)) %>%
  mutate(area_m2 = st_area(geometry),
         area_km2 = as.numeric(area_m2) / 1e6) %>%  # convert to km^2
  st_drop_geometry() %>%
  select(iso3,area_m2, area_km2)

saveRDS(country_area,"input/data_prep/country_area.rds")

