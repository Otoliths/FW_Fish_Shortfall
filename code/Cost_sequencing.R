# The sequencing cost benchmark of $0.006 per Mb, as reported by NHGRI, was calibrated
# for each country based on GDP per capita. This correction accounts for regional
# economic disparities affecting laboratory infrastructure, reagent procurement,
# and technical capacity. The adjusted cost was estimated by scaling the U.S.-based
# reference cost inversely to national GDP per capita, providing a more realistic approximation of sequencing expenses across countries.

library(dplyr)
library(sf)
library(tidyr)
library(stringr)

country_iso <- readxl::read_excel("input/raw/country_iso.xls")
basin <- readRDS("input/raw/basin/basin_sf_v1.rds") %>%
  st_drop_geometry() %>%
  mutate(country = str_split(country, ";")) %>%     # split multiple basins
  mutate(country = lapply(country, unique)) %>%     # remove duplicate basin names
  mutate(country = lapply(country, sort)) %>%       # sort alphabetically for consistency
  unnest(country)  %>% 
  left_join(country_iso[,1:2], by = "country")


gdp_final <- read.csv("input/raw/gdp_2024.csv")
any(is.na(basin$iso3))
basin <- basin %>% left_join(gdp_final[,c(2,4)], by = "iso3")
base_ref <- filter(gdp_final,iso3 == "USA") %>% pull(gdp)
basin$ratio <- base_ref/basin$gdp

cost_s <- basin %>% 
  group_by(basin_id) %>%
  summarise(cost_seq_per = mean(0.006*ratio,na.rm = T))

cost_s <- cost_s %>% left_join(darwinian[,c(1,6)], by = "basin_id") %>%
  mutate(cost_sequencing = cost_seq_per*SRnoseq)

saveRDS(cost_s,"input/processed/basin_cost_sequencing.rds")

################################################################################
library(dplyr)
library(sf)

gdp <- read.csv("input/raw/gdp_2024.csv")
darwin <- read.csv("output/tables/country_darwinian_shortfall.csv")
country <- readRDS("input/processed/country_fix.rds") |> st_drop_geometry()

usa_gdp <- gdp |> filter(iso3 == "USA") |> pull(gdp) |> first()

cost_s <- country |>
  left_join(gdp |> select(iso3, gdp), by = "iso3") |>
  left_join(darwin |> select(iso3, SRnoseq), by = "iso3") |>
  mutate(
    ratio            = usa_gdp / gdp,
    cost_seq_per     = 0.006 * ratio,
    cost_sequencing  = cost_seq_per * SRnoseq
  )

saveRDS(cost_s, "input/processed/country_cost_sequencing.rds")
