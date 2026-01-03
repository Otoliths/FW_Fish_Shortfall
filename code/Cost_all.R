library(dplyr)

min_max_normalize_safe <- function(x, epsilon = 1e-6) {
  min_val <- min(x,na.rm = T)
  max_val <- max(x,na.rm = T)
  range_val <- max_val - min_val + epsilon
  return((x - min_val) / range_val)
}


cost_t <- readRDS("input/processed/basin_cost_taxonomy.rds")
cost_f <- readRDS("input/processed/basin_cost_field_sampling.rds")
cost_s <- readRDS("input/processed/basin_cost_sequencing.rds")

cost <- cost_t %>% select(basin_id = basin_id,cost_t = cost_years) %>%
  left_join(cost_f %>% select(basin_id = basin_id,cost_f = median_TravelTime_Hours), by = "basin_id") %>%
  left_join(cost_s %>% select(basin_id = basin_id,cost_s = cost_sequencing), by = "basin_id")
  
cost$cost_t_norm <- rank(cost$cost_t, ties.method="average") / max(rank(cost$cost_t))
cost$cost_f_norm <- rank(cost$cost_f, ties.method="average") / max(rank(cost$cost_f))
cost$cost_s_norm <- rank(cost$cost_s, ties.method="average") / max(rank(cost$cost_s))

cost_all <- cost %>%
  mutate(cost_scale = (cost_t_norm+cost_f_norm+cost_s_norm)/3)

range(cost_all$cost_scale)

saveRDS(cost_all,"input/processed/cost_basin_all.rds")
################################################################################
library(dplyr)

min_max_normalize_safe <- function(x, epsilon = 1e-6) {
  min_val <- min(x,na.rm = T)
  max_val <- max(x,na.rm = T)
  range_val <- max_val - min_val + epsilon
  return((x - min_val) / range_val)
}


cost_t <- readRDS("input/processed/country_cost_taxonomy.rds")
cost_f <- readRDS("input/processed/country_cost_field_sampling.rds")
cost_s <- readRDS("input/processed/country_cost_sequencing.rds")

cost <- cost_t %>% select(iso3 = iso3,cost_t = cost_years) %>%
  left_join(cost_f %>% select(iso3 = Country_ISO,cost_f = median_TravelTime_Hours), by = "iso3") %>%
  left_join(cost_s %>% select(iso3 = iso3,cost_s = cost_sequencing), by = "iso3")

cost$cost_t_norm <- rank(cost$cost_t, ties.method="average") / max(rank(cost$cost_t))
cost$cost_f_norm <- rank(cost$cost_f, ties.method="average") / max(rank(cost$cost_f))
cost$cost_s_norm <- rank(cost$cost_s, ties.method="average") / max(rank(cost$cost_s))

cost_all <- cost %>%
  mutate(cost_scale = (cost_t_norm+cost_f_norm+cost_s_norm)/3)

range(cost_all$cost_scale)

saveRDS(cost_all,"input/processed/cost_country_all.rds")
