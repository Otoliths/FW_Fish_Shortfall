library(tibble)
library(scales)
library(ggplot2)
library(dplyr)        # Data manipulation verbs
library(brms)         # Bayesian multilevel models interface
library(arrow)


# ============================================================
# Load basin-level data and construct species × realm table
# ============================================================

# Basin-level dataset with at least:
# - valid_name           : species name (or other unique taxon ID)
# - biogeographic_realm  : realm code / factor
base_df <- read_parquet("input/data_prep/linnaean_basin_single.parquet")
# Species-level posterior probabilities
data_L <- readRDS("output/stan_survdata_basin_linnaean_single.rds")
fit <- readRDS("output/model/basin_linnaean_weibull_full.rds")
source("code/functions/xxx.r")
lin_res <- compute_Linnaean_prob(
  fit          = fit,
  year         = 2025,
  data         = data_L,
  include_data = FALSE,
  re_formula   = NULL,
  draws        = NULL,
  return_draws = TRUE
)

prob_df <- cbind(data_L,lin_res$summary)


realm_richness_weighted <- base_df %>%
  distinct(valid_name, biogeographic_realm) %>%
  add_count(valid_name, name = "n_realms") %>%
  mutate(weight = 1 / n_realms) %>%
  group_by(biogeographic_realm) %>%
  summarise(
    richness_weighted = sum(weight),
    richness_weighted_round = round(richness_weighted),
    .groups = "drop"
  )

# weights for (valid_name × realm), used to aggregate probabilities consistently
w_realm <- base_df %>%
  distinct(valid_name, biogeographic_realm) %>%
  add_count(valid_name, name = "n_realms") %>%
  transmute(valid_name, biogeographic_realm, weight = 1 / n_realms)

# realm posterior probabilities (weighted.mean)
realm_prob <- prob_df %>%
  select(valid_name, biogeographic_realm, prob_undesc_median, prob_undesc_lower, prob_undesc_upper) %>%
  left_join(w_realm, by = c("valid_name", "biogeographic_realm")) %>%
  dplyr::filter(!is.na(weight)) %>%
  group_by(biogeographic_realm) %>%
  summarise(
    prob_median  = weighted.mean(prob_undesc_median,  w = weight),
    prob_lower = weighted.mean(prob_undesc_lower, w = weight),
    prob_upper = weighted.mean(prob_undesc_upper, w = weight),
    .groups = "drop"
  )

realm_dat <- realm_richness_weighted %>%
  left_join(realm_prob, by = "biogeographic_realm")

eps <- 1e-12
realm_SR_summary <- realm_dat %>%
  mutate(
    prob_median  = pmax(pmin(prob_median,  1 - eps), 0),
    prob_lower = pmax(pmin(prob_lower, 1 - eps), 0),
    prob_upper = pmax(pmin(prob_upper, 1 - eps), 0),
    SRdesc    = richness_weighted / (1 - prob_median)  - richness_weighted,
    SRdesc_ll = richness_weighted / (1 - prob_lower) - richness_weighted,
    SRdesc_ul = richness_weighted / (1 - prob_upper) - richness_weighted,
    SRdesc    = round(pmax(SRdesc,    0)),
    SRdesc_ll = round(pmax(SRdesc_ll, 0)),
    SRdesc_ul = round(pmax(SRdesc_ul, 0))
  )

global_summary <- realm_SR_summary %>%
  summarise(
    region            = "Global",
    richness_weighted = sum(richness_weighted),
    SRdesc            = sum(SRdesc),
    SRdesc_ll         = sum(SRdesc_ll),
    SRdesc_ul         = sum(SRdesc_ul),
    .groups = "drop"
  ) %>%
  mutate(
    richness_weighted_round = round(richness_weighted),
    prob_median  = SRdesc    / (richness_weighted + SRdesc),
    prob_lower = SRdesc_ll / (richness_weighted + SRdesc_ll),
    prob_upper = SRdesc_ul / (richness_weighted + SRdesc_ul)
  )

final_realm_global <- dplyr::bind_rows(
  global_summary %>%
    dplyr::select(region, richness_weighted, prob_median, prob_lower, prob_upper, SRdesc, SRdesc_ll, SRdesc_ul),
  realm_SR_summary %>%
    dplyr::mutate(region = biogeographic_realm) %>%
    dplyr::select(region, richness_weighted, prob_median, prob_lower, prob_upper, SRdesc, SRdesc_ll, SRdesc_ul)
) %>%
  mutate(
    prob_display   = sprintf("%.3f (%.3f–%.3f)", prob_median, prob_lower, prob_upper),
    SRdesc_display = sprintf("%d (%d–%d)", SRdesc, SRdesc_ll, SRdesc_ul)
  ) %>%
  dplyr::select(
    Region = region,
    richness_weighted,
    Undescribed_probability = prob_display,
    Undescribed_species     = SRdesc_display
  )

write.csv(final_realm_global, paste0("output/tables/biogeographic_realm_linnaean_shortfall.csv"), row.names = FALSE)
rm(realm_richness_weighted,realm_prob,realm_dat,realm_SR_summary,
   global_summary)
################################################################################
#Basin-level probability of remaining undescribed
basin_richness_weighted <- base_df %>%
  distinct(valid_name, basin_id) %>%
  add_count(valid_name, name = "n_basin_id") %>%
  mutate(weight = 1 / n_basin_id) %>%
  group_by(basin_id) %>%
  summarise(
    richness_weighted = sum(weight),
    richness_weighted_round = round(richness_weighted),
    .groups = "drop"
  )

# weights for (valid_name × basin_id), used to aggregate probabilities consistently
w_basin <- base_df %>%
  distinct(valid_name, basin_id) %>%
  add_count(valid_name, name = "n_basin_id") %>%
  transmute(valid_name, basin_id, weight = 1 / n_basin_id)

# Basin posterior probabilities (weighted.mean)
basin_prob <- prob_df %>%
  select(valid_name, basin_id, prob_undesc_median, prob_undesc_lower, prob_undesc_upper) %>%
  left_join(w_basin, by = c("valid_name", "basin_id")) %>%
  dplyr::filter(!is.na(weight)) %>%
  group_by(basin_id) %>%
  summarise(
    prob_undesc_median  = weighted.mean(prob_undesc_median,  w = weight),
    prob_undesc_lower = weighted.mean(prob_undesc_lower, w = weight),
    prob_undesc_upper = weighted.mean(prob_undesc_upper, w = weight),
    .groups = "drop"
  )

# Join with scaled basin richness
basin_dat <- basin_richness_weighted %>%
  left_join(basin_prob, by = "basin_id")

eps <- 1e-12
basin_SR_summary <- basin_dat %>%
  mutate(
    prob_undesc_median  = pmax(pmin(prob_undesc_median,  1 - eps), 0),
    prob_undesc_lower = pmax(pmin(prob_undesc_lower, 1 - eps), 0),
    prob_undesc_upper = pmax(pmin(prob_undesc_upper, 1 - eps), 0),
    SRdesc_raw    = richness_weighted / (1 - prob_undesc_median)  - richness_weighted,
    SRdesc_ll_raw = richness_weighted / (1 - prob_undesc_lower) - richness_weighted,
    SRdesc_ul_raw = richness_weighted / (1 - prob_undesc_upper) - richness_weighted,
    SRdesc    = round(pmax(SRdesc_raw,    0)),
    SRdesc_ll = round(pmax(SRdesc_ll_raw, 0)),
    SRdesc_ul = round(pmax(SRdesc_ul_raw, 0))
  )

basin_SR_summary

write.csv(basin_SR_summary, paste0("output/tables/basin_linnaean_shortfall.csv"), row.names = FALSE)

######Drainage basins with x number of undescribed species----------------------

calculate_species_percentages <- function(data) {
  # Load required libraries
  library(dplyr)
  
  # Remove NA values from SRdesc
  SRdesc_values <- data$SRdesc[!is.na(data$SRdesc)]
  
  # Calculate total number of drainage basins (only non-NA SRdesc values are considered)
  total_basin <- n_distinct(data$basin_id[!is.na(data$SRdesc)])
  
  # Define fixed thresholds
  threshold_low <- 2   # Fixed lower threshold
  threshold_high <- 5  # Fixed upper threshold
  
  # Count the number of basins below the lower threshold
  basins_lt_threshold_low <- sum(data$SRdesc < threshold_low, na.rm = TRUE)
  
  # Count the number of basins above the higher threshold
  basins_gt_threshold_high <- sum(data$SRdesc > threshold_high, na.rm = TRUE)
  
  # Calculate the percentages based on total basin count
  percentage_lt_threshold_low <- round((basins_lt_threshold_low / total_basin) * 100, 2)
  percentage_gt_threshold_high <- round((basins_gt_threshold_high / total_basin) * 100, 2)
  
  # Print results
  cat(sprintf("\n%% Drainage basins with x < %.0f of undescribed species: %.2f%%",
              threshold_low, percentage_lt_threshold_low))
  cat(sprintf("\n%% Drainage basins with x > %.0f of undescribed species: %.2f%%\n",
              threshold_high, percentage_gt_threshold_high))
}

realm <- read.csv("input/raw/biogeographic_list.csv")
basin_SR_summary %>% 
  left_join(realm, by = "basin_id") %>%
  filter(!is.na(biogeographic_realm)) %>%         
  group_by(biogeographic_realm) %>%
  group_walk(~ {
    cat("\n========== Realm:", .y$biogeographic_realm, "==========\n")
    calculate_species_percentages(.x)
  })

calculate_species_percentages(basin_SR_summary)

rm(basin_richness_weighted,basin_prob,basin_dat,basin_SR_summary,realm)