library(tibble)
library(scales)
library(ggplot2)
library(dplyr)        # Data manipulation verbs
library(brms)         # Bayesian multilevel models interface
library(arrow)

################################################################################
# ------------------------------------------------------------------------------
# Estimating Linnaean shortfalls (number of undescribed species) at
# continent, country and family levels
# ------------------------------------------------------------------------------

# ============================================================
# Load country-level data and construct species × realm table
# ============================================================

# country-level dataset with at least:
# - valid_name           : species name (or other unique taxon ID)
# - continent  : code / factor
base_df <- read_parquet("input/data_prep/linnaean_country_single.parquet")
# Species-level posterior probabilities
data_L <- readRDS("output/stan_survdata_country_linnaean_single.rds")
source("code/functions/xxx.r")
lin_res <- compute_Linnaean_prob_country_gompertz(
  fit          = fit,
  year         = 2025,
  data         = data_L,
  include_data = F,
  re_formula   = NULL,
  draws        = NULL,
  return_draws = TRUE
)

prob_df <- cbind(data_L,lin_res$summary)

# ---- 1) species×continent weights (1/k) ----
w_cont <- base_df %>%
  distinct(valid_name, continent) %>%
  add_count(valid_name, name = "n_continent") %>%
  transmute(valid_name, continent, weight = 1 / n_continent)

continent_richness_weighted <- w_cont %>%
  group_by(continent) %>%
  summarise(
    richness_weighted = sum(weight),
    richness_weighted_round = round(richness_weighted),
    .groups = "drop"
  )

# ---- 2) continent posterior probabilities (weighted.mean) ----
continent_prob <- prob_df %>%
  select(valid_name, continent, prob_undesc_median, prob_undesc_lower, prob_undesc_upper) %>%
  left_join(w_cont, by = c("valid_name", "continent")) %>%
  filter(!is.na(weight)) %>%   # 
  group_by(continent) %>%
  summarise(
    prob_median  = weighted.mean(prob_undesc_median,  w = weight),
    prob_lower = weighted.mean(prob_undesc_lower, w = weight),
    prob_upper = weighted.mean(prob_undesc_upper, w = weight),
    .groups = "drop"
  )

# ---- 3) shortfall summary (SR = D/(1-p) - D) ----
eps <- 1e-12
continent_SR_summary <- continent_richness_weighted %>%
  left_join(continent_prob, by = "continent") %>%
  mutate(
    prob_median  = pmin(prob_median,  1 - eps),
    prob_lower = pmin(prob_lower, 1 - eps),
    prob_upper = pmin(prob_upper, 1 - eps),
    SRdesc    = round(pmax(richness_weighted / (1 - prob_median)  - richness_weighted, 0)),
    SRdesc_ll = round(pmax(richness_weighted / (1 - prob_lower) - richness_weighted, 0)),
    SRdesc_ul = round(pmax(richness_weighted / (1 - prob_upper) - richness_weighted, 0))
  )

# ---- 4) global ----
global_summary <- continent_SR_summary %>%
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

# ---- 5) final table ----
final_continent_global <- bind_rows(
  global_summary %>%
    transmute(region, richness_weighted, prob_median, prob_lower, prob_upper, SRdesc, SRdesc_ll, SRdesc_ul),
  continent_SR_summary %>%
    transmute(
      region = continent,
      richness_weighted, prob_median, prob_lower, prob_upper,
      SRdesc, SRdesc_ll, SRdesc_ul
    )
) %>%
  mutate(
    prob_display   = sprintf("%.3f (%.3f–%.3f)", prob_median, prob_lower, prob_upper),
    SRdesc_display = sprintf("%d (%d–%d)", SRdesc, SRdesc_ll, SRdesc_ul)
  ) %>%
  transmute(
    Region = region,
    richness_weighted,
    Undescribed_probability = prob_display,
    Undescribed_species     = SRdesc_display
  )

final_continent_global

write.csv(final_continent_global, paste0("output/tables/continent_linnaean_shortfall.csv"), row.names = FALSE)
rm(w_cont,continent_richness_weighted,continent_prob,continent_SR_summary,
   global_summary,final_continent_global)
################################################################################
# country-level probability of remaining undescribed (delay-derived)

# ---- 1) weights (species × country) + country richness ----
w_country <- base_df %>%
  distinct(valid_name, iso3) %>%
  add_count(valid_name, name = "n_iso3") %>%
  transmute(valid_name, iso3, weight = 1 / n_iso3)

country_richness_weighted <- w_country %>%
  group_by(iso3) %>%
  summarise(
    richness_weighted = sum(weight),
    richness_weighted_round = round(richness_weighted),
    .groups = "drop"
  )

# ---- 2) country posterior probabilities (weighted.mean) ----
country_prob <- prob_df %>%
  select(valid_name, iso3, prob_undesc_median, prob_undesc_lower, prob_undesc_upper) %>%
  left_join(w_country, by = c("valid_name", "iso3")) %>%
  filter(!is.na(weight)) %>%
  group_by(iso3) %>%
  summarise(
    prob_undesc_median  = weighted.mean(prob_undesc_median,  w = weight),
    prob_undesc_lower = weighted.mean(prob_undesc_lower, w = weight),
    prob_undesc_upper = weighted.mean(prob_undesc_upper, w = weight),
    .groups = "drop"
  )

# ---- 3) Join ----
country_dat <- country_richness_weighted %>%
  left_join(country_prob, by = "iso3")

# ---- 4) Compute SRdesc (your correct formula) ----
eps <- 1e-12
country_SR_summary <- country_dat %>%
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

country_SR_summary

write.csv(country_SR_summary, paste0("output/tables/country_linnaean_shortfall.csv"), row.names = FALSE)

######Drainage countrys with x number of undescribed species----------------------

calculate_species_percentages <- function(data) {
  # Load required libraries
  library(dplyr)
  
  # Remove NA values from SRdesc
  SRdesc_values <- data$SRdesc[!is.na(data$SRdesc)]
  
  # Calculate total number of drainage countrys (only non-NA SRdesc values are considered)
  total_country <- n_distinct(data$iso3[!is.na(data$SRdesc)])
  
  # Define fixed thresholds
  threshold_low <- 11   # Fixed lower threshold
  threshold_high <- 50  # Fixed upper threshold
  
  # Count the number of countrys below the lower threshold
  countrys_lt_threshold_low <- sum(data$SRdesc < threshold_low, na.rm = TRUE)
  
  # Count the number of countrys above the higher threshold
  countrys_gt_threshold_high <- sum(data$SRdesc > threshold_high, na.rm = TRUE)
  
  # Calculate the percentages based on total country count
  percentage_lt_threshold_low <- round((countrys_lt_threshold_low / total_country) * 100, 2)
  percentage_gt_threshold_high <- round((countrys_gt_threshold_high / total_country) * 100, 2)
  
  # Print results
  cat(sprintf("\n%% Drainage countrys with x < %.0f of undescribed species: %.2f%%",
              threshold_low, percentage_lt_threshold_low))
  cat(sprintf("\n%% Drainage countrys with x > %.0f of undescribed species: %.2f%%\n",
              threshold_high, percentage_gt_threshold_high))
}

country <- read.csv("input/raw/country_list.csv")
country_SR_summary %>% 
  left_join(country, by = "iso3") %>%
  filter(!is.na(continent)) %>%          
  group_by(continent) %>%
  group_walk(~ {
    cat("\n========== Realm:", .y$continent, "==========\n")
    calculate_species_percentages(.x)
  })

calculate_species_percentages(country_SR_summary)

rm(country_richness_weighted,country_prob,country_dat,country_SR_summary)