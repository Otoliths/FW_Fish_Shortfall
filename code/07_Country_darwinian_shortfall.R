library(tibble)
library(scales)
library(ggplot2)
library(dplyr)        # Data manipulation verbs
library(brms)         # Bayesian multilevel models interface
library(arrow)

################################################################################
# ============================================================
# Load country-level data and construct species × continent table
# ============================================================
source("code/functions/xxx.r")
source("code/functions/gompertz_family.R")
fit_D <- readRDS("output/model/country_darwinian_gompertz_full.rds")
data_D <- readRDS("output/stan_survdata_country_darwinian_single.rds")

dar_res <- compute_Darwinian_prob(
  fit          = fit_D,
  year         = 2025,
  data         = data_D,
  year_var     = "year_description",
  include_data = F,
  re_formula   = NULL,
  draws        = NULL,
  return_draws = TRUE
)

darwinian_df <- cbind(data_D,dar_res$summary)
################################################################################
# Compute weighted counts of non-sequenced species per continent
continent_darwinian <- darwinian_df %>%
  dplyr::filter(event == 0) %>%          
  dplyr::distinct(valid_name, continent, .keep_all = TRUE) %>%
  dplyr::add_count(valid_name, name = "n_continent") %>%
  dplyr::mutate(weight = 1 / n_continent) %>%
  dplyr::group_by(continent) %>%
  dplyr::summarise(
    n_noseq = sum(weight),
    prob_noseq_median  = sum(weight * prob_noseq_median)  / sum(weight),
    prob_noseq_lower = sum(weight * prob_noseq_lower) / sum(weight),
    prob_noseq_upper = sum(weight * prob_noseq_upper) / sum(weight),
    SRnoseq    = sum(weight * prob_noseq_median),
    SRnoseq_ll = sum(weight * prob_noseq_lower),
    SRnoseq_ul = sum(weight * prob_noseq_upper),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    SRnoseq    = round(pmax(SRnoseq,    0)),
    SRnoseq_ll = round(pmax(SRnoseq_ll, 0)),
    SRnoseq_ul = round(pmax(SRnoseq_ul, 0))
  )

global_darwinian <- darwinian_df %>%
  dplyr::filter(event == 0) %>%  
  dplyr::distinct(valid_name, continent, .keep_all = TRUE) %>%
  dplyr::add_count(valid_name, name = "n_continent") %>%
  dplyr::mutate(weight = 1 / n_continent) %>%
  dplyr::summarise(
    n_noseq = sum(weight),
    prob_noseq_median  = sum(weight * prob_noseq_median)  / sum(weight),
    prob_noseq_lower = sum(weight * prob_noseq_lower) / sum(weight),
    prob_noseq_upper = sum(weight * prob_noseq_upper) / sum(weight),
    SRnoseq     = sum(weight * prob_noseq_median),
    SRnoseq_ll  = sum(weight * prob_noseq_lower),
    SRnoseq_ul  = sum(weight * prob_noseq_upper),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    continent      = "Global",
    SRnoseq    = round(pmax(SRnoseq,    0)),
    SRnoseq_ll = round(pmax(SRnoseq_ll, 0)),
    SRnoseq_ul = round(pmax(SRnoseq_ul, 0))
  )


final_continent_global <- dplyr::bind_rows(
  global_darwinian,
  continent_darwinian
) %>%
  mutate(
    prob_display   = sprintf("%.3f (%.3f–%.3f)", prob_noseq_median,prob_noseq_lower,prob_noseq_upper),
    SRnoseq_display = sprintf("%d (%d–%d)", SRnoseq, SRnoseq_ll, SRnoseq_ul)
  ) %>%
  select(
    Region = continent,
    n_noseq = n_noseq,
    Unsequenced_probability = prob_display,
    Unsequenced_species = SRnoseq_display
  )


write.csv(final_continent_global, paste0("output/tables/continent_darwinian_shortfall.csv"), row.names = FALSE)
rm(continent_darwinian,global_darwinian,final_continent_global)
################################################################################
#country-level probability of remaining unsequenced
country_darwinian <- darwinian_df %>%
  dplyr::filter(event == 0) %>%          
  dplyr::distinct(valid_name, iso3, .keep_all = TRUE) %>%
  dplyr::add_count(valid_name, name = "n_iso3") %>%
  dplyr::mutate(weight = 1 / n_iso3) %>%
  dplyr::group_by(iso3) %>%
  dplyr::summarise(
    n_noseq = sum(weight),
    prob_noseq_median  = sum(weight * prob_noseq_median)  / sum(weight),
    prob_noseq_lower = sum(weight * prob_noseq_lower) / sum(weight),
    prob_noseq_upper = sum(weight * prob_noseq_upper) / sum(weight),
    SRnoseq    = sum(weight * prob_noseq_median),
    SRnoseq_ll = sum(weight * prob_noseq_lower),
    SRnoseq_ul = sum(weight * prob_noseq_upper),
    
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    SRnoseq    = round(pmax(SRnoseq,    0)),
    SRnoseq_ll = round(pmax(SRnoseq_ll, 0)),
    SRnoseq_ul = round(pmax(SRnoseq_ul, 0))
  )

write.csv(country_darwinian, paste0("output/tables/country_darwinian_shortfall.csv"), row.names = FALSE)
######Drainage country with x number of undescribed species----------------------

calculate_species_percentages <- function(data) {
  # Load required libraries
  library(dplyr)
  
  # Remove NA values from SRnoseq
  SRnoseq_values <- data$SRnoseq[!is.na(data$SRnoseq)]
  
  # Calculate total number of drainage country (only non-NA SRnoseq values are considered)
  total_country <- n_distinct(data$iso3[!is.na(data$SRnoseq)])
  
  # Define fixed thresholds
  threshold_low <- 11   # Fixed lower threshold
  threshold_high <- 50  # Fixed upper threshold
  
  # Count the number of country below the lower threshold
  country_lt_threshold_low <- sum(data$SRnoseq < threshold_low, na.rm = TRUE)
  
  # Count the number of country above the higher threshold
  country_gt_threshold_high <- sum(data$SRnoseq > threshold_high, na.rm = TRUE)
  
  # Calculate the percentages based on total country count
  percentage_lt_threshold_low <- round((country_lt_threshold_low / total_country) * 100, 2)
  percentage_gt_threshold_high <- round((country_gt_threshold_high / total_country) * 100, 2)
  
  # Print results
  cat(sprintf("\n%% Drainage country with x < %.0f of undescribed species: %.2f%%",
              threshold_low, percentage_lt_threshold_low))
  cat(sprintf("\n%% Drainage country with x > %.0f of undescribed species: %.2f%%\n",
              threshold_high, percentage_gt_threshold_high))
}
country <- read.csv("input/raw/country_list.csv")
country_darwinian %>% 
  left_join(country, by = "iso3") %>%
  filter(!is.na(continent)) %>%          
  group_by(continent) %>%
  group_walk(~ {
    cat("\n========== continent:", .y$continent, "==========\n")
    calculate_species_percentages(.x)
  })

calculate_species_percentages(country_darwinian)


rm(country_darwinian,country)