library(tibble)
library(scales)
library(ggplot2)
library(dplyr)        # Data manipulation verbs
library(brms)         # Bayesian multilevel models interface
library(arrow)

# ============================================================
# Load country-level data and construct species × continent table
# ============================================================
source("code/functions/xxx.r")
source("code/functions/gompertz_family.R")
data_W <- readRDS("output/stan_survdata_country_wallacean_single.rds")
fit <- readRDS("output/model/country_wallacean_gompertz_full.rds")
wal_res <- compute_Wallacean_prob(
  fit          = fit,
  year         = 2025,
  data         = data_W,
  year_var     = "year_description",
  include_data = F,
  re_formula   = NULL,
  draws        = NULL,
  return_draws = TRUE
)

wallace_df <- cbind(data_W,wal_res$summary)

################################################################################
# Compute weighted counts of non-geolocated species per continent
continent_wallace <- wallace_df %>%
  dplyr::filter(event == 0) %>%
  dplyr::distinct(valid_name, continent, .keep_all = TRUE) %>%
  dplyr::add_count(valid_name, name = "n_continent") %>%
  dplyr::mutate(weight = 1 / n_continent) %>%
  dplyr::group_by(continent) %>%
  dplyr::summarise(
    n_nongeoloc = sum(weight),
    prob_nongeoloc_median  = sum(weight * prob_nongeoloc_median)  / sum(weight),
    prob_nongeoloc_lower = sum(weight * prob_nongeoloc_lower) / sum(weight),
    prob_nongeoloc_upper = sum(weight * prob_nongeoloc_upper) / sum(weight),
    SRnoloc    = sum(weight * prob_nongeoloc_median),
    SRnoloc_ll = sum(weight * prob_nongeoloc_lower),
    SRnoloc_ul = sum(weight * prob_nongeoloc_upper),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    SRnoloc    = round(pmax(SRnoloc,    0)),
    SRnoloc_ll = round(pmax(SRnoloc_ll, 0)),
    SRnoloc_ul = round(pmax(SRnoloc_ul, 0))
  )

global_wallace <- wallace_df %>%
  dplyr::filter(event == 0) %>%  
  dplyr::distinct(valid_name, continent, .keep_all = TRUE) %>%
  dplyr::add_count(valid_name, name = "n_continent") %>%
  dplyr::mutate(weight = 1 / n_continent) %>%
  dplyr::summarise(
    n_nongeoloc = sum(weight),
    prob_nongeoloc_median  = sum(weight * prob_nongeoloc_median)  / sum(weight),
    prob_nongeoloc_lower = sum(weight * prob_nongeoloc_lower) / sum(weight),
    prob_nongeoloc_upper = sum(weight * prob_nongeoloc_upper) / sum(weight),
    SRnoloc     = sum(weight * prob_nongeoloc_median),
    SRnoloc_ll  = sum(weight * prob_nongeoloc_lower),
    SRnoloc_ul  = sum(weight * prob_nongeoloc_upper),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    continent      = "Global",
    SRnoloc    = round(pmax(SRnoloc,    0)),
    SRnoloc_ll = round(pmax(SRnoloc_ll, 0)),
    SRnoloc_ul = round(pmax(SRnoloc_ul, 0))
  )


final_continent_global <- dplyr::bind_rows(
  global_wallace,
  continent_wallace
) %>%
  mutate(
    prob_display   = sprintf("%.3f (%.3f–%.3f)", prob_nongeoloc_median,prob_nongeoloc_lower,prob_nongeoloc_upper),
    SRgeoloc_display = sprintf("%d (%d–%d)", SRnoloc, SRnoloc_ll, SRnoloc_ul)
  ) %>%
  select(
    Region = continent,
    n_nongeoloc = n_nongeoloc,
    Ungeolocated_probability = prob_display,
    Ungeolocated_species = SRgeoloc_display
  )


write.csv(final_continent_global, paste0("output/tables/continent_wallacean_shortfall.csv"), row.names = FALSE)
rm(continent_wallace,global_wallace,final_continent_global)
################################################################################
#country-level probability of remaining ungeolocated

country_wallace <- wallace_df %>%
  dplyr::filter(event == 0) %>%
  dplyr::distinct(valid_name, iso3, .keep_all = TRUE) %>%
  dplyr::add_count(valid_name, name = "n_iso3") %>%
  dplyr::mutate(weight = 1 / n_iso3) %>%
  dplyr::group_by(iso3) %>%
  dplyr::summarise(
    n_nongeoloc = sum(weight),
    prob_nongeoloc_median  = sum(weight * prob_nongeoloc_median)  / sum(weight),
    prob_nongeoloc_lower = sum(weight * prob_nongeoloc_lower) / sum(weight),
    prob_nongeoloc_upper = sum(weight * prob_nongeoloc_upper) / sum(weight),
    SRnoloc    = sum(weight * prob_nongeoloc_median),
    SRnoloc_ll = sum(weight * prob_nongeoloc_lower),
    SRnoloc_ul = sum(weight * prob_nongeoloc_upper),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    SRnoloc    = round(pmax(SRnoloc,    0)),
    SRnoloc_ll = round(pmax(SRnoloc_ll, 0)),
    SRnoloc_ul = round(pmax(SRnoloc_ul, 0))
  )


write.csv(country_wallace, paste0("output/tables/country_wallacean_shortfall.csv"), row.names = FALSE)

######Drainage countrys with x number of undescribed species----------------------

calculate_species_percentages <- function(data) {
  # Load required libraries
  library(dplyr)
  
  # Remove NA values from SRnoloc
  # SRnoloc_values <- data$SRnoloc[!is.na(data$SRnoloc)]
  
  # Calculate total number of drainage countrys (only non-NA SRnoloc values are considered)
  total_country_mean <- n_distinct(data$iso3[!is.na(data$SRnoloc)])
  total_country_ll <- n_distinct(data$iso3[!is.na(data$SRnoloc_ll)])
  total_country_ul <- n_distinct(data$iso3[!is.na(data$SRnoloc_ul)])
  
  # Define fixed thresholds
  threshold_low <- 11   # Fixed lower threshold
  threshold_high <- 50  # Fixed upper threshold
  
  # Count the number of countrys below the lower threshold
  countrys_lt_threshold_low_mean <- sum(data$SRnoloc < threshold_low, na.rm = TRUE)
  countrys_lt_threshold_low_ll <- sum(data$SRnoloc_ll < threshold_low, na.rm = TRUE)
  countrys_lt_threshold_low_ul <- sum(data$SRnoloc_ul < threshold_low, na.rm = TRUE)
  
  # Count the number of countrys above the higher threshold
  countrys_gt_threshold_high_mean <- sum(data$SRnoloc > threshold_high, na.rm = TRUE)
  countrys_gt_threshold_high_ll <- sum(data$SRnoloc_ll > threshold_high, na.rm = TRUE)
  countrys_gt_threshold_high_ul <- sum(data$SRnoloc_ul > threshold_high, na.rm = TRUE)
  
  # Calculate the percentages based on total country count
  percentage_lt_threshold_low_mean <- round((countrys_lt_threshold_low_mean / total_country_mean) * 100, 2)
  percentage_lt_threshold_low_ll <- round((countrys_lt_threshold_low_ll / total_country_ll) * 100, 2)
  percentage_lt_threshold_low_ul <- round((countrys_lt_threshold_low_ul / total_country_ul) * 100, 2)
  
  percentage_gt_threshold_high_mean <- round((countrys_gt_threshold_high_mean / total_country_mean) * 100, 2)
  percentage_gt_threshold_high_ll <- round((countrys_gt_threshold_high_ll / total_country_mean) * 100, 2)
  percentage_gt_threshold_high_ul <- round((countrys_gt_threshold_high_ul / total_country_mean) * 100, 2)
  
  # Print results
  cat(sprintf("\n%% Drainage countrys with x < %.0f of undescribed species(mean): %.2f%%",
              threshold_low, percentage_lt_threshold_low_mean))
  cat(sprintf("\n%% Drainage countrys with x < %.0f of undescribed species(ll): %.2f%%",
              threshold_low, percentage_lt_threshold_low_ll))
  cat(sprintf("\n%% Drainage countrys with x < %.0f of undescribed species(ul): %.2f%%\n",
              threshold_low, percentage_lt_threshold_low_ul))
  cat(sprintf("\n%% Drainage countrys with x > %.0f of undescribed species(mean): %.2f%%",
              threshold_high, percentage_gt_threshold_high_mean))
  cat(sprintf("\n%% Drainage countrys with x > %.0f of undescribed species(ll): %.2f%%",
              threshold_high, percentage_gt_threshold_high_ll))
  cat(sprintf("\n%% Drainage countrys with x > %.0f of undescribed species(ul): %.2f%%",
              threshold_high, percentage_gt_threshold_high_ul))
}

country <- read.csv("input/raw/country_list.csv")
country_wallace %>% 
  left_join(country, by = "iso3") %>%
  filter(!is.na(continent)) %>%          
  group_by(continent) %>%
  group_walk(~ {
    cat("\n========== continent:", .y$continent, "==========\n")
    calculate_species_percentages(.x)
  })

calculate_species_percentages(country_wallace)

rm(country_wallace,country)