library(tibble)
library(scales)
library(ggplot2)
library(dplyr)        # Data manipulation verbs
library(brms)         # Bayesian multilevel models interface
library(arrow)

# ============================================================
# Load basin-level data and construct species × realm table
# ============================================================
source("code/functions/xxx.r")
source("code/functions/gompertz_family.R")
fit_D <- readRDS("output/model/basin_darwinian_gompertz_full.rds")
data_D <- readRDS("output/stan_survdata_basin_darwinian_single.rds")
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
# Compute weighted counts of non-sequenced species per realm
realm_darwinian <- darwinian_df %>%
  dplyr::filter(event == 0) %>%          
  dplyr::distinct(valid_name, biogeographic_realm, .keep_all = TRUE) %>%
  dplyr::add_count(valid_name, name = "n_biogeographic_realm") %>%
  dplyr::mutate(weight = 1 / n_biogeographic_realm) %>%
  dplyr::group_by(biogeographic_realm) %>%
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
  dplyr::distinct(valid_name, biogeographic_realm, .keep_all = TRUE) %>%
  dplyr::add_count(valid_name, name = "n_biogeographic_realm") %>%
  dplyr::mutate(weight = 1 / n_biogeographic_realm) %>%
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
    biogeographic_realm      = "Global",
    SRnoseq    = round(pmax(SRnoseq,    0)),
    SRnoseq_ll = round(pmax(SRnoseq_ll, 0)),
    SRnoseq_ul = round(pmax(SRnoseq_ul, 0))
  )


final_realm_global <- dplyr::bind_rows(
  global_darwinian,
  realm_darwinian
) %>%
  mutate(
    prob_display   = sprintf("%.3f (%.3f–%.3f)", prob_noseq_median,prob_noseq_lower,prob_noseq_upper),
    SRnoseq_display = sprintf("%d (%d–%d)", SRnoseq, SRnoseq_ll, SRnoseq_ul)
  ) %>%
  select(
    Region = biogeographic_realm,
    n_noseq = n_noseq,
    Unsequenced_probability = prob_display,
    Unsequenced_species = SRnoseq_display
  )


write.csv(final_realm_global, paste0("output/tables/biogeographic_realm_darwinian_shortfall.csv"), row.names = FALSE)
rm(realm_darwinian,global_darwinian,final_realm_global)
################################################################################
#Basin-level probability of remaining unsequenced
basin_darwinian <- darwinian_df %>%
  dplyr::filter(event == 0) %>%          
  dplyr::distinct(valid_name, basin_id, .keep_all = TRUE) %>%
  dplyr::add_count(valid_name, name = "n_basin_id") %>%
  dplyr::mutate(weight = 1 / n_basin_id) %>%
  dplyr::group_by(basin_id) %>%
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

write.csv(basin_darwinian, paste0("output/tables/basin_darwinian_shortfall.csv"), row.names = FALSE)
######Drainage basins with x number of undescribed species----------------------

calculate_species_percentages <- function(data) {
  # Load required libraries
  library(dplyr)
  
  # Remove NA values from SRnoseq
  SRnoseq_values <- data$SRnoseq[!is.na(data$SRnoseq)]
  
  # Calculate total number of drainage basins (only non-NA SRnoseq values are considered)
  total_basin <- n_distinct(data$basin_id[!is.na(data$SRnoseq)])
  
  # Define fixed thresholds
  threshold_low <- 2   # Fixed lower threshold
  threshold_high <- 5  # Fixed upper threshold
  
  # Count the number of basins below the lower threshold
  basins_lt_threshold_low <- sum(data$SRnoseq < threshold_low, na.rm = TRUE)
  
  # Count the number of basins above the higher threshold
  basins_gt_threshold_high <- sum(data$SRnoseq > threshold_high, na.rm = TRUE)
  
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
basin_darwinian %>% 
  left_join(realm, by = "basin_id") %>%
  filter(!is.na(biogeographic_realm)) %>%          
  group_by(biogeographic_realm) %>%
  group_walk(~ {
    cat("\n========== Realm:", .y$biogeographic_realm, "==========\n")
    calculate_species_percentages(.x)
  })

calculate_species_percentages(basin_darwinian)


rm(basin_darwinian,realm)