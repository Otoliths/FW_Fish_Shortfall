library(tibble)
library(scales)
library(ggplot2)
library(dplyr)        # Data manipulation verbs
library(brms)         # Bayesian multilevel models interface
library(arrow)

fit <- readRDS("output/model/basin_darwinian_all.rds") 
comparison_loo <- loo_compare(fit$lognormal, 
                              fit$gamma,
                              fit$exponential,
                              fit$weibull,
                              model_names = names(fit), 
                              criterion = "loo") %>%
  as.data.frame()%>%
  tibble::rownames_to_column(var = "model")
comparison_waic <- loo_compare(fit$lognormal,
                               fit$gamma,
                               fit$exponential,
                               fit$weibull,
                               model_names = names(fit),
                               criterion = "waic")%>%
  as.data.frame()%>%
  tibble::rownames_to_column(var = "model")

bayes_r2_results <- lapply(fit, bayes_R2)
bayes_r2_df <- do.call(rbind, lapply(names(bayes_r2_results), function(model) {
  r2_data <- bayes_r2_results[[model]]
  data.frame(
    model = model,
    bayes_r2_estimate = r2_data["R2", "Estimate"],
    bayes_r2_se = r2_data["R2", "Est.Error"],
    bayes_r2_2.5 = r2_data["R2", "Q2.5"],
    bayes_r2_97.5 = r2_data["R2", "Q97.5"]
  )
}))

model_results <- comparison_loo %>%
  left_join(bayes_r2_df, by = "model") %>%
  mutate(shortfall = "Darwinian")

write.csv(model_results, paste0("output/tables/basin_darwinian_model_comparison.csv"), row.names = FALSE)
rm(comparison_loo,bayes_r2_df,bayes_r2_results,model_results)

################################################################################
# ------------------------------------------------------------------------------
# Estimating Darwinian shortfalls (number of species lacking genetic sequences)
# at realm, basin, and family levels
# ------------------------------------------------------------------------------
# Darwinian shortfalls quantify the expected number of already-described species
# that will still lack any genetic sequence (e.g., COI, Cytb, 16S) by a target
# year (here, 2025). Unlike Linnaean shortfalls, this metric does not require
# estimating total species richness, because the analysis is restricted to
# species that are already formally described.
#
# For each species i, a lognormal time-to-event model provides the posterior
# probability:
#
#       p_i = Pr(T_i > T*, meaning the species has not yet obtained its first
#             genetic sequence by year T* = 2025).
#
# Let U denote the set of currently described species that lack any genetic
# sequence at present:
#
#       U = { i : species i has no available genetic sequence }.
#
# For any taxonomic or geographic level L (realm, basin, or family), define:
#
#       n_L = |U_L|  = number of described species in level L that currently
#                      lack a genetic sequence.
#
#       p_mean,L   = mean posterior probability (over i ∈ U_L) that these
#                     species will still remain unsequenced by 2025.
#
#       p_lower,L , p_upper,L = corresponding uncertainty bounds.
#
# The expected Darwinian shortfall at level L is:
#
#       DS_L      = sum_{i ∈ U_L} p_i
#                 ≈ n_L * p_mean,L
#
# with uncertainty intervals:
#
#       DS_ll,L   = n_L * p_lower,L
#       DS_ul,L   = n_L * p_upper,L
#
# Because negative values are impossible, all estimates are truncated at zero
# and rounded to integers. This formulation applies consistently to:
#   - biogeographic realms,
#   - drainage basins,
#   - taxonomic families.
# ------------------------------------------------------------------------------

# ============================================================
# Load basin-level data and construct species × realm table
# ============================================================
source("code/functions/xxx.r")
fit_D <- fit$gamma
data_D <- readRDS("output/stan_survdata_basin_darwinian.rds")
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
    prob_noseq_mean  = sum(weight * prob_noseq_mean)  / sum(weight),
    prob_noseq_lower = sum(weight * prob_noseq_lower) / sum(weight),
    prob_noseq_upper = sum(weight * prob_noseq_upper) / sum(weight),
    SRnoseq    = sum(weight * prob_noseq_mean),
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
    prob_noseq_mean  = sum(weight * prob_noseq_mean)  / sum(weight),
    prob_noseq_lower = sum(weight * prob_noseq_lower) / sum(weight),
    prob_noseq_upper = sum(weight * prob_noseq_upper) / sum(weight),
    SRnoseq     = sum(weight * prob_noseq_mean),
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
    prob_display   = sprintf("%.3f (%.3f–%.3f)", prob_noseq_mean,prob_noseq_lower,prob_noseq_upper),
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
    prob_noseq_mean  = sum(weight * prob_noseq_mean)  / sum(weight),
    prob_noseq_lower = sum(weight * prob_noseq_lower) / sum(weight),
    prob_noseq_upper = sum(weight * prob_noseq_upper) / sum(weight),
    SRnoseq    = sum(weight * prob_noseq_mean),
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
