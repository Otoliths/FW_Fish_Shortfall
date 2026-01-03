library(tibble)
library(scales)
library(ggplot2)
library(dplyr)        # Data manipulation verbs
library(brms)         # Bayesian multilevel models interface
library(arrow)

fit <- readRDS("output/model/basin_wallacean_all.rds") 
comparison_loo <- loo_compare(fit$lognormal, 
                              fit$gamma,
                              fit$exponential,
                              fit$weibull,
                              model_names = names(fit), 
                              criterion = "loo") %>%
  as.data.frame()%>%
  tibble::rownames_to_column(var = "model")
# comparison_waic <- loo_compare(fit$lognormal, 
#                                fit$gamma,
#                                fit$exponential,
#                                fit$weibull,
#                                model_names = names(fit), 
#                                criterion = "waic")%>%
#   as.data.frame()%>%
#   tibble::rownames_to_column(var = "model")

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
  mutate(shortfall = "Wallacean")

write.csv(model_results, paste0("output/tables/basin_wallacean_model_comparison.csv"), row.names = FALSE)
rm(comparison_loo,bayes_r2_df,bayes_r2_results,model_results)

################################################################################
# ------------------------------------------------------------------------------
# Estimating Wallacean shortfalls (number of species lacking georeferenced
# occurrence records) at realm, basin, and family levels
# ------------------------------------------------------------------------------
# Wallacean shortfalls quantify the expected number of already-described species
# that will still lack a georeferenced occurrence record by a target year
# (here, 2025). Unlike Linnaean shortfalls, this metric does not require
# estimating total species richness, because the analysis is restricted to
# species that are already formally described.
#
# For each species i, a lognormal time-to-event model provides the posterior
# probability:
#
#       p_i = Pr(T_i > T*, the species has not yet obtained its first
#             georeferenced occurrence by year T* = 2025).
#
# Let U denote the set of currently described species that lack any
# georeferenced occurrence record:
#
#       U = { i : species i has no geolocation record at present }.
#
# For any taxonomic or geographic level L (realm, basin, or family), let:
#
#       n_L = | U_L |  = number of described species in level L that currently
#                        lack a georeferenced occurrence.
#
#       p_mean,L  = mean posterior probability (over i ∈ U_L) that these species
#                    will still lack a georeferenced record by 2025.
#
#       p_lower,L , p_upper,L = corresponding uncertainty bounds.
#
# The expected Wallacean shortfall at level L is therefore:
#
#       WS_L      = sum_{i ∈ U_L} p_i
#                 ≈ n_L * p_mean,L
#
# with uncertainty bounds:
#
#       WS_ll,L   = n_L * p_lower,L
#       WS_ul,L   = n_L * p_upper,L
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
data_W <- readRDS("output/stan_survdata_basin_wallacean.rds")
wal_res <- compute_Wallacean_prob(
  fit          = fit$lognormal,
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
# Compute weighted counts of non-geolocated species per realm
realm_wallace <- wallace_df %>%
  dplyr::filter(event == 0) %>%
  dplyr::distinct(valid_name, biogeographic_realm, .keep_all = TRUE) %>%
    dplyr::add_count(valid_name, name = "n_biogeographic_realm") %>%
  dplyr::mutate(weight = 1 / n_biogeographic_realm) %>%
  dplyr::group_by(biogeographic_realm) %>%
  dplyr::summarise(
    n_nongeoloc = sum(weight),
    prob_nongeoloc_mean  = sum(weight * prob_nongeoloc_mean)  / sum(weight),
    prob_nongeoloc_lower = sum(weight * prob_nongeoloc_lower) / sum(weight),
    prob_nongeoloc_upper = sum(weight * prob_nongeoloc_upper) / sum(weight),
    SRnoloc    = sum(weight * prob_nongeoloc_mean),
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
  dplyr::distinct(valid_name, biogeographic_realm, .keep_all = TRUE) %>%
  dplyr::add_count(valid_name, name = "n_biogeographic_realm") %>%
  dplyr::mutate(weight = 1 / n_biogeographic_realm) %>%
  dplyr::summarise(
    n_nongeoloc = sum(weight),
    prob_nongeoloc_mean  = sum(weight * prob_nongeoloc_mean)  / sum(weight),
    prob_nongeoloc_lower = sum(weight * prob_nongeoloc_lower) / sum(weight),
    prob_nongeoloc_upper = sum(weight * prob_nongeoloc_upper) / sum(weight),
    SRnoloc     = sum(weight * prob_nongeoloc_mean),
    SRnoloc_ll  = sum(weight * prob_nongeoloc_lower),
    SRnoloc_ul  = sum(weight * prob_nongeoloc_upper),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    biogeographic_realm      = "Global",
    SRnoloc    = round(pmax(SRnoloc,    0)),
    SRnoloc_ll = round(pmax(SRnoloc_ll, 0)),
    SRnoloc_ul = round(pmax(SRnoloc_ul, 0))
  )


final_realm_global <- dplyr::bind_rows(
  global_wallace,
  realm_wallace
) %>%
  mutate(
    prob_display   = sprintf("%.3f (%.3f–%.3f)", prob_nongeoloc_mean,prob_nongeoloc_lower,prob_nongeoloc_upper),
    SRgeoloc_display = sprintf("%d (%d–%d)", SRnoloc, SRnoloc_ll, SRnoloc_ul)
  ) %>%
  select(
    Region = biogeographic_realm,
    n_nongeoloc = n_nongeoloc,
    Ungeolocated_probability = prob_display,
    Ungeolocated_species = SRgeoloc_display
  )


write.csv(final_realm_global, paste0("output/tables/biogeographic_realm_wallacean_shortfall.csv"), row.names = FALSE)
rm(realm_wallace,global_wallace,final_realm_global)
################################################################################
#Basin-level probability of remaining ungeolocated

basin_wallace <- wallace_df %>%
  dplyr::filter(event == 0) %>%
  dplyr::distinct(valid_name, basin_id, .keep_all = TRUE) %>%
  dplyr::add_count(valid_name, name = "n_basin_id") %>%
  dplyr::mutate(weight = 1 / n_basin_id) %>%
  dplyr::group_by(basin_id) %>%
  dplyr::summarise(
    n_nongeoloc = sum(weight),
    prob_nongeoloc_mean  = sum(weight * prob_nongeoloc_mean)  / sum(weight),
    prob_nongeoloc_lower = sum(weight * prob_nongeoloc_lower) / sum(weight),
    prob_nongeoloc_upper = sum(weight * prob_nongeoloc_upper) / sum(weight),
    SRnoloc    = sum(weight * prob_nongeoloc_mean),
    SRnoloc_ll = sum(weight * prob_nongeoloc_lower),
    SRnoloc_ul = sum(weight * prob_nongeoloc_upper),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    SRnoloc    = round(pmax(SRnoloc,    0)),
    SRnoloc_ll = round(pmax(SRnoloc_ll, 0)),
    SRnoloc_ul = round(pmax(SRnoloc_ul, 0))
  )


write.csv(basin_wallace, paste0("output/tables/basin_wallacean_shortfall.csv"), row.names = FALSE)

######Drainage basins with x number of undescribed species----------------------

calculate_species_percentages <- function(data) {
  # Load required libraries
  library(dplyr)
  
  # Remove NA values from SRnoloc
  # SRnoloc_values <- data$SRnoloc[!is.na(data$SRnoloc)]
  
  # Calculate total number of drainage basins (only non-NA SRnoloc values are considered)
  total_basin_mean <- n_distinct(data$basin_id[!is.na(data$SRnoloc)])
  total_basin_ll <- n_distinct(data$basin_id[!is.na(data$SRnoloc_ll)])
  total_basin_ul <- n_distinct(data$basin_id[!is.na(data$SRnoloc_ul)])
  
  # Define fixed thresholds
  threshold_low <- 2   # Fixed lower threshold
  threshold_high <- 5  # Fixed upper threshold
  
  # Count the number of basins below the lower threshold
  basins_lt_threshold_low_mean <- sum(data$SRnoloc < threshold_low, na.rm = TRUE)
  basins_lt_threshold_low_ll <- sum(data$SRnoloc_ll < threshold_low, na.rm = TRUE)
  basins_lt_threshold_low_ul <- sum(data$SRnoloc_ul < threshold_low, na.rm = TRUE)
  
  # Count the number of basins above the higher threshold
  basins_gt_threshold_high_mean <- sum(data$SRnoloc > threshold_high, na.rm = TRUE)
  basins_gt_threshold_high_ll <- sum(data$SRnoloc_ll > threshold_high, na.rm = TRUE)
  basins_gt_threshold_high_ul <- sum(data$SRnoloc_ul > threshold_high, na.rm = TRUE)
  
  # Calculate the percentages based on total basin count
  percentage_lt_threshold_low_mean <- round((basins_lt_threshold_low_mean / total_basin_mean) * 100, 2)
  percentage_lt_threshold_low_ll <- round((basins_lt_threshold_low_ll / total_basin_ll) * 100, 2)
  percentage_lt_threshold_low_ul <- round((basins_lt_threshold_low_ul / total_basin_ul) * 100, 2)
  
  percentage_gt_threshold_high_mean <- round((basins_gt_threshold_high_mean / total_basin_mean) * 100, 2)
  percentage_gt_threshold_high_ll <- round((basins_gt_threshold_high_ll / total_basin_mean) * 100, 2)
  percentage_gt_threshold_high_ul <- round((basins_gt_threshold_high_ul / total_basin_mean) * 100, 2)
  
  # Print results
  cat(sprintf("\n%% Drainage basins with x < %.0f of undescribed species(mean): %.2f%%",
              threshold_low, percentage_lt_threshold_low_mean))
  cat(sprintf("\n%% Drainage basins with x < %.0f of undescribed species(ll): %.2f%%",
              threshold_low, percentage_lt_threshold_low_ll))
  cat(sprintf("\n%% Drainage basins with x < %.0f of undescribed species(ul): %.2f%%\n",
              threshold_low, percentage_lt_threshold_low_ul))
  cat(sprintf("\n%% Drainage basins with x > %.0f of undescribed species(mean): %.2f%%",
              threshold_high, percentage_gt_threshold_high_mean))
  cat(sprintf("\n%% Drainage basins with x > %.0f of undescribed species(ll): %.2f%%",
              threshold_high, percentage_gt_threshold_high_ll))
  cat(sprintf("\n%% Drainage basins with x > %.0f of undescribed species(ul): %.2f%%",
              threshold_high, percentage_gt_threshold_high_ul))
}

realm <- read.csv("input/raw/biogeographic_list.csv")
basin_wallace %>% 
  left_join(realm, by = "basin_id") %>%
  filter(!is.na(biogeographic_realm)) %>%          
  group_by(biogeographic_realm) %>%
  group_walk(~ {
    cat("\n========== Realm:", .y$biogeographic_realm, "==========\n")
    calculate_species_percentages(.x)
  })

calculate_species_percentages(basin_wallace)

rm(basin_wallace,realm)
