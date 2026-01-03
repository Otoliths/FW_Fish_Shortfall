library(tibble)
library(scales)
library(ggplot2)
library(dplyr)        # Data manipulation verbs
library(brms)         # Bayesian multilevel models interface
library(arrow)

fit <- readRDS("output/model/basin_linnaean_all.rds") 
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
  mutate(shortfall = "Linnaean")

write.csv(model_results, paste0("output/tables/basin_linnaean_model_comparison.csv"), row.names = FALSE)
rm(comparison_loo,bayes_r2_df,bayes_r2_results,model_results)

################################################################################
# ------------------------------------------------------------------------------
# Estimating Linnaean shortfalls (number of undescribed species) at
# realm, basin and family levels
# ------------------------------------------------------------------------------
# This script quantifies Linnaean shortfalls by combining:
#   (i) observed described species richness at a given level L
#       (L ∈ {realm, basin, family}), and
#   (ii) posterior probabilities that species remain undescribed by a target year
#        (here: year = 2025) from a lognormal discovery-time model.
#
# For each level L, we define:
#   S_obs,L      = observed number of described species at level L
#                  (possibly rescaled so that the sum across levels
#                   matches the known global total of 18,821 species).
#   p_mean,L     = expected fraction of species that remain undescribed
#                  at level L, obtained as the richness-weighted mean
#                  posterior probability across species.
#   p_lower,L    = lower bound of this fraction (e.g. based on the
#                  mean of species-level lower posterior bounds).
#   p_upper,L    = upper bound of this fraction (analogously defined).
#
# Given these quantities, the total (described + undescribed) richness at
# level L is approximated by:
#
#       S_tot,L = S_obs,L / (1 - p_L)
#
# where p_L is p_mean,L, p_lower,L, or p_upper,L respectively.
#
# The expected number of undescribed species (Linnaean shortfall) at level L is:
#
#       SR_desc,L = S_tot,L - S_obs,L
#                 = S_obs,L * ( 1 / (1 - p_L) - 1 ).
#
# In practice, we compute:
#   SR_desc,L    using p_mean,L  (point estimate),
#   SR_desc_ll,L using p_lower,L (lower bound),
#   SR_desc_ul,L using p_upper,L (upper bound),
# and then truncate negative values to zero and round to integers.
#
# The same mathematical formulation is applied consistently at:
#   - realm level: L = biogeographic realm,
#   - basin level: L = drainage basin,
#   - family level: L = taxonomic family.
# ------------------------------------------------------------------------------

# ============================================================
# Load basin-level data and construct species × realm table
# ============================================================

# Basin-level dataset with at least:
# - valid_name           : species name (or other unique taxon ID)
# - biogeographic_realm  : realm code / factor
base_df <- read_parquet("input/data_prep/linnaean_basin.parquet")
# Species-level posterior probabilities
data_L <- readRDS("output/stan_survdata_basin_linnaean.rds")
source("code/functions/xxx.r")
lin_res <- compute_Linnaean_prob(
  fit          = fit$lognormal,
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
  select(valid_name, biogeographic_realm, prob_undesc_mean, prob_undesc_lower, prob_undesc_upper) %>%
  left_join(w_realm, by = c("valid_name", "biogeographic_realm")) %>%
  dplyr::filter(!is.na(weight)) %>%
  group_by(biogeographic_realm) %>%
  summarise(
    prob_mean  = weighted.mean(prob_undesc_mean,  w = weight),
    prob_lower = weighted.mean(prob_undesc_lower, w = weight),
    prob_upper = weighted.mean(prob_undesc_upper, w = weight),
    .groups = "drop"
  )

realm_dat <- realm_richness_weighted %>%
  left_join(realm_prob, by = "biogeographic_realm")

eps <- 1e-12
realm_SR_summary <- realm_dat %>%
  mutate(
    prob_mean  = pmax(pmin(prob_mean,  1 - eps), 0),
    prob_lower = pmax(pmin(prob_lower, 1 - eps), 0),
    prob_upper = pmax(pmin(prob_upper, 1 - eps), 0),
    SRdesc    = richness_weighted / (1 - prob_mean)  - richness_weighted,
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
    prob_mean  = SRdesc    / (richness_weighted + SRdesc),
    prob_lower = SRdesc_ll / (richness_weighted + SRdesc_ll),
    prob_upper = SRdesc_ul / (richness_weighted + SRdesc_ul)
  )

final_realm_global <- dplyr::bind_rows(
  global_summary %>%
    dplyr::select(region, richness_weighted, prob_mean, prob_lower, prob_upper, SRdesc, SRdesc_ll, SRdesc_ul),
  realm_SR_summary %>%
    dplyr::mutate(region = biogeographic_realm) %>%
    dplyr::select(region, richness_weighted, prob_mean, prob_lower, prob_upper, SRdesc, SRdesc_ll, SRdesc_ul)
) %>%
  mutate(
    prob_display   = sprintf("%.3f (%.3f–%.3f)", prob_mean, prob_lower, prob_upper),
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
  select(valid_name, basin_id, prob_undesc_mean, prob_undesc_lower, prob_undesc_upper) %>%
  left_join(w_basin, by = c("valid_name", "basin_id")) %>%
  dplyr::filter(!is.na(weight)) %>%
  group_by(basin_id) %>%
  summarise(
    prob_undesc_mean  = weighted.mean(prob_undesc_mean,  w = weight),
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
    prob_undesc_mean  = pmax(pmin(prob_undesc_mean,  1 - eps), 0),
    prob_undesc_lower = pmax(pmin(prob_undesc_lower, 1 - eps), 0),
    prob_undesc_upper = pmax(pmin(prob_undesc_upper, 1 - eps), 0),
    SRdesc_raw    = richness_weighted / (1 - prob_undesc_mean)  - richness_weighted,
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
