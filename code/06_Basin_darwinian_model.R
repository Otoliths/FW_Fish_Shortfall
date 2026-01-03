# ------------------------------------------------------------
# Load required libraries for data processing and modeling
# ------------------------------------------------------------
library(dplyr)        # Data manipulation verbs
library(magrittr)     # Pipe operator support (%>% and .)
library(data.table)   # Fast data handling for large tables
library(rstanarm)     # Bayesian regression modeling (not used directly here, but loaded)
library(brms)         # Interface for Bayesian multilevel models using Stan
library(scales)       # Scaling helper functions (e.g., rescale)
library(bestNormalize) # Normalization / transformation utilities (not used directly here)
library(arrow)        # Read/write parquet and other columnar data formats

# Use cmdstanr backend for brms and enable debugging output
options(brms.backend = "cmdstanr", brms.debug = TRUE)

# ------------------------------------------------------------
# 1. Load input data
# ------------------------------------------------------------

# Import raw species data from Excel file (CAS freshwater list)
cas <- openxlsx::read.xlsx("input/raw/cas_freshwater_v1.xlsx")

# Load pre-processed basin-level data from parquet
base_df <- read_parquet("input/data_prep/darwinian_basin.parquet")

# Merge CAS species attributes (e.g. valid_name, year of description) into basin data
# Column indices 5,2,and 8 in 'cas' should correspond to valid_name ,year_description, and family
data <- base_df %>% 
  left_join(cas[, c(5,2,8)], by = "valid_name")


basin_discharge_avg <- read_parquet("input/data_prep/basin_discharge_avg.parquet")
basin_watertemp_avg <- read_parquet("input/data_prep/basin_watertemp_avg.parquet")
basin_rarity_avg <- read_parquet("input/data_prep/basin_rarity_avg.parquet")
basin_population_density_avg <- read_parquet("input/data_prep/basin_population_density_avg.parquet")
basin_elevation_avg <- read_parquet("input/data_prep/basin_elevation_avg.parquet")
basin_latitude_avg <- read_parquet("input/data_prep/basin_latitude_avg.parquet")

data_filled <- data %>%
  # ---- 1) discharge ----
left_join(
  basin_discharge_avg %>%
    select(basin_id, discharge_avg = discharge),
  by = c("basin_id")
) %>%
  mutate(
    discharge = if_else(
      is.na(discharge) & is.na(year_sequence),
      discharge_avg,
      discharge
    )
  ) %>%
  select(-discharge_avg) %>%
  
  # ---- 2) watertemp ----
left_join(
  basin_watertemp_avg %>%
    select(basin_id, watertemp_avg = watertemp),
  by = c("basin_id")
) %>%
  mutate(
    watertemp = if_else(
      is.na(watertemp) & is.na(year_sequence),
      watertemp_avg,
      watertemp
    )
  ) %>%
  select(-watertemp_avg) %>%
  
  # ---- 3) range_rarity ----
left_join(
  basin_rarity_avg %>%
    select(basin_id, valid_name,range_rarity_avg = range_rarity_mean),
  by = c("basin_id","valid_name")
) %>%
  mutate(
    range_rarity = if_else(
      is.na(range_rarity) & is.na(year_sequence),
      range_rarity_avg,
      range_rarity
    )
  ) %>%
  select(-range_rarity_avg) %>%
  
  # ---- 4) population_density ----
left_join(
  basin_population_density_avg %>%
    select(basin_id, population_density_avg = population_density),
  by = c("basin_id")
) %>%
  mutate(
    population_density = if_else(
      is.na(population_density) & is.na(year_sequence),
      population_density_avg,
      population_density
    )
  ) %>%
  select(-population_density_avg)  %>%
  # ---- 5) elevation ----
left_join(
  basin_elevation_avg %>%
    select(basin_id, elevation_avg = elevation) %>%
    filter(!is.na(elevation_avg)),
  by = c("basin_id")
) %>%
  mutate(
    elevation = if_else(
      is.na(elevation) & is.na(year_sequence),
      elevation_avg,
      elevation
    )
  ) %>%
  select(-elevation_avg) %>%
  # ---- 6) latitude ----
left_join(
  basin_latitude_avg %>%
    select(basin_id, latitude_avg = latitude) %>%
    filter(!is.na(latitude_avg)),
  by = c("basin_id")
) %>%
  mutate(
    latitude = if_else(
      is.na(latitude) & is.na(year_sequence),
      latitude_avg,
      latitude
    )
  ) %>%
  select(-latitude_avg)

data_filled$range_size[is.na(data_filled$range_size)] <- 0
data_filled$preserved_specimen[is.na(data_filled$preserved_specimen)]  <- 0
data_filled$sequencing_effort[is.na(data_filled$sequencing_effort)]  <- 0
data_filled$sampling_effort[is.na(data_filled$sampling_effort)]  <- 0

rm(basin_discharge_avg,basin_watertemp_avg,basin_rarity_avg,basin_elevation_avg,
   basin_population_density_avg,cas,data,base_df,basin_latitude_avg)
# ------------------------------------------------------------
# 2. Basic filtering and construction of survival time
# ------------------------------------------------------------

# We now treat the availability of a molecular sequence as a time-to-event
# process at the species level.
#
# - year_description: year when the species was formally described
# - year_sequence   : earliest year when a molecular sequence (e.g. DNA barcode
#                     or genomic record) became available for that species
#
# For species with at least one sequence record (year_sequence not NA),
#   time  = year_sequence - year_description
#   event = 1  (the "sequencing event" has occurred)
#
# For species without any sequence record (year_sequence is NA),
#   we treat them as right-censored at a fixed cutoff year (e.g. 2024):
#   time  = cutoff_year - year_description
#   event = 0  (no sequencing event observed up to cutoff_year)
#
# This definition allows us to fit a standard survival model with both
# events (sequenced species) and right-censored observations (unsequenced
# species), quantifying the Darwinian shortfall in a coherent framework.

cutoff_year <- 2024  # analysis cutoff; adjust if needed

data_surv <- data_filled %>%
  # Keep species with a valid description year
  filter(!is.na(year_description)) %>%
  mutate(
    # Indicator for whether the species has at least one sequence record
    has_sequence = !is.na(year_sequence),
    
    # Survival time in years:
    # - if sequenced: lag between description and first sequence
    # - if unsequenced: time from description to cutoff_year (right-censored)
    time = if_else(
      has_sequence,
      year_sequence - year_description,
      cutoff_year   - year_description
    ),
    
    # Event indicator:
    # 1 = sequencing event observed, 0 = right-censored at cutoff_year
    event = if_else(has_sequence, 1L, 0L)
  ) %>%
  # Remove species with non-positive or pathological times
  filter(time > 0)

vars <- c(
  "log_length_max", "taxonmic_effort", "taxonomic_activity",
  "watershed_area", "range_size", "elevation", "latitude",
  "discharge", "watertemp", "preserved_specimen",
  "sequencing_effort", "sampling_effort",
  "range_rarity", "population_density"
)


data_surv <- data_surv[, c("basin_id",
                           "valid_name",
                   "biogeographic_realm",
                   vars,
                   "family",
                   "event","time","year_description")] %>% .[complete.cases(.), ] 


table(data_surv$event)
# 0     1 
# 84871 10585 
rm(data_filled)
# ------------------------------------------------------------
# 3. Define predictor set and transformation helpers
# ------------------------------------------------------------

# VIF < 5 is typically considered acceptable
usdm::vif(data_surv[, vars] %>% as.data.frame())

# ------------------------------------------------------------
# 4. Apply transformations: create "z_*" columns
# ------------------------------------------------------------
df_transformed <- data_surv %>%
  mutate(across(all_of(vars),
                ~ as.numeric(scale(.x)),
                .names = "z_{.col}"))


usdm::vif(df_transformed[,22:35] %>% as.data.frame())
# ------------------------------------------------------------
# 5. Assemble survival modeling dataset
# ------------------------------------------------------------

# Keep only the variables required for the survival analysis:
# - basin_id: drainage basin identifier (nested within biogeographic realms)
# - biogeographic_realm: one of the six major global freshwater biogeographic regions
# - names(preds): transformed predictor variables (all "log_*" columns)
# - family: taxonomic family used as an additional random-effect grouping factor
# - time: sequence lag time (years between species description and its first sequenced record)
# - event: 1 = sequenced, 0 = not yet sequenced by cutoff_year
stan_survdata <- df_transformed %>%
  dplyr::transmute(
    basin_id            = factor(basin_id),
    valid_name          = valid_name,
    biogeographic_realm = factor(biogeographic_realm),
    
    body_size           = z_log_length_max,
    taxonmic_effort     = z_taxonmic_effort,
    taxonomic_activity  = z_taxonomic_activity,
    watershed_area      = z_watershed_area,
    range_size          = z_range_size,
    elevation           = z_elevation,
    latitude            = z_latitude,
    discharge           = z_discharge,
    watertemp           = z_watertemp,
    preserved_specimen  = z_preserved_specimen,
    sequencing_effort   = z_sequencing_effort,
    sampling_effort     = z_sampling_effort,
    range_rarity        = z_range_rarity,
    population_density  = z_population_density,
    
    family_group        = factor(family),
    time                = time,
    event               = as.integer(event),
    year_description    = year_description
  ) %>%
  .[complete.cases(.), ] %>%
  as.data.frame()%>%
  distinct()

rm(data_surv,df_transformed);gc()
#saveRDS(stan_survdata,"stan_survdata_basin_darwinian.rds")
# ------------------------------------------------------------
# 6. Specify Bayesian survival model formula
# ------------------------------------------------------------

# Darwinian sequencing model with global + realm-level slopes
# - time: years from species description to first sequence (or censoring)
# - event: 1 = sequenced, 0 = not yet sequenced by cutoff_year
# - cens(1 - event): 0 = exact event time, 1 = right-censored

# Vector of 13 standardized predictors
preds <- c(
  "body_size",
  "taxonmic_effort",
  "taxonomic_activity",
  "watershed_area",
  "range_size",
  "elevation",
  "latitude",
  "discharge",
  "watertemp",
  "preserved_specimen",
  #"sequencing_effort", #removed sequencing effort avoid circularity in Darwinian modeling
  "sampling_effort", 
  "range_rarity",
  "population_density"
)

# Collapse predictors into a single string "x1 + x2 + ... + x14"
pred_str <- paste(preds, collapse = " + ")

# ------------------------------------------------------------
# Hierarchical survival model
#
#  • Global fixed effects for all 13 predictors
#  • Realm-level random slopes for the same 13 predictors
#    (partial pooling avoids numerical instability in data-poor realms such as Oceania)
#  • Basin- and family-level random intercepts
# ------------------------------------------------------------

form <- bf(
  as.formula(
    paste(
      # survival likelihood with censoring indicator
      "time | cens(1 - event) ~ 1 +",
      
      # fixed-effect interactions: all predictors × realm
      paste0("(", pred_str, ") * biogeographic_realm"),
      
      # hierarchical random intercepts
      "+ (1 | basin_id)",
      "+ (1 | family_group)"
    )
  )
)


# ------------------------------------------------------------
# 7. Global MCMC configuration (shared across distributions)
# ------------------------------------------------------------

priors <- c(
  prior(normal(0, 0.5), class = "b"),          # Global slopes and all other fixed effects
  prior(normal(0, 0.5), class = "Intercept"),  # Global intercept
  prior(student_t(3, 0, 1), class = "sd")      # SDs of realm/basin/family random effects
)




common_args <- list(
  formula = form,
  data    = stan_survdata,
  chains  = 8,
  cores   = 8,
  seed    = 2025,
  backend = "cmdstanr",
  stan_model_args = list(
    cpp_options = list(
      "STAN_OPENCL=TRUE",           # Enable OpenCL (GPU) support if available
      "OPENCL_DEVICE=0",            # Select OpenCL device
      "OPENCL_PLATFORM=0",          # Select OpenCL platform
      "CXXFLAGS += -O3 -Wno-overloaded-virtual"  # Compiler optimization flags
    )
  ),
  threads   = threading(10, static = TRUE),  # Within-chain parallelization
  normalize = FALSE,                         # Use raw scale (we already normalized predictors)
  thin      = 10,                            # Thinning to reduce autocorrelation / disk usage
  save_pars = save_pars(all = TRUE),         # Save all parameters (useful for post-processing)
  init      = 0.1,                          # Small positive initial values
  control = list(adapt_delta = 0.999,         # High target acceptance to reduce divergences
                 max_treedepth = 15          # ↑ Allow deeper NUTS trees; avoids hitting limit at 10
  ),
  prior     = priors 
)

# ------------------------------------------------------------
# 8. Model configurations: candidate survival distributions
# ------------------------------------------------------------

# Each list element specifies a different parametric survival family
# and the corresponding number of iterations / warmup samples.
model_configs <- list(
  list(model_name = "lognormal",   family = brmsfamily("lognormal"),   iter = 4000, warmup = 2000),
  list(model_name = "gamma",       family = brmsfamily("Gamma"),      iter = 4000, warmup = 2000),
  list(model_name = "exponential", family = brmsfamily("exponential"), iter = 4000, warmup = 2000),
  list(model_name = "weibull",     family = brmsfamily("weibull"),    iter = 4000, warmup = 2000)
)


# ------------------------------------------------------------
# 9. Loop over model configurations and fit each distribution
# ------------------------------------------------------------

# Create output directory (for per-model and combined saves)
output_dir <- "output/model"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

fits <- list()

for (config in model_configs) {
  model_name  <- config$model_name
  family_name <- config$family$family
  
  # Build file paths for raw and final model objects
  raw_path   <- file.path(output_dir, sprintf("basin_darwinian_%s_raw.rds",  model_name))
  final_path <- file.path(output_dir, sprintf("basin_darwinian_%s.rds",      model_name))
  
  tryCatch({
    # --------------------------------------------------------
    # 0) Check if model files already exist (skip re-fitting)
    # --------------------------------------------------------
    
    if (file.exists(final_path)) {
      message(sprintf(
        "\n[SKIP] %s model (%s) already fitted. Loading from %s",
        toupper(model_name),
        family_name,
        final_path
      ))
      
      fit_loaded <- readRDS(final_path)
      fits[[model_name]] <- fit_loaded
      next
    }
    
    # If the final model does not exist, check for a raw version
    if (file.exists(raw_path)) {
      message(sprintf(
        "\n[RESUME] Found raw %s model (%s). Loading from %s and adding criteria...",
        toupper(model_name),
        family_name,
        raw_path
      ))
      fit_tmp <- readRDS(raw_path)
    } else {
      # --------------------------------------------------------
      # 1) Fit a new Bayesian model from scratch
      # --------------------------------------------------------
      start_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      message(sprintf(
        "\n[FITTING] %s model (%s family) starting at %s",
        toupper(model_name),
        family_name,
        start_time
      ))
      
      # Combine global MCMC settings with model-specific arguments
      args <- modifyList(common_args, config[!names(config) %in% "model_name"])
      
      fit_tmp <- do.call(brm, args)
      
      # Save the raw model
      saveRDS(fit_tmp, raw_path)
      if (file.exists(raw_path)) {
        message(sprintf("[OK] Raw model saved to %s", raw_path))
      } else {
        warning(sprintf("[WARNING] Raw model for %s not found after saveRDS.", model_name))
      }
    }
    
    # Store the raw model (either newly fitted or loaded)
    fits[[model_name]] <- fit_tmp
    
    # --------------------------------------------------------
    # 2) Try to add LOO/WAIC criteria
    #    If this fails, keep the raw model
    # --------------------------------------------------------
    fit_tmp2 <- tryCatch({
      message(sprintf("[INFO] Adding LOO/WAIC criteria for model %s ...", model_name))
      add_criterion(fit_tmp, c("loo", "waic"))
    }, error = function(e2) {
      message(sprintf(
        "[WARNING] add_criterion failed for model %s: %s",
        model_name, e2$message
      ))
      fit_tmp
    })
    
    # Save updated model (with criteria if successful)
    fits[[model_name]] <- fit_tmp2
    
    saveRDS(fit_tmp2, final_path)
    if (file.exists(final_path)) {
      message(sprintf("[OK] Final model saved to %s", final_path))
    } else {
      warning(sprintf("[WARNING] Final model for %s not found after saveRDS.", model_name))
    }
    
  }, error = function(e) {
    message(sprintf("[ERROR] Model %s failed: %s", model_name, e$message))
  })
}

# Save the combined list of all fitted models
if (length(fits) > 0L) {
  all_path <- file.path(output_dir, "basin_darwinian_all.rds")
  saveRDS(fits, all_path)
  if (file.exists(all_path)) {
    message(sprintf("[OK] All models list saved to %s", all_path))
  } else {
    warning("[WARNING] basin_linnaean_all.rds not found after saveRDS.")
  }
} else {
  warning("[WARNING] No models were successfully fitted; basin_darwinian_all.rds not written.")
}




# Model Comparison & Reporting -------------------------------------------------
# Generate unified comparison report containing both metrics
comparison_loo <- loo_compare(fits$lognormal, 
                              fits$gamma,
                              fits$exponential,
                              fits$weibull,
                              model_names = names(fits), 
                              criterion = "loo")
comparison_waic <- loo_compare(fits$lognormal, 
                               fits$gamma,
                               fits$exponential,
                               fits$weibull,
                               model_names = names(fits), 
                               criterion = "waic")

# Create combined output document
report_content <- capture.output({
  cat("=== Bayesian Model Evaluation Report ===\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M"), "\n\n")
  
  # LOO Results Section
  cat("--- Leave-One-Out Cross Validation (LOOIC) ---\n")
  print(comparison_loo, simplify = FALSE)
  cat("\nLOO Interpretation:\n")
  cat("* Models with LOOIC difference < 2SE are comparable\n")
  cat("* Negative elpd_diff indicates worse performance\n\n")
  
  # WAIC Results Section
  cat("\n--- Widely Applicable Information Criterion (WAIC) ---\n")
  print(comparison_waic, simplify = FALSE)
  cat("\nWAIC Interpretation:\n")
  cat("* Useful complement to LOO but more parametric\n")
  cat("* Consistent results increase confidence\n")
  
  # Consensus Analysis
  cat("\n--- Final Recommendation ---\n")
  best_loo <- rownames(comparison_loo)[1]
  best_waic <- rownames(comparison_waic)[1]
  
  if(best_loo == best_waic) {
    cat("Consensus Model:", best_loo, "\n")
  } else {
    cat("Divergent Results:\n")
    cat("- LOO Recommends:", best_loo, "\n")
    cat("- WAIC Recommends:", best_waic, "\n")
    cat("Action: Prioritize LOO with WAIC consistency check\n")
  }
})

# Save unified output
output_dir <- "output/logs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_path <- file.path(output_dir, "basin_darwinian_model_comparison")
writeLines(report_content, paste0(output_path, ".txt"))
writeLines(c("# Model Evaluation Report\n", gsub("=== ", "## ", report_content)), 
           paste0(output_path, ".md"))

# Terminal confirmation
message(sprintf(
  "\nAnalysis complete. Unified report saved to:\n- TXT: %s\n- MD: %s",
  paste0(output_path, ".txt"),
  paste0(output_path, ".md")
))
