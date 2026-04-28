# ------------------------------------------------------------
# Load required libraries for data processing and modeling
# ------------------------------------------------------------
library(dplyr)        # Data manipulation verbs
library(magrittr)     # Pipe operator support (%>% and .)
library(data.table)   # Fast data handling for large tables
library(brms)         # Interface for Bayesian multilevel models using Stan
library(scales)       # Scaling helper functions (e.g., rescale)
library(arrow)        # Read/write parquet and other columnar data formats

# Use cmdstanr backend for brms and enable debugging output
options(brms.backend = "cmdstanr", brms.debug = TRUE)

# ------------------------------------------------------------
# 1. Load input data
# ------------------------------------------------------------

# Import raw species data from Excel file (CAS freshwater list)
cas <- openxlsx::read.xlsx("input/raw/cas_freshwater_v1.xlsx")

# Load pre-processed country-level data from parquet
base_df <- read_parquet("input/data_prep/linnaean_country_single.parquet")

# Merge CAS species attributes (e.g. valid_name, year of description) into country data
# Column indices 5 and 8 in 'cas' should correspond to valid_name and family
data <- base_df %>% 
  left_join(cas[, c(5, 8)], by = "valid_name")


# ------------------------------------------------------------
# 2. Basic filtering and construction of survival time
# ------------------------------------------------------------

# - Keep only complete cases for the selected columns
# - Compute discovery time since Linnaeus (year 1758)
# - Remove non-positive times (pre-Linnaean or problematic records)
# - Create event indicator (=1 for all records in this fully observed dataset)
data <- data %>% 
  filter(year_description > 1758) %>%
  mutate(time = year_description - 1758) %>%  # Discovery time since Linnaeus (in years)
  filter(time > 0) %>%                        # Remove records with non-positive discovery time
  mutate(event = 1) %>%                        # Event indicator for survival model (no censoring here)
  .[complete.cases(.), ]

data$country_area <- log10(data$country_area+1)
data$range_size <- log10(data$range_size+1)  
data$elevation <- log10(data$elevation+300)
data$discharge <- log10(data$discharge+1)
data$population_density <- log10(data$population_density+1)
data$preserved_specimen <- log10(data$preserved_specimen+1)
data$sequencing_effort <- log10(data$sequencing_effort+1)
data$sampling_effort <- log10(data$sampling_effort+1)
data$taxonomic_effort <- log10(data$taxonomic_effort+1)
data$taxonomic_activity <- log10(data$taxonomic_activity+1)
data$watertemp <- log10(data$watertemp+1)
data$range_rarity <- log10(data$range_rarity+1)
  
  # ------------------------------------------------------------
# 3. Define predictor set and transformation helpers
# ------------------------------------------------------------

# List of continuous predictors to be normalized or transformed
# 'log_length_max' is assumed to be already on log scale
vars <- c(
  "log_length_max", "taxonomic_effort", "taxonomic_activity",
  "country_area", "range_size", "elevation", "latitude",
  "discharge", "watertemp", "preserved_specimen",
  "range_rarity", "population_density",
  "sequencing_effort", "sampling_effort"
)
# data$latitude <- abs(data$latitude)
# Convert to plain data.frame for usdm::vif
# VIF < 5 is typically considered acceptable
usdm::vif(data[, vars] %>% as.data.frame())

# ------------------------------------------------------------
# 4. Apply transformations: create "z_*" columns
# ------------------------------------------------------------
df_transformed <- data %>%
  mutate(across(all_of(vars),
                ~ as.numeric(scale(.x)),
                .names = "z_{.col}"))


usdm::vif(df_transformed[,23:34] %>% as.data.frame())

# ------------------------------------------------------------
# 5. Assemble survival modeling dataset
# ------------------------------------------------------------

# Keep only the fields needed for survival modeling
# - country_id: drainage country identifier
# - continent: 6 major continent
# - names(preds): transformed predictors (all "z_*" columns)
# - family: taxonomic family
# - time: discovery time since 1758

stan_survdata <- df_transformed %>%
  dplyr::transmute(
    iso3                = factor(iso3),
    valid_name          = valid_name,
    continent           = factor(continent),
    
    body_size           = z_log_length_max,
    taxonmic_effort     = z_taxonmic_effort,
    taxonomic_activity  = z_taxonomic_activity,
    country_area        = z_country_area,
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
    time                = time
  ) %>%
  .[complete.cases(.), ] %>%
  as.data.frame() %>%
  distinct()


rm(base_df,cas,data,df_transformed);gc()
#saveRDS(stan_survdata,"stan_survdata_country_linnaean.rds")
# ------------------------------------------------------------
# 6. Specify Bayesian survival model formula
# ------------------------------------------------------------

# Linnaean discovery model with global + continent-level slopes
# - time: years since Linnaeus (1758) until species description
# - Global fixed effects + continent-level random slopes
# - Random intercepts: country_id, family_group

# Vector of 14 standardized predictors
preds <- c(
  "body_size",
  "taxonmic_effort",
  "taxonomic_activity",
  "country_area",
  "range_size",
  "elevation",
  "latitude",
  "discharge",
  "watertemp",
  "preserved_specimen",
  "sequencing_effort",
  "sampling_effort",
  "range_rarity",
  "population_density"
)

# Collapse predictors into a single string "x1 + x2 + ... + x14"
pred_str <- paste(preds, collapse = " + ")

# ------------------------------------------------------------
# Hierarchical survival model
#
#  • Global fixed effects for all 14 predictors
#  • continent-level random slopes for the same 14 predictors
#    (partial pooling avoids numerical instability in data-poor realms such as Oceania)
#  • country- and family-level random intercepts
# ------------------------------------------------------------

form <- bf(
  as.formula(
    paste(
      # survival likelihood with censoring indicator
      "time ~ 1 +",
      
      # fixed-effect interactions: all predictors × continent
      paste0("(", pred_str, ") * continent"),
      
      # hierarchical random intercepts
      "+ (1 | iso3)",
      "+ (1 | family_group)"
    )
  )
)


# ------------------------------------------------------------
# 7. MCMC configuration (shared across distributions)
# ------------------------------------------------------------

# These arguments are shared across all survival distributions and then
# combined with model-specific settings in 'model_configs'.


priors <- c(
  prior(normal(0, 1.5), class = "b"),          # Global slopes and all other fixed effects
  prior(normal(0, 1.5), class = "Intercept"),  # Global intercept
  prior(student_t(3, 0, 2), class = "sd")      # SDs of realm/basin/family random effects
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
  init      = 0.1,                           # Small positive initial values
  control = list(adapt_delta = 0.98,        # High target acceptance to reduce divergences
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
  list(model_name = "lognormal",   family = brmsfamily("lognormal"),   iter = 6000, warmup = 3000),
  list(model_name = "gamma",       family = brmsfamily("Gamma"),      iter = 6000, warmup = 3000),
  list(model_name = "exponential", family = brmsfamily("exponential"), iter = 6000, warmup = 3000),
  list(model_name = "weibull",     family = brmsfamily("weibull"),    iter = 6000, warmup = 3000)
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
  raw_path   <- file.path(output_dir, sprintf("country_linnaean_%s_raw.rds",  model_name))
  final_path <- file.path(output_dir, sprintf("country_linnaean_%s.rds",      model_name))
  
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
  all_path <- file.path(output_dir, "country_linnaean_all.rds")
  saveRDS(fits, all_path)
  if (file.exists(all_path)) {
    message(sprintf("[OK] All models list saved to %s", all_path))
  } else {
    warning("[WARNING] country_linnaean_all.rds not found after saveRDS.")
  }
} else {
  warning("[WARNING] No models were successfully fitted; country_linnaean_all.rds not written.")
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
output_path <- file.path(output_dir, "country_linnaean_model_comparison")
writeLines(report_content, paste0(output_path, ".txt"))
writeLines(c("# Model Evaluation Report\n", gsub("=== ", "## ", report_content)), 
           paste0(output_path, ".md"))

# Terminal confirmation
message(sprintf(
  "\nAnalysis complete. Unified report saved to:\n- TXT: %s\n- MD: %s",
  paste0(output_path, ".txt"),
  paste0(output_path, ".md")
))

################################################################################
source("code/functions/gompertz_family.R")  # updated 2025-10-08

form_gompertz <- bf(
  as.formula(
    paste(
      # survival likelihood with censoring indicator
      "time ~ 1 +",
      
      # fixed-effect interactions: all predictors × continent
      paste0("(", pred_str, ") * continent"),
      
      # hierarchical random intercepts
      "+ (1 | iso3)",
      "+ (1 | family_group)"
    )
  ),
  gamma ~ 1
)

fit_gomp <- brm(
  formula = form_gompertz,
  data    = stan_survdata,
  chains  = 8,
  cores   = 8,
  iter = 6000,
  warmup = 3000,
  seed    = 2025,
  backend = "cmdstanr",
  stan_model_args = list(
    cpp_options = list(
      "STAN_THREADS=TRUE",
      # "STAN_OPENCL=TRUE",           # Enable OpenCL (GPU) support if available
      # "OPENCL_DEVICE=0",            # Select OpenCL device
      # "OPENCL_PLATFORM=0",          # Select OpenCL platform
      "CXXFLAGS += -O3 -Wno-overloaded-virtual"  # Compiler optimization flags
    )
  ),
  threads   = threading(10, static = TRUE),  # Within-chain parallelization
  normalize = TRUE,                         #  keep full likelihood for LOO comparison
  thin      = 1,                            # Thinning to reduce autocorrelation / disk usage
  save_pars = save_pars(all = TRUE),         # Save all parameters (useful for post-processing)
  init      = 0.1,                          # Small positive initial values
  control = list(adapt_delta = 0.99,         # High target acceptance to reduce divergences
                 max_treedepth = 15          # ↑ Allow deeper NUTS trees; avoids hitting limit at 10
  ),
  prior     = priors_linnaean,
  family = gompertz_family,
  stanvars = stanvar(scode = stan_funs, block = "functions")
)

saveRDS(fit_gomp, file.path("output/model", "country_linnaean_gompertz.rds"))


# stan_survdata <- readRDS("output/stan_survdata_country_linnaean_multiple.rds")
# brm_multiple(
#   formula = form_gompertz,
#   data    = stan_survdata,
#   chains  = 8,
#   cores   = 8,
#   iter = 6000,
#   warmup = 3000,
#   seed    = 2025,
#   backend = "cmdstanr",
#   stan_model_args = list(
#     cpp_options = list(
#       "STAN_THREADS=TRUE",
#       # "STAN_OPENCL=TRUE",           # Enable OpenCL (GPU) support if available
#       # "OPENCL_DEVICE=0",            # Select OpenCL device
#       # "OPENCL_PLATFORM=0",          # Select OpenCL platform
#       "CXXFLAGS += -O3 -Wno-overloaded-virtual"  # Compiler optimization flags
#     )
#   ),
#   threads   = threading(10, static = TRUE),  # Within-chain parallelization
#   normalize = TRUE,                         #  keep full likelihood for LOO comparison
#   thin      = 1,                            # Thinning to reduce autocorrelation / disk usage
#   save_pars = save_pars(all = TRUE),         # Save all parameters (useful for post-processing)
#   init      = 0.1,                          # Small positive initial values
#   control = list(adapt_delta = 0.99,         # High target acceptance to reduce divergences
#                  max_treedepth = 15          # ↑ Allow deeper NUTS trees; avoids hitting limit at 10
#   ),
#   prior     = priors,
#   family = gompertz_family,
#   stanvars = stanvar(scode = stan_funs, block = "functions"),
#   combine = TRUE                           #fitted models should be combined into a single fitted model
# )

