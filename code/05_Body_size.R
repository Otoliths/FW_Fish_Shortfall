library(rfishbase)   # version 5.0.2
library(dplyr)
library(Rphylopars)

## -----------------------------------------
## 1. Download / prepare FishBase trait data
## -----------------------------------------

# 1.1 Get FishBase species table (has SpecCode, Length, etc.)
sp <- species()

# 1.2 Get FishBase taxonomic names (SpecCode <-> Species)
sp_names <- load_taxa()

# 1.3 Save raw FishBase objects if needed for reproducibility
# save.image("input/raw/Fishbase-species.RData")  # not ideal but you used it originally
# load("input/raw/Fishbase-species.RData")

# 1.4 Merge Length info and species names
trait <- sp %>%
  dplyr::select(SpecCode, Length) %>%        # keep only SpecCode + max length
  left_join(sp_names, by = "SpecCode")       # brings in Species name

## trait now has: SpecCode, Length, Species, plus possibly extra taxonomic columns


## -----------------------------------------
## 2. Load CAS freshwater species list
## -----------------------------------------

# CAS freshwater list of valid species names
cas <- openxlsx::read.xlsx("input/raw/cas_freshwater_v1.xlsx")

# Join CAS names to FishBase max body length (in cm)
df <- cas %>%
  dplyr::select(valid_name) %>%
  left_join(trait, by = c("valid_name" = "Species"))

# Optional: coverage summary of Length completeness
df %>%
  summarise(
    total      = n(),
    na_count   = sum(is.na(Length)),
    na_percent = mean(is.na(Length)) * 100
  )
# Example you observed:
# total na_count na_percent
# 18821     3486   18.52186



## ------------------------------------------------------------
## Phylogeny-informed multiple imputation of fish body size
##
## Purpose:
##   - Impute missing fish maximum length values using Rphylopars
##   - Propagate phylogenetic uncertainty across 100 trees
##   - Generate multiple imputed datasets (MI), rather than a
##     single tree-averaged point estimate
##
## Assumptions:
##   - df has columns:
##       * valid_name  : species name in "Genus species" format
##       * Length      : maximum length in cm
##   - Tree files:
##       input/processed/tree_fish/tree_fish_18k_n100_1.rds
##       ...
##       input/processed/tree_fish/tree_fish_18k_n100_100.rds
##
## Outputs:
##   - mi_list         : list of multiply imputed datasets
##   - final_single    : single tree-averaged dataset (for checking only)
##   - saved RDS files : one RDS per imputed dataset
## ------------------------------------------------------------

## ------------------------------------------------------------
## 0. Load required packages
## ------------------------------------------------------------
library(dplyr)
library(Rphylopars)

## ------------------------------------------------------------
## 1. User-defined settings
## ------------------------------------------------------------

# Number of phylogenetic trees
n_tree <- 100

# Number of multiple imputations to generate
m <- 20

# Random seed for reproducibility
set.seed(123)

# Directory to save MI datasets
outdir <- "input/data_prep/mi_body_size"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## ------------------------------------------------------------
## 2. Build vector of tree files
## ------------------------------------------------------------
tree_files <- sprintf(
  "input/processed/tree_fish/tree_fish_18k_n100_%d.rds",
  1:n_tree
)

# Check that all tree files exist
missing_tree_files <- tree_files[!file.exists(tree_files)]
if (length(missing_tree_files) > 0) {
  stop(
    "The following tree files are missing:\n",
    paste(missing_tree_files, collapse = "\n")
  )
}

## ------------------------------------------------------------
## 3. Prepare species names for tree matching
## ------------------------------------------------------------

# Convert species names from "Genus species" to "Genus_species"
df <- df %>%
  mutate(
    tree_name = gsub(" ", "_", valid_name)
  )

## ------------------------------------------------------------
## 4. Use the first tree to identify species present in phylogeny
## ------------------------------------------------------------
tre1 <- readRDS(tree_files[1])
tre1$tip.label <- gsub("[*]*$", "", tre1$tip.label)

in_tree_idx <- which(df$tree_name %in% tre1$tip.label)

# Species in the phylogeny
trait_data <- df[in_tree_idx, c("tree_name", "Length")]
names(trait_data)[1] <- "species"  # Required by Rphylopars

# Species not in the phylogeny cannot be phylogenetically imputed
remain_data <- df[-in_tree_idx, ]

# Sanity check
stopifnot("Length" %in% colnames(trait_data))

## ------------------------------------------------------------
## 5. Work on log scale for body size imputation
##
## Reason:
##   - Maximum length is typically strongly right-skewed
##   - Log transformation better matches Gaussian assumptions
## ------------------------------------------------------------
trait_data <- trait_data %>%
  mutate(
    log_Length = ifelse(is.na(Length), NA_real_, log10(Length + 1))
  ) %>%
  select(species, log_Length)

# Identify observed and missing species within the tree
obs_idx <- which(!is.na(trait_data$log_Length))
mis_idx <- which(is.na(trait_data$log_Length))

## ------------------------------------------------------------
## 6. Fit Rphylopars across all trees
##
## For each tree, store:
##   - Predicted tip values for all species
##   - Residual SD estimated from observed species
##
## Note:
##   This residual SD is used later to generate stochastic draws
##   for multiple imputation.
## ------------------------------------------------------------
n_spp <- nrow(trait_data)

pred_mat  <- matrix(NA_real_, nrow = n_spp, ncol = n_tree)
sigma_vec <- rep(NA_real_, n_tree)

rownames(pred_mat) <- trait_data$species
colnames(pred_mat) <- paste0("tree_", seq_len(n_tree))

for (i in seq_along(tree_files)) {
  cat("Processing tree", i, "of", n_tree, "\n")
  
  tre <- readRDS(tree_files[i])
  tre$tip.label <- gsub("[*]*$", "", tre$tip.label)
  
  # Fit phylogenetic trait model on the current tree
  PPE <- suppressWarnings(
    phylopars(
      trait_data = trait_data[, c("species", "log_Length")],
      tree       = tre
    )
  )
  
  # Extract reconstructed values for tips and ancestors
  sim_data <- as.data.frame(PPE[["anc_recon"]])
  sim_data$species <- rownames(sim_data)
  
  # Keep only tip predictions for the species in trait_data
  sim_tips <- sim_data[
    sim_data$species %in% trait_data$species,
    c("species", "log_Length")
  ]
  
  # Reorder predictions to match the original species order
  sim_tips <- sim_tips[match(trait_data$species, sim_tips$species), ]
  
  # Store predicted values
  pred_mat[, i] <- sim_tips$log_Length
  
  # Estimate tree-specific residual SD using observed species only
  obs_pred <- sim_tips$log_Length[obs_idx]
  obs_true <- trait_data$log_Length[obs_idx]
  
  sigma_vec[i] <- sd(obs_true - obs_pred, na.rm = TRUE)
}

# Replace zero or missing sigma values with a global fallback
global_sigma <- median(sigma_vec[is.finite(sigma_vec) & sigma_vec > 0], na.rm = TRUE)

sigma_vec[!is.finite(sigma_vec) | sigma_vec <= 0] <- global_sigma

## ------------------------------------------------------------
## 7. Create a single tree-averaged dataset (for inspection only)
##
## Important:
##   This is NOT the formal MI result.
##   It is only useful for quick checks / descriptive summaries.
## ------------------------------------------------------------
mean_pred <- rowMeans(pred_mat, na.rm = TRUE)

impute_df_single <- data.frame(
  species            = rownames(pred_mat),
  log_Length_imputed = mean_pred,
  stringsAsFactors   = FALSE
) %>%
  mutate(
    valid_name = gsub("_", " ", species)
  )

final_single <- df %>%
  left_join(
    impute_df_single[, c("valid_name", "log_Length_imputed")],
    by = "valid_name"
  ) %>%
  mutate(
    log_length_max = ifelse(
      is.na(Length) & !is.na(log_Length_imputed),
      log_Length_imputed,
      log10(Length + 1)
    ),
    imputed_flag = ifelse(
      is.na(Length) & !is.na(log_Length_imputed),
      TRUE,
      FALSE
    )
  ) %>%
  select(
    valid_name,
    log_length_max,
    imputed_flag
  )
anyNA(final_single)
saveRDS(final_single,"input/data_prep/body_size_final_single.rds")

## ------------------------------------------------------------
## 8. Generate multiple imputed datasets
##
## Strategy:
##   - For each missing species, randomly select one tree
##   - Use that tree-specific predicted mean
##   - Add stochastic noise using the tree-specific residual SD
##
## This yields m complete datasets for downstream MI analysis.
## ------------------------------------------------------------
mi_list <- vector("list", length = m)

# Species-specific SD across trees
row_sd_pred <- apply(pred_mat, 1, sd, na.rm = TRUE)

# Global fallback SD
global_sd_pred <- median(row_sd_pred[is.finite(row_sd_pred) & row_sd_pred > 0], na.rm = TRUE)
if (!is.finite(global_sd_pred)) {
  stop("global_sd_pred is NA or non-finite; check pred_mat.")
}

# Replace bad species-specific SDs
row_sd_pred[!is.finite(row_sd_pred) | row_sd_pred <= 0] <- global_sd_pred

for (j in seq_len(m)) {
  cat("Generating imputed dataset", j, "of", m, "\n")
  
  # Start from the original log-scale trait vector
  imputed_values <- trait_data$log_Length
  
  if (length(mis_idx) > 0) {
    # Randomly choose one tree for each missing species
    sampled_tree <- sample(seq_len(n_tree), size = length(mis_idx), replace = TRUE)
    
    # Tree-specific predicted means for missing species
    mu_hat <- pred_mat[cbind(mis_idx, sampled_tree)]
    
    # Fallback 1: replace NA tree-specific means with row means across trees
    row_mean_pred <- rowMeans(pred_mat, na.rm = TRUE)
    bad_mu <- is.na(mu_hat) | !is.finite(mu_hat)
    if (any(bad_mu)) {
      mu_hat[bad_mu] <- row_mean_pred[mis_idx[bad_mu]]
    }
    
    # Use species-specific SD across trees as uncertainty scale
    sigma <- row_sd_pred[mis_idx]
    
    # Replace bad sigma values with the global fallback SD
    bad_sigma <- is.na(sigma) | !is.finite(sigma) | sigma <= 0
    if (any(bad_sigma)) {
      sigma[bad_sigma] <- global_sd_pred
    }
    
    # Generate draws only for valid entries
    ok <- !is.na(mu_hat) & is.finite(mu_hat) &
      !is.na(sigma) & is.finite(sigma) & sigma > 0
    
    draws <- rep(NA_real_, length(mis_idx))
    if (any(ok)) {
      draws[ok] <- rnorm(sum(ok), mean = mu_hat[ok], sd = sigma[ok])
      draws[ok] <- pmax(draws[ok], 0)
    }
    
    # Fill missing entries
    imputed_values[mis_idx] <- draws
  }
  
  tmp <- data.frame(
    species = trait_data$species,
    log_Length_orig = trait_data$log_Length,
    log_Length_mi = imputed_values,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      valid_name = gsub("_", " ", species),
      imputed_flag = is.na(log_Length_orig) & !is.na(log_Length_mi)
    ) %>%
    select(valid_name, log_Length_mi, imputed_flag)
  
  final_mi <- df %>%
    select(valid_name, Length) %>%
    left_join(tmp, by = "valid_name") %>%
    mutate(
      log_length_max = ifelse(
        !is.na(Length),
        log10(Length + 1),
        log_Length_mi
      ),
      imputed_flag = ifelse(is.na(imputed_flag), FALSE, imputed_flag)
    ) %>%
    select(valid_name, log_length_max, imputed_flag)
  
  mi_list[[j]] <- final_mi
  
  saveRDS(
    final_mi,
    file = file.path(outdir, sprintf("body_size_mi_%02d.rds", j))
  )
}



## ------------------------------------------------------------
## 9. Optional summary checks
## ------------------------------------------------------------

# Quick preview of the first MI dataset
head(mi_list[[1]])

# Compare number of imputed species in each dataset
sapply(mi_list, function(x) sum(x$imputed_flag, na.rm = TRUE))

# Summary of tree-specific residual SD
summary(sigma_vec)

# Number of species that could not be phylogenetically imputed
# because they were absent from the phylogeny
n_not_in_tree <- nrow(remain_data)
cat("Species not present in the phylogeny:", n_not_in_tree, "\n")

## ------------------------------------------------------------
## 10. Objects returned in session
##
## - pred_mat      : tree-specific predicted log body size
## - sigma_vec     : tree-specific residual SD
## - final_single  : single tree-averaged dataset (inspection only)
## - mi_list       : list of m multiply imputed datasets
## ------------------------------------------------------------

# FishBase body size
# │
# missing values
# │
# Rphylopars across 100 trees
# │
# generate 20 MI datasets
# │
# run model on each dataset
# │
# Rubin’s rules
# │
# final parameter estimates



m <- nrow(fit_results)

# pooled estimate
Q_bar <- mean(fit_results$beta)

# within-imputation variance
U_bar <- mean(fit_results$se^2)

# between-imputation variance
B <- var(fit_results$beta)

# total variance
T_var <- U_bar + (1 + 1/m) * B

# pooled standard error
T_se <- sqrt(T_var)

# 95% CI
lower <- Q_bar - 1.96 * T_se
upper <- Q_bar + 1.96 * T_se

data.frame(
  estimate = Q_bar,
  se = T_se,
  CI_low = lower,
  CI_high = upper
)

# To account for uncertainty in missing body-size traits, we generated 20 multiple-imputed datasets using 
# phylogenetic predictions and stochastic draws from a normal distribution. 
# All stochastic procedures were conducted with a fixed random seed to ensure full reproducibility of the imputation process.

# Missing body-size values were handled using phylogeny-informed multiple imputation. 
# We generated multiple imputed datasets and fitted the same Bayesian model to each dataset using brms::brm_multiple. 
# Following the Bayesian multiple-imputation workflow implemented in brms, pooled inference was obtained by combining 
# posterior draws across submodels, thereby propagating imputation uncertainty into final parameter estimates.

## ------------------------------------------------------------
## Phylogenetic imputation across 100 trees + imputation flag
##
## Assumes:
##   - df has columns: valid_name (e.g. "Genus species") and Length (mm)
##   - Tree files: tree_fish_18k_n100_1.rds ... tree_fish_18k_n100_100.rds
## ------------------------------------------------------------

## 3. Prepare data for phylogenetic imputation -----------------

# 3.1 Build vector of tree files (1–100)
tree_files <- sprintf(
  "input/processed/tree_fish/tree_fish_18k_n100_%d.rds",
  1:100
)

# 3.2 Make sure species naming style matches the tree
# Tree uses "Genus_species"; df$valid_name is "Genus species"
df <- df %>%
  mutate(
    tree_name = gsub(" ", "_", valid_name)   # for matching tree tip labels
  )

# 3.3 Use the first tree to identify which species are in the phylogeny
tre1 <- readRDS(tree_files[1])
tre1$tip.label <- gsub("[*]*$", "", tre1$tip.label)

in_tree_idx <- which(df$tree_name %in% tre1$tip.label)

trait_data <- df[in_tree_idx, c("tree_name", "Length")]
names(trait_data)[1] <- "species"  # Rphylopars expects first col = species

# Species not present in the tree cannot be phylogenetically imputed
remain_data <- df[-in_tree_idx, ]

# Quick check: make sure the trait column name is "Length"
# If your trait column has a different name, update it here accordingly.
stopifnot("Length" %in% colnames(trait_data))

## 4. Run phylogenetic model across 100 trees ------------------

# Matrix to store predictions:
# rows = species in tree, cols = 100 phylogenetic replicates
n_spp  <- nrow(trait_data)
n_tree <- length(tree_files)

pred_mat <- matrix(NA_real_, nrow = n_spp, ncol = n_tree)
rownames(pred_mat) <- trait_data$species

for (i in seq_along(tree_files)) {
  tre <- readRDS(tree_files[i])
  tre$tip.label <- gsub("[*]*$", "", tre$tip.label)
  
  # Fit Rphylopars on this tree
  PPE <- suppressWarnings(
    phylopars(
      trait_data = trait_data[, c("species", "Length")],
      tree       = tre
    )
  )
  
  # Extract reconstructed trait values (tips + ancestors)
  sim_data <- as.data.frame(PPE[["anc_recon"]])
  sim_data$species <- rownames(sim_data)
  
  # Keep only tips corresponding to our species set
  sim_tips <- sim_data[sim_data$species %in% trait_data$species,
                       c("species", "Length")]
  
  # Align order to trait_data$species
  sim_tips <- sim_tips[match(trait_data$species, sim_tips$species), ]
  
  # Store predictions for this tree
  pred_mat[, i] <- sim_tips$Length
}

# Mean prediction across 100 trees (species-level mean)
mean_pred <- rowMeans(pred_mat, na.rm = TRUE)

impute_df <- data.frame(
  species        = rownames(pred_mat),
  Length_imputed = mean_pred,
  stringsAsFactors = FALSE
)

# Convert "Genus_species" back to "Genus species" for joining
impute_df <- impute_df %>%
  mutate(valid_name = gsub("_", " ", species))

## 5. Combine observed + imputed Length + imputed_flag ---------

final <- df %>%
  # Join imputed values back to full species list
  left_join(impute_df[, c("valid_name", "Length_imputed")],
            by = "valid_name") %>%
  mutate(
    # If original Length is missing but we have an imputed value → use imputed
    Length_final = ifelse(
      is.na(Length) & !is.na(Length_imputed),
      Length_imputed,
      Length
    ),
    # Flag imputed records: TRUE = originally missing and filled by phylogenetic imputation
    imputed_flag = ifelse(
      is.na(Length) & !is.na(Length_imputed),
      TRUE,
      FALSE
    ),
    # Compute log-transformed body size for downstream analyses
    log_length_max = log(Length_final + 1)
  ) %>%
  # Keep only the key columns you need; adjust as needed
  select(
    valid_name,
    log_length_max,
    imputed_flag
  )

# Preview
head(final)

# Save output table for downstream use (including imputation flag)
write.csv(final, "input/data_prep/body_size_flag.csv",row.names = F)
saveRDS(final[,1:2], "input/data_prep/body_size.rds")
