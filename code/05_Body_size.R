library(rfishbase)   # version 5.0.1
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

# Join CAS names to FishBase Length
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
# 18821     3950    20.9872



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
