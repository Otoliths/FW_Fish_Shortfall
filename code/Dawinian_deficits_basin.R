# PD Deficit Calculation Script
# This script computes Dawinian (Phylogenetic Diversity, PD) Deficits for freshwater fish assemblages
# across flexible spatial grouping levels, following Nakamura et al. (2022). 
# **FishPhyloMaker:::PD_defict**


# 1. Core PD Deficit Formula
# The PD deficit quantifies the relative contribution of imputed vs. non-imputed species to phylogenetic diversity.
# It ranges from 0 to 1:
#   - Values near 0: Most diversity from non-imputed species (high completeness)
#   - Values near 1: Most diversity from imputed species (low completeness)
# Formula:
#   PD_deficit = PD_imputed / (PD_imputed + PD_non_imputed)
# Where:
#   - PD_imputed: Phylogenetic diversity of species grafted at the genus or family level (graft = 1)
#   - PD_non_imputed: Phylogenetic diversity of species natively present in the megatree (graft = 0)


# 2. Supported Grouping Types
# The script adapts to 3 analysis modes via config$grouping_col:
#   - "biogeographic_realm": Calculate for biogeographic_realm (e.g., Afrotropic)
#   - "basin_id": Calculate for drainage basins (via unique basin_id codes)
#   - "none": Calculate a single global PD deficit for all species combined



# ==============================================================================
# 1. Package Setup (Same as Before)
# ==============================================================================
required_packages <- c("dplyr", "tidyr", "stringr", "openxlsx", "geiger", "ape", 
                       "purrr", "future", "furrr", "caper")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("Installing missing packages:", paste(missing_packages, collapse = ", "), "\n")
  install.packages(missing_packages, dependencies = TRUE)
}
lapply(required_packages, library, character.only = TRUE)

# Enable parallel computing (adjust workers based on your CPU)
plan(multisession, workers = 2)


# ==============================================================================
# 2. Centralized Configuration (with "no grouping" support)
# ==============================================================================
config <- list(
  # File paths
  data_path = "input/raw/cas_freshwater_v1.xlsx",
  biogeo_path = "input/raw/biogeographic_list.csv",
  tree_dir = "input/processed/tree_fish/",
  tree_prefix = "tree_fish_18k_n100_",
  tree_suffix = ".rds",
  output_csv = "output/tables/pd_deficits_all_basin.csv",  # Default for no grouping
  
  # Analysis parameters
  num_trees = 100,
  pd_root_edge = TRUE,  # Include root edge in PD calculations
  expected_total_species = 18821,  # Validation check: known total from raw data
  
  # Flexible grouping: 
  # - Use "biogeographic_realm" or "basin_id" for grouped analysis
  # - Use "none" to analyze all species together (no grouping)
  grouping_col = "basin_id"  # Current setting: no grouping (all species in one group)
)

# Auto-update output filename based on grouping type
if (config$grouping_col == "none") {
  config$output_csv <- "output/tables/pd_deficits_all_basin.csv"  # Explicit name for global analysis
} else {
  config$output_csv <- sprintf("output/tables/pd_deficits_%s.csv", config$grouping_col)  # Group-specific name
}


# ==============================================================================
# 3. Preload Trees (Same as Before—No Changes)
# ==============================================================================
cat("=== Preloading Trees ===\n")
tree_paths <- sprintf(
  "%s%s%s%s", 
  config$tree_dir, 
  config$tree_prefix, 
  1:config$num_trees,
  config$tree_suffix
)

# Load trees in parallel
trees <- future_map(tree_paths, function(path) {
  if (file.exists(path)) readRDS(path) else NULL
}, .progress = TRUE)
#names(trees) <- paste0("Tree_", 1:config$num_trees)
trees <- discard(trees, is.null)
cat(sprintf("Preloaded %d valid trees\n", length(trees)))


# ==============================================================================
# 4. Preprocess Main Data (Same as Before—No Changes)
# ==============================================================================
cat("\n=== Preprocessing Main Data ===\n")
raw_data <- openxlsx::read.xlsx(config$data_path)
biogeo_lookup <- read.csv(config$biogeo_path)

# Clean biogeographic data and merge with basin/biogeographic_realm info
df <- raw_data %>%
  mutate(basin = str_split(basin, ";")) %>%
  mutate(basin = lapply(basin, unique)) %>%
  mutate(basin = lapply(basin, sort)) %>%
  unnest(basin) %>%
  dplyr::select(valid_name, basin) %>%
  left_join(biogeo_lookup, by = "basin") %>%
  mutate(valid_name = gsub(" ", "_", valid_name)) %>%
  {if (config$grouping_col != "none")  # Check if we're using a real grouping column
    drop_na(., all_of(config$grouping_col))  # Drop NA for the grouping column (e.g., "basin_id")
    else .}  # No grouping: Keep all rows (no column to check for NA)

# Add graft status (from reference tree)
reference_tree <- trees[[1]]
graft_status <- reference_tree$graft_status %>%
  mutate(graft = ifelse(str_detect(tip_label, "[*]+$"), 1, 0)) %>%
  drop_na(species, graft) %>%
  dplyr::select(species, graft)

# Merge graft status and remove duplicates
df <- df %>%
  dplyr::select(
    valid_name, 
    if (config$grouping_col != "none") all_of(config$grouping_col)  # Conditional grouping column
  ) %>%
  left_join(graft_status, by = c("valid_name" = "species")) %>%
  drop_na(graft) %>%
  distinct()  # Avoid duplicate species-group pairs

cat(sprintf("Final data: %d species-%s entries\n", nrow(df), config$grouping_col))


# ==============================================================================
# 5. Preprocessing Function (supports grouped and ungrouped analysis)
# ==============================================================================
preprocess_group_data <- function(df, group_col) {
  # Case 1: No grouping ("none") – treat all species as a single global group
  if (group_col == "none") {
    # Create a dummy group label "All" to standardize structure with grouped data
    df_with_dummy_group <- df %>% mutate(Group = "All")
    
    return(df_with_dummy_group %>%
             group_by(Group) %>%  # Single group: "All"
             group_nest(.key = "species_data") %>%  # Nest all species data under the dummy group
             mutate(
               # Extract unique species + graft status (avoids duplicates in PD calculations)
               species_graft = map(species_data, ~distinct(.x, valid_name, graft)),
               # Calculate richness (total unique species in the global group)
               Richness = map_int(species_graft, nrow),
               # Precompute lists of imputed (graft=1) and non-imputed (graft=0) species
               tips_imp = map(species_graft, ~filter(.x, graft == 1)$valid_name),
               tips_non_imp = map(species_graft, ~filter(.x, graft == 0)$valid_name),
               # Remove raw nested data to save memory
               species_data = NULL
             ) %>%
             filter(Richness > 0))  # Ensure the global group has species (sanity check)
  }
  
  # Case 2: Grouped analysis (biogeographic_realm or basin_id)
  df %>%
    group_by(!!sym(group_col)) %>%  # Dynamically group by the specified column
    group_nest(.key = "species_data") %>%  # Nest species data within each group
    mutate(
      # Extract unique species + graft status per group
      species_graft = map(species_data, ~distinct(.x, valid_name, graft)),
      # Richness = number of unique species in the group
      Richness = map_int(species_graft, nrow),
      # Precompute imputed/non-imputed species lists for fast PD lookup
      tips_imp = map(species_graft, ~filter(.x, graft == 1)$valid_name),
      tips_non_imp = map(species_graft, ~filter(.x, graft == 0)$valid_name),
      # Clean up unused data
      species_data = NULL
    ) %>%
    filter(Richness > 0) %>%  # Skip groups with no species (avoid empty calculations)
    rename(Group = !!sym(group_col))  # Standardize group column name to "Group" for consistency
}
# ==============================================================================
# Example usage (integrates with downstream analysis)
# ==============================================================================
# Preprocess data based on config (works for grouped or ungrouped)
cat(sprintf("\n=== Preprocessing Data (Grouping: %s) ===\n", 
            ifelse(config$grouping_col == "none", "Global (no grouping)", config$grouping_col)))

group_data <- preprocess_group_data(df, config$grouping_col)

# Verify preprocessing output
cat(sprintf("Preprocessing complete: %d group(s) detected\n", nrow(group_data)))
if (config$grouping_col == "none") {
  cat(sprintf("Global group richness: %d species\n", group_data$Richness))
} else {
  cat(sprintf("Range of group richness: %d – %d species\n", 
              min(group_data$Richness), max(group_data$Richness)))
}


# ==============================================================================
# 6. FLEXIBLE PD Calculation Helper (Works for Any Group)
# ==============================================================================
calc_pd_for_group_tree <- function(group_row, tree, tree_idx, group_col_name) {
  # --------------------------
  # 1. Extract inputs and initialize debug info
  # --------------------------
  group_name <- group_row$Group  # Group identifier (e.g., "All", specific basin ID)
  tips_imp <- unlist(group_row$tips_imp)  # List of grafted/imputed species (ensure vector format)
  tips_non_imp <- unlist(group_row$tips_non_imp)  # List of natively present/non-imputed species
  richness <- group_row$Richness
  
  # Debug header (improved readability)
  cat(sprintf("\n\n=== Processing group: %s (Tree ID: %d) ===", group_name, tree_idx))
  cat(sprintf("\n- Input imputed species count: %d", length(tips_imp)))
  cat(sprintf("\n- Input non-imputed species count: %d", length(tips_non_imp)))
  cat(sprintf("\n- Overlapping species in lists (anomaly): %d", length(intersect(tips_imp, tips_non_imp))))
  
  # --------------------------
  # 2. Check species-tree matches
  # --------------------------
  tree_tips <- tree$tip.label
  shared_imp <- intersect(tips_imp, tree_tips)  # Imputed species present in the tree
  shared_non_imp <- intersect(tips_non_imp, tree_tips)  # Non-imputed species present in the tree
  total_shared <- length(shared_imp) + length(shared_non_imp)
  
  cat(sprintf("\n- Matched imputed species in tree: %d", length(shared_imp)))
  cat(sprintf("\n- Matched non-imputed species in tree: %d", length(shared_non_imp)))
  
  # --------------------------
  # 3. Handle extreme case: No matched species
  # --------------------------
  if (total_shared == 0) {
    cat("\n- Extreme case: No matched species → PD_deficit = NA")
    return(tibble(
      !!sym(ifelse(group_col_name == "none", "Group", group_col_name)) := group_name,
      Tree = tree_idx,
      Richness = richness,
      PD_deficit = NA_real_  # Explicit numeric NA
    ))
  }
  
  # --------------------------
  # 4. Prune tree and calculate PD
  # --------------------------
  # Retain only matched species in the tree (prevent empty tree)
  cleaned_tree <- keep.tip(tree, union(shared_imp, shared_non_imp))
  
  # Check if pruned tree is valid (extreme case: tree has no branches after pruning)
  if (length(cleaned_tree$tip.label) == 0 || sum(cleaned_tree$edge.length) == 0) {
    cat("\n- Extreme case: Pruned tree has no valid branches → PD_deficit = NA")
    return(tibble(
      !!sym(ifelse(group_col_name == "none", "Group", group_col_name)) := group_name,
      Tree = tree_idx,
      Richness = richness,
      PD_deficit = NA_real_
    ))
  }
  
  # Calculate clade matrix and PD values
  clmat <- caper::clade.matrix(cleaned_tree)
  pd_imp <- pd.calc(clmat, tip.subset = shared_imp, root.edge = config$pd_root_edge)
  pd_non_imp <- pd.calc(clmat, tip.subset = shared_non_imp, root.edge = config$pd_root_edge)
  
  # Check for abnormal PD results (e.g., NA or negative values)
  if (is.na(pd_imp) || is.na(pd_non_imp) || pd_imp < 0 || pd_non_imp < 0) {
    cat(sprintf("\n- Extreme case: Abnormal PD values (imputed=%.2f, non-imputed=%.2f) → PD_deficit = NA", pd_imp, pd_non_imp))
    return(tibble(
      !!sym(ifelse(group_col_name == "none", "Group", group_col_name)) := group_name,
      Tree = tree_idx,
      Richness = richness,
      PD_deficit = NA_real_
    ))
  }
  
  cat(sprintf("\n- Calculated PD_imputed: %.2f", pd_imp))
  cat(sprintf("\n- Calculated PD_non_imputed: %.2f", pd_non_imp))
  
  # --------------------------
  # 5. Calculate PD Deficit (refined extreme cases)
  # --------------------------
  pd_total <- pd_imp + pd_non_imp
  
  # Only non-imputed species present → deficit = 0 (no contribution from grafted species)
  if (length(shared_imp) == 0) {
    pd_deficit <- 0
    cat("\n- Extreme case: Only non-imputed species → PD_deficit = 0")
  }
  # Only imputed species present → deficit = 1 (all PD from grafted species)
  else if (length(shared_non_imp) == 0) {
    pd_deficit <- 1
    cat("\n- Extreme case: Only imputed species → PD_deficit = 1")
  }
  # Total PD = 0 (theoretically rare; avoid division by zero)
  else if (pd_total == 0) {
    pd_deficit <- NA_real_
    cat("\n- Extreme case: Total PD = 0 → PD_deficit = NA")
  }
  # Normal calculation
  else {
    pd_deficit <- pd_imp / pd_total
    cat(sprintf("\n- Normal calculation: PD_deficit = %.4f", pd_deficit))
  }
  
  # --------------------------
  # 6. Return results (unified format)
  # --------------------------
  tibble(
    !!sym(ifelse(group_col_name == "none", "Group", group_col_name)) := group_name,
    Tree = tree_idx,
    Richness = richness,
    PD_deficit = pd_deficit
  )
}

# ==============================================================================
# 7. Parallel Calculation (Flexible for Any Group)
# ==============================================================================
cat(sprintf("\n=== Running Parallel PD Calculation for %s ===\n", config$grouping_col))
# Create all (group, tree) pairs to process
pair_data <- expand_grid(
  group_row = transpose(group_data),  # Convert group data to list rows
  tree_info = transpose(tibble(
    tree = trees,
    tree_idx = as.integer(str_remove(names(trees), "Tree_"))
  ))
)

# Process pairs in parallel
final_results <- future_map2_dfr(
  .x = pair_data$group_row,
  .y = pair_data$tree_info,
  .f = function(group_row, tree_info) {
    calc_pd_for_group_tree(
      group_row = group_row,
      tree = tree_info$tree,
      tree_idx = tree_info$tree_idx,
      group_col_name = config$grouping_col  # Pass grouping column name to helper
    )
  },
  .progress = TRUE  # Show progress bar
)


# ==============================================================================
# 8. Final Output (Flexible to Grouping Column)
# ==============================================================================
cat(sprintf("\n=== Compiling Results for %s ===\n", config$grouping_col))
# Save results to CSV (filename matches grouping)
write.csv(final_results, config$output_csv, row.names = FALSE)
cat(sprintf("Results saved to: %s\n", config$output_csv))

# Print summary stats
# Generate summary statistics (works for grouped AND ungrouped analysis)
results_summary <- final_results %>%
  # Dynamically select grouping column:
  # - Use the specified column (e.g., "basin_id", "biogeographic_realm") for grouped analysis
  # - Use the dummy "Group" column (value = "All") for ungrouped/global analysis
  group_by(
    !!sym(ifelse(
      config$grouping_col == "none",  # Check if we're in "no grouping" mode
      "Group",                        # Ungrouped: Use dummy "Group" column
      config$grouping_col             # Grouped: Use specified column (e.g., "basin_id")
    ))
  ) %>%
  summarise(
    N_Trees = n(),  # Number of trees analyzed per group (should = config$num_trees for valid groups)
    Mean_PD_Deficit = round(mean(PD_deficit, na.rm = TRUE), 4),  # Avg PD deficit (ignore NAs from missing species)
    Richness = first(Richness),  # Richness is constant per group (safe to use `first()`)
    .groups = "drop"  # Clean up grouping metadata to avoid unintended behavior
  ) 
# results_summary <- results_summary %>% left_join(biogeo_lookup, by = "basin_id")
# Show top 10 groups (for quick inspection)
cat(sprintf("\nTop 10 %s by Mean PD Deficit:\n", config$grouping_col))
print(head(results_summary, 10))


# ==============================================================================
# 9. Clean Up
# ==============================================================================
plan(sequential)  # Disable parallelism
rm(trees, group_data, pair_data, reference_tree)
gc()
cat(sprintf("\n=== Analysis Complete for %s ===\n", config$grouping_col))