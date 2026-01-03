# ================================
#   Taxonomic Effort Index (TAE)
#   Part 1: Simple Example
#   Part 2: Full Dataset (18,821 species)
# ================================

# --- Load required packages ---
library(httr)      # 1.4.7
library(xml2)      # 1.3.6
library(dplyr)     # 1.1.4
library(purrr)     # 1.0.2
library(stringr)   # 1.5.1
library(cli)       # 3.6.3
library(tibble)    # 3.2.1
library(tidyr)     # 1.3.1
source("code/functions/TAE_function.R")  # updated 2025-10-08


# ========================================
# Part 1. Simple Example: Single Species
# ========================================

# --- 1. Fetch reference metadata ---
# Retrieve metadata for five example references of *Amatitlania kanna*.
ids <- c(29310, 29317, 33194, 33605, 34357)
refs_raw <- get_cas_ref(query = ids, quiet = FALSE)

# Inspect retrieved fields: authorship, year, reference text, and URL.
print(refs_raw)

# --- 2. Clean authorship strings ---
# Standardize name order and separators for consistent author counts.
refs_clean <- clean_refs(refs_raw, ref_authorship)

# --- 3. Compute TAE for this species ---
# TAE integrates temporal span, number of unique authors, and temporal decay.
tae <- calculate_TAE(refs_clean, lambda = log(2)/20, author_sep = ",")
tae  # expected ≈ 1.88434


# ========================================
# Part 2. Full Dataset: 18,821 Species
# ========================================

# --- 4. Compute TAE for all species ---

# Load reference metadata and species–reference mapping
cas_reference <- readRDS("input/raw/cas_reference.rds")
df <- openxlsx::read.xlsx("input/raw/cas_freshwater_v1.xlsx")

# Prepare reference IDs: split, deduplicate, sort, unnest, remove duplicates
df <- df %>%
  select(valid_name, authorship, references) %>%
  mutate(references = str_split(references, ", ")) %>%
  mutate(references = lapply(references, unique)) %>%
  mutate(references = lapply(references, sort)) %>%
  unnest(references) %>%
  distinct()

# Check for missing reference IDs
any(is.na(df$references))

# Convert IDs to numeric and join with CAS reference metadata
df$references <- as.numeric(df$references)
dff <- df %>% left_join(cas_reference, by = c("references" = "ref_id"))

# Remove records with missing publication years
dff_clean <- dff %>% filter(!is.na(ref_year))

# Calculate TAE per species
TAE_per_species <- dff_clean %>%
  group_by(valid_name) %>%
  summarise(
    taxonmic_effort = calculate_TAE(pick(everything()))
  )

# Inspect species with missing TAE values (if any)
TAE_per_species %>%
  filter(is.na(taxonmic_effort))

TAE_per_species %>%
  summarise(
    max_val = max(taxonmic_effort, na.rm = TRUE),
    min_val = min(taxonmic_effort, na.rm = TRUE),
    mean_val = mean(taxonmic_effort, na.rm = TRUE),
    median_val = median(taxonmic_effort)
  )

# Save the results
saveRDS(TAE_per_species, "input/data_prep/taxonmic_effort.rds")
