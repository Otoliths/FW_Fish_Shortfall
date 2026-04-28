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
cas_reference <- openxlsx::read.xlsx("input/raw/cas_reference.xlsx")
df <- openxlsx::read.xlsx("input/raw/cas_freshwater_v1.xlsx")

# Prepare reference IDs: split, deduplicate, sort, unnest, remove duplicates
dff_clean <- df %>%
  select(valid_name, authorship, references) %>%
  mutate(references = str_split(references, ",\\s*")) %>%
  mutate(references = lapply(references, unique)) %>%
  mutate(references = lapply(references, sort)) %>%
  unnest(references) %>%
  distinct() %>%
  mutate(references = suppressWarnings(as.numeric(references))) %>%
  left_join(cas_reference, by = c("references" = "ref_id")) %>%
  filter(!is.na(ref_year))


calculate_TAE_series <- function(df,
                                 year_seq = 1758:2024,
                                 lambda = log(2) / 20,
                                 author_sep = ",") {
  
  # -----------------------------
  # 0. Basic input checks
  # -----------------------------
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
    return(rep(NA_real_, length(year_seq)))
  }
  
  if (!"ref_year" %in% names(df)) {
    warning("Column 'ref_year' not found.")
    return(rep(NA_real_, length(year_seq)))
  }
  
  if (!"ref_authorship" %in% names(df)) {
    df$ref_authorship <- NA_character_
  }
  
  # -----------------------------
  # 1. Keep only valid reference years
  # -----------------------------
  df2 <- df[!is.na(df$ref_year), , drop = FALSE]
  if (nrow(df2) == 0) return(rep(NA_real_, length(year_seq)))
  
  df2$ref_year <- suppressWarnings(as.integer(df2$ref_year))
  df2 <- df2[!is.na(df2$ref_year), , drop = FALSE]
  if (nrow(df2) == 0) return(rep(NA_real_, length(year_seq)))
  
  # -----------------------------
  # 2. Optional cleaning
  # -----------------------------
  if (exists("clean_refs", mode = "function")) {
    df2 <- clean_refs(df2, "ref_authorship")
  }
  
  df2$ref_authorship[is.na(df2$ref_authorship)] <- ""
  
  # -----------------------------
  # 3. Precompute Ai once
  # -----------------------------
  split_auth <- strsplit(
    df2$ref_authorship,
    paste0("\\s*", author_sep, "\\s*"),
    perl = TRUE
  )
  
  Ai <- vapply(split_auth, function(x) {
    x <- trimws(x)
    x <- x[nzchar(x)]
    length(unique(x))
  }, integer(1))
  
  Ai[is.na(Ai) | Ai < 0] <- 0L
  
  # log(1 + Ai), reused later
  logAi <- log1p(Ai)
  ref_year <- df2$ref_year
  
  # sort by year for cumulative inclusion
  ord <- order(ref_year)
  ref_year <- ref_year[ord]
  logAi <- logAi[ord]
  
  n_ref <- length(ref_year)
  n_year <- length(year_seq)
  
  out <- rep(NA_real_, n_year)
  
  # years before first reference remain NA
  first_ref_year <- ref_year[1]
  
  valid_idx <- which(year_seq >= first_ref_year)
  if (length(valid_idx) == 0) return(out)
  
  # -----------------------------
  # 4. Build year x reference age matrix
  #    only for years >= first reference year
  # -----------------------------
  years_valid <- year_seq[valid_idx]
  
  # age_mat[r, c] = years_valid[r] - ref_year[c]
  age_mat <- outer(years_valid, ref_year, "-")
  
  # only references already available by year t
  included_mat <- age_mat >= 0
  
  # exponential decay
  weight_mat <- exp(-lambda * pmax(age_mat, 0))
  
  # multiply by log(1 + Ai), and set future refs to 0
  contrib_mat <- weight_mat * 
    matrix(logAi, nrow = nrow(weight_mat), ncol = ncol(weight_mat), byrow = TRUE) *
    included_mat
  
  # cumulative number of references available by each year
  N_t <- rowSums(included_mat)
  
  # earliest available reference year at each t is fixed = first_ref_year
  D_t <- years_valid - first_ref_year + 1L
  
  # final TAE
  out_valid <- sqrt(D_t) * (rowSums(contrib_mat) / N_t)
  
  out[valid_idx] <- out_valid
  out
}

species_nested <- dff_clean %>%
  dplyr::group_by(valid_name) %>%
  tidyr::nest()

year_seq <- 1758:2024

TAE_yearly_list <- vector("list", nrow(species_nested))

for (i in seq_len(nrow(species_nested))) {
  sp <- species_nested$valid_name[i]
  sp_df <- species_nested$data[[i]]
  
  cat(sprintf("[%d/%d] Computing TAE for %s\n", i, nrow(species_nested), sp))
  
  tae_vec <- calculate_TAE_series(
    df = sp_df,
    year_seq = year_seq
  )
  
  TAE_yearly_list[[i]] <- data.frame(
    valid_name = sp,
    year = year_seq,
    TAE = tae_vec,
    stringsAsFactors = FALSE
  )
}

TAE_yearly_df <- dplyr::bind_rows(TAE_yearly_list)
saveRDS(TAE_yearly_df, "input/data_prep/taxonomic_effort_year.rds")
