# ============================================================
# Taxonomic activity pipeline
# ------------------------------------------------------------
# Goal
# ----
# We estimate *taxonomic activity* following the logic of
# Moura & Jetz (2021, Nat. Ecol. Evol. 5:631–639).
#
# For each species, we want to quantify how "active" taxonomy
# was in its family at the time it was described. We do this by:
#
# 1. Parsing historical nomenclatural information (valid names,
#    synonyms, alternative combinations) from free-text
#    catalog records (CAS-style species accounts).
#
# 2. Extracting and normalizing authorship strings, including
#    multiple coauthors (e.g. "Smith & Jones", "Kottelat, Freyhof & Li").
#
# 3. For each family-year combination:
#      - Combine all unique taxonomists who described species
#        in that family in that year (across all species).
#      - Count how many distinct people that represents.
#      - Count how many unique species were described in that
#        family in that year.
#
# 4. Define taxonomic activity for that family-year as:
#        (# distinct taxonomists in that family/year) /
#        (# distinct species described in that family/year)
#
# This gives us something like "average number of active taxonomists
# per species" for that family and year — a proxy for research effort
# and naming intensity at the time of description.
#
# Output
# ------
# A data frame with one row per species, including:
# - all parsed authorships;
# - all contributing authors (cleaned and de-duplicated);
# - per-species author count;
# - per-family-year taxonomic activity;
# which we save to disk as an .rds file.
#
# Input data assumptions
# ----------------------
# We assume the source spreadsheet has at least:
#   cas_info          : long free-text CAS entry for the species
#   valid_name        : currently accepted species name
#   year_description  : publication year of the focal species
#   authorship        : canonical authorship string for the focal species
#                       (usually "Genus species Author 1979" or
#                        "Genus species (Author 1979)")
#   family            : taxonomic family of the focal species
#
# ============================================================

# -------------------------------
# Load required packages
# -------------------------------
library(dplyr)       # data wrangling pipes / mutate / group_by / summarise
library(stringr)     # regex extract / clean text
library(tidyr)       # unnesting / list-column helpers
library(purrr)       # map / pmap for list-column processing
library(openxlsx)    # read.xlsx for Excel import

# -------------------------------
# 1. Load input data
# -------------------------------
raw <- openxlsx::read.xlsx("input/raw/cas_freshwater_v1.xlsx")

# Keep only the columns we need downstream
df <- raw[, c("cas_info",
              "valid_name",
              "year_description",
              "authorship",
              "family")]

# -------------------------------
# 2. Helper to clean and collapse
#    extracted historical names
# -------------------------------
# str_extract_all() returns a vector of strings like
#   "Genus species Author 1979"
#   "Genus species (Author 1979)"
# We:
# - strip parentheses around authors/years;
# - trim whitespace;
# - drop duplicates;
# - join them back into a single semicolon-separated string.
collapse_names <- function(x) {
  if (length(x) == 0) return(NA_character_)
  x %>%
    gsub("[()]", "", .) %>%     # remove literal parentheses
    stringr::str_squish() %>%   # collapse multiple spaces
    unique() %>%                # drop repeated variants
    paste0(collapse = ";")       # join with ";"
}

# -------------------------------
# 3. Extract historical usage
#    ("Valid as ..." / "Synonym of ...")
# -------------------------------
# The CAS text often includes lines like:
#   •Valid as Aphyosemion etzeli (Berkenkamp 1979) --
#   •Synonym of Rhodeus sericeus (Pallas 1776) --
#   •Valid as Roloffia etzeli Berlenkamp 1979 --
#
# We capture:
#   <Genus> <species> <Author Year>
# allowing both "(Author 1979)" and "Author 1979".
#
# Notes on regex:
# - We use a fixed-length lookbehind for "•Valid as " or "•Synonym of "
#   because ICU regex in stringr requires look-behind with bounded length.
# - We then allow two patterns for the author+year block:
#      " (Author 1979)"  [with parentheses]
#      " Author 1979"    [without parentheses]
df <- df %>%
  mutate(
    past_valid_names_raw = stringr::str_extract_all(
      cas_info,
      "(?<=•Valid as )[A-Z][a-z]+\\s[a-z]+(?:\\s\\([^)]*?\\d{4}\\)|[^()]*?\\d{4})"
    ),
    past_synonym_names_raw = stringr::str_extract_all(
      cas_info,
      "(?<=•Synonym of )[A-Z][a-z]+\\s[a-z]+(?:\\s\\([^)]*?\\d{4}\\)|[^()]*?\\d{4})"
    )
  ) %>%
  mutate(
    # Clean text and collapse to a single string per row
    past_valid_names   = purrr::map_chr(past_valid_names_raw,   collapse_names),
    past_synonym_names = purrr::map_chr(past_synonym_names_raw, collapse_names)
  ) %>%
  select(-past_valid_names_raw, -past_synonym_names_raw)

# After this step, for a row we might have:
#   past_valid_names   = "Aphyosemion etzeli Berkenkamp 1979;Roloffia etzeli Berlenkamp 1979"
#   past_synonym_names = "Aphyosemion geryi Lambert 1958;Scriptaphyosemion roloffi Roloff 1936"

# -------------------------------
# 4. Build a unified "authorships"
#    string for each focal species
# -------------------------------
# We merge:
#   - the focal authorship (authorship column),
#   - all past "Valid as ..." names,
#   - all "Synonym of ..." names,
# into a single semicolon-separated string.
#
# We also:
# - drop NAs / empty strings;
# - trim whitespace;
# - deduplicate.
df <- df %>%
  mutate(
    authorships = purrr::pmap_chr(
      list(authorship, past_valid_names, past_synonym_names),
      ~ c(..1, ..2, ..3) %>%                # collect the 3 sources
        unlist() %>%                        # flatten
        stats::na.omit() %>%                # drop NA values
        stringr::str_trim() %>%             # trim leading/trailing whitespace
        purrr::discard(~ .x == "") %>%      # drop empty strings
        unique() %>%                        # ensure uniqueness
        paste0(collapse = ";")               # join with ";"
    )
  )
# We intentionally keep "authorship", "past_valid_names",
# and "past_synonym_names" in df for traceability. You can drop later.

# Example authorships for one species might look like:
# "Pampus candidus (Cuvier 1829);Pampus candidus Cuvier 1829;Pampus argenteus Euphrasen 1788"

# -------------------------------
# 5. Parse authorships to extract
#    individual human authors
# -------------------------------
# Goal:
# - For each focal species (row),
#   turn "authorships" into a set of unique author names.
#
# Why we need this:
# - The same species can have multiple historical combinations,
#   sometimes with multiple coauthors.
# - We want a cleaned list of *people*, not "Genus species ...".
#
# Steps:
#   (a) Split authorships by ";" into "taxon strings"
#       e.g. "Pampus candidus (Cuvier 1829)"
#   (b) Extract the author block between the species epithet
#       and the 4-digit year.
#       This handles:
#          "Genus species (Author 1979)"
#          "Genus species Author 1979"
#          "Genus species Smith, Jones & Brown 2001"
#   (c) Split multi-author strings on ",", "&", or " and "
#       to get individual names.
#   (d) Deduplicate within the row.
#
# We'll store:
#   - "authors": human-readable list, like "Cuvier; Euphrasen"
#   - "authors_count_species": how many unique authors worked
#                              on this focal species (across all combos)
df <- df %>%
  mutate(
    authors_info = purrr::map(authorships, function(s) {
      # If no authorships, return empty
      if (is.na(s) || s == "") {
        return(list(authors_vec = character(0)))
      }
      
      # (a) split into individual "Genus species Author Year" entries
      items <- stringr::str_split(s, ";\\s*")[[1]]
      items <- stringr::str_squish(items)
      
      # (b) pull out the "Author ... Year" block
      # Regex explanation:
      #   ^[A-Z][a-z]+\\s[a-z]+     = Genus species
      #   \\s*\\(?                  = optional "(" after a space
      #   ([^()]*?)\\s\\d{4}        = capture authors until the year
      #   \\)?$                     = optional ")"
      #
      # The capture group ([^()]*?) is what we want:
      # e.g. "Cuvier", "Euphrasen",
      #      "Kottelat & Freyhof", "Smith, Jones & Brown"
      m <- stringr::str_match(
        items,
        "^[A-Z][a-z]+\\s[a-z]+\\s*\\(?([^()]*?)\\s\\d{4}\\)?$"
      )
      
      raw_authors <- m[, 2]                     # take captured group 1
      raw_authors <- raw_authors[!is.na(raw_authors)]
      raw_authors <- stringr::str_squish(raw_authors)
      
      # (c) split on ',', '&', or ' and ' to isolate individuals
      indiv_authors <- unlist(
        stringr::str_split(
          raw_authors,
          "\\s*(?:,|&| and )\\s*",
          simplify = FALSE
        )
      )
      
      indiv_authors <- stringr::str_squish(indiv_authors)
      indiv_authors <- indiv_authors[indiv_authors != ""]
      indiv_authors <- unique(indiv_authors)
      
      # Return as a list so we keep it structured per-row
      list(authors_vec = indiv_authors)
    })
  ) %>%
  mutate(
    # Flatten the authors_vec list-column into readable strings
    authors = purrr::map_chr(
      authors_info,
      ~ if (length(.x$authors_vec) == 0) NA_character_
      else paste0(.x$authors_vec, collapse = ";")),
    
    # Count distinct authors for this focal species
    authors_count_species = purrr::map_int(
      authors_info,
      ~ length(.x$authors_vec)
    )) %>% 
  dplyr::select(-authors_info)

# After this step:
#   authors might be "Cuvier; Euphrasen"
#   authors_count_species might be 2

# -------------------------------
# 6. Compute family-year level
#    "taxonomic activity"
# -------------------------------
# We now move from species-level to family-year-level.
#
# For each combination of (family, year_description):
#   - Gather ALL authors from ALL species in that family-year.
#   - Take the unique set of authors.
#   - Count how many unique authors were active.
#   - Count how many distinct species were described.
#
# Then we define:
#   taxonomic_activity_family_year =
#       (# unique authors in that family-year) /
#       (# species described in that family-year)
#
# This matches the verbal definition: "number of unique taxonomists
# who described species within the same family and year,
# standardized by the number of species in that family-year."
family_year_stats <- df %>%
  filter(!is.na(family), !is.na(year_description)) %>%
  group_by(family, year_description) %>%
  summarize(
    species_count = n_distinct(valid_name),
    unique_authors_family_year = list({
      all_authors <- unlist(str_split(authors[!is.na(authors)], ";\\s*"))
      all_authors <- str_squish(all_authors)
      unique(all_authors[all_authors != ""])
    }),
    n_authors_family_year = lengths(unique_authors_family_year),
    taxonomic_activity_family_year = n_authors_family_year / species_count,
    .groups = "drop"
  )



# Explanation:
#   species_count               = how many different species were named
#                                 in this family in this year
#   n_authors_family_year       = how many unique people contributed
#                                 to naming those species
#   taxonomic_activity_family_year
#                               = average #taxonomists per species

# -------------------------------
# 7. Attach static TAA back to each species
# -------------------------------
taxonomic_activity <- df %>%
  left_join(
    family_year_stats %>%
      select(
        family,
        year_description,
        species_count,
        n_authors_family_year,
        taxonomic_activity_family_year
      ),
    by = c("family", "year_description")
  ) %>%
  mutate(
    # Keep per-species author count (optional diagnostic)
    authors_count = authors_count_species,
    
    # This is the final per-species taxonomic activity metric,
    # i.e. the activity level of the species' family-year context.
    taxonomic_activity = taxonomic_activity_family_year
  )

# Optional QA: species with missing static TAA
taxonomic_activity %>%
  filter(is.na(taxonomic_activity))
# ==================================================
# 8. Time-aware but non-cumulative TAA function
# ==================================================
# Logic:
# - TAA is NOT cumulative
# - It is a temporally localized indicator
# - TAA(t) = the most recent available family-year activity up to t
# - Therefore, TAA only changes when a new family-year activity appears

calculate_TAA <- function(df_species,
                          family_year_stats,
                          t_eval = NULL) {
  
  # -----------------------------
  # 0. Basic input checks
  # -----------------------------
  if (is.null(df_species) || !is.data.frame(df_species) || nrow(df_species) == 0) {
    return(NA_real_)
  }
  
  if (is.null(family_year_stats) || !is.data.frame(family_year_stats) || nrow(family_year_stats) == 0) {
    return(NA_real_)
  }
  
  if (!"family" %in% names(df_species)) {
    warning("Column 'family' not found in df_species.")
    return(NA_real_)
  }
  
  required_cols <- c("family", "year_description", "taxonomic_activity_family_year")
  if (!all(required_cols %in% names(family_year_stats))) {
    warning("Required columns missing in family_year_stats.")
    return(NA_real_)
  }
  
  # -----------------------------
  # 1. Extract focal family
  # -----------------------------
  fam <- unique(df_species$family[!is.na(df_species$family)])
  if (length(fam) == 0) {
    return(NA_real_)
  }
  fam <- fam[1]
  
  # -----------------------------
  # 2. Prepare family-level time series
  # -----------------------------
  fy <- family_year_stats %>%
    dplyr::filter(.data$family == fam, !is.na(.data$year_description)) %>%
    dplyr::mutate(
      year_description = suppressWarnings(as.integer(.data$year_description)),
      taxonomic_activity_family_year = as.numeric(.data$taxonomic_activity_family_year)
    ) %>%
    dplyr::filter(
      !is.na(.data$year_description),
      !is.na(.data$taxonomic_activity_family_year)
    ) %>%
    dplyr::arrange(.data$year_description) %>%
    dplyr::distinct(.data$family, .data$year_description, .keep_all = TRUE)
  
  if (nrow(fy) == 0) {
    return(NA_real_)
  }
  
  # -----------------------------
  # 3. Define evaluation year
  # -----------------------------
  # Default: use the species' own description year
  if (is.null(t_eval)) {
    if (!"year_description" %in% names(df_species)) {
      warning("Column 'year_description' not found in df_species.")
      return(NA_real_)
    }
    
    t_eval <- unique(df_species$year_description[!is.na(df_species$year_description)])
    if (length(t_eval) == 0) {
      return(NA_real_)
    }
    t_eval <- t_eval[1]
  }
  
  t_eval <- suppressWarnings(as.integer(t_eval))
  if (!is.finite(t_eval)) {
    warning("Invalid 't_eval'.")
    return(NA_real_)
  }
  
  # -----------------------------
  # 4. Temporal correction
  # -----------------------------
  t_min_all <- min(fy$year_description, na.rm = TRUE)
  
  # Before the first available family-year activity, TAA is undefined
  if (t_eval < t_min_all) {
    return(NA_real_)
  }
  
  # Keep all family-year records up to t_eval
  fy_sub <- fy[fy$year_description <= t_eval, , drop = FALSE]
  if (nrow(fy_sub) == 0) {
    return(NA_real_)
  }
  
  # Use the most recent available family-year activity up to t_eval
  taa_value <- fy_sub$taxonomic_activity_family_year[nrow(fy_sub)]
  
  return(as.numeric(taa_value))
}

calculate_TAA(
  df_species = df %>% dplyr::filter(valid_name == "Amatitlania kanna"),
  family_year_stats = family_year_stats
)
calculate_TAA(
  df_species = df %>% dplyr::filter(valid_name == "Amatitlania kanna"),
  family_year_stats = family_year_stats,
  t_eval = 2011
)


# ==================================================
# Batch compute yearly TAA for each species: 1758:2024
# ==================================================
# species_nested <- df %>%
#   dplyr::group_by(valid_name) %>%
#   tidyr::nest()
# 
# year_seq <- 1758:2024
# TAA_yearly_list <- vector("list", nrow(species_nested))
# 
# for (i in seq_len(nrow(species_nested))) {
#   sp <- species_nested$valid_name[i]
#   sp_df <- species_nested$data[[i]]
#   
#   cat(sprintf("[%d/%d] Computing TAA for %s\n", i, nrow(species_nested), sp))
#   
#   taa_vec <- vapply(
#     year_seq,
#     function(y) calculate_TAA(
#       df_species = sp_df,
#       family_year_stats = family_year_stats,
#       t_eval = y
#     ),
#     numeric(1)
#   )
#   
#   TAA_yearly_list[[i]] <- data.frame(
#     valid_name = sp,
#     year = year_seq,
#     TAA = taa_vec,
#     stringsAsFactors = FALSE
#   )
# }
# 
# TAA_yearly <- dplyr::bind_rows(TAA_yearly_list)
# 
# saveRDS(TAA_yearly, "input/data_prep/taxonomic_activity_year.rds")


################################################################################
build_family_TAA_lookup <- function(family_year_stats,
                                    year_seq = 1758:2024) {
  required_cols <- c("family", "year_description", "taxonomic_activity_family_year")
  
  fy <- family_year_stats[, required_cols, drop = FALSE]
  fy$year_description <- suppressWarnings(as.integer(fy$year_description))
  fy$taxonomic_activity_family_year <- suppressWarnings(as.numeric(fy$taxonomic_activity_family_year))
  
  fy <- fy[!is.na(fy$family) &
             !is.na(fy$year_description) &
             !is.na(fy$taxonomic_activity_family_year), , drop = FALSE]
  
  fy <- fy[order(fy$family, fy$year_description), , drop = FALSE]
  fy <- fy[!duplicated(fy[c("family", "year_description")]), , drop = FALSE]
  
  family_list <- split(fy, fy$family)
  out_list <- vector("list", length(family_list))
  fam_names <- names(family_list)
  
  for (i in seq_along(family_list)) {
    fam <- fam_names[i]
    dat <- family_list[[i]]
    
    idx <- findInterval(year_seq, dat$year_description)
    taa_vec <- rep(NA_real_, length(year_seq))
    ok <- idx > 0
    taa_vec[ok] <- dat$taxonomic_activity_family_year[idx[ok]]
    
    out_list[[i]] <- data.frame(
      family = fam,
      year = year_seq,
      TAA = taa_vec,
      stringsAsFactors = FALSE
    )
  }
  
  dplyr::bind_rows(out_list)
}

year_seq <- 1758:2024

family_TAA_lookup <- build_family_TAA_lookup(
  family_year_stats = family_year_stats,
  year_seq = year_seq
)

species_family <- df %>%
  dplyr::filter(!is.na(valid_name), !is.na(family)) %>%
  dplyr::group_by(valid_name) %>%
  dplyr::summarise(family = dplyr::first(family), .groups = "drop")

TAA_yearly_df <- species_family %>%
  dplyr::left_join(
    family_TAA_lookup,
    by = "family",
    relationship = "many-to-many"
  ) %>%
  dplyr::select(valid_name, year, TAA) %>%
  dplyr::arrange(valid_name, year)

saveRDS(TAA_yearly_df, "input/data_prep/taxonomic_activity_year.rds")
library(ggplot2)
ggplot(data = TAA_yearly_df,aes(x= year, y = TAA,group = as.factor(valid_name)))+
  geom_line(na.rm = T,linewidth = 0.2)


ggplot(data = family_year_stats[which(family_year_stats$family %in% "Acestrorhamphidae"),],aes(x= year_description, y = taxonomic_activity_family_year,group = as.factor(family)))+
  geom_line()

ggplot(data = family_year_stats,aes(x= year_description, y = taxonomic_activity_family_year,group = as.factor(family)))+
  geom_line(linewidth = 0.2)+
  scale_y_continuous(limits = c(0,3))


# ==================================================
# QA checks
# ==================================================
# Preview
head(TAA_yearly)

# Missingness
TAA_yearly %>%
  summarise(
    n_total = n(),
    n_na = sum(is.na(TAA)),
    prop_na = mean(is.na(TAA))
  )

# Summary stats
TAA_yearly %>%
  summarise(
    max_val = max(TAA, na.rm = TRUE),
    min_val = min(TAA, na.rm = TRUE),
    mean_val = mean(TAA, na.rm = TRUE),
    median_val = median(TAA, na.rm = TRUE)
  )

# Example: inspect one species trajectory
TAA_yearly %>%
  filter(valid_name == species_names[1]) %>%
  head(30)

# ==================================================
# 9. Save outputs
# ==================================================
# dir.create("output", showWarnings = FALSE, recursive = TRUE)
# 
# saveRDS(family_year_stats, "output/family_year_stats.rds")
# saveRDS(taxonomic_activity, "output/taxonomic_activity_static_by_species.rds")
# saveRDS(TAA_per_species, "output/TAA_per_species.rds")
saveRDS(TAA_yearly, "input/data_prep/taxonomic_activity_yearly.rds")

write.csv(TAA_per_species, "output/TAA_per_species.csv", row.names = FALSE)
write.csv(TAA_yearly, "output/TAA_yearly_1758_2024.csv", row.names = FALSE)


# -------------------------------
# 9. Save result
# -------------------------------
saveRDS(taxonomic_activity[,c(2,15)],"input/data_prep/taxonomic_activity.rds")

# End of pipeline
