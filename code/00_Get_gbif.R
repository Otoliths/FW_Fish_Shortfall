# ================================================================
#  GBIF Occurrence Download Script
#  
#  Description:
#  - Reads a species list from an Excel file
#  - Retrieves taxon keys for each species from GBIF
#  - Submits a single occurrence download request for all species
#  - Saves the download request key as an RDS file for later use
#  
#  Note:
#  - You need a registered GBIF account to submit download requests.
#  - Sign up at: https://www.gbif.org/user/profile
#  - Retrieving taxon keys may take ~2 hours depending on list size.
#  - Large occurrence downloads (e.g., >30M records) may take several hours to complete.
#
#  - Author: Dr. Liuyong Ding
#  - Institute: Institute of Hydrobiology, Chinese Academy of Sciences
#  - Contact: ly_ding@126.com
#  - Date: 2025-10-08
# ================================================================

# === Load required libraries ===
library(rgbif)     # Interface to the GBIF API for taxon keys and occurrence downloads
library(readxl)    # Read species list from Excel files (.xlsx)
library(stringr)   # String manipulation (e.g., cleaning species names)
library(dplyr)     # Data wrangling and tidy data operations

# === Read species list ===
# The Excel file should contain a column named 'valid_name'
cas <- read_excel("input/raw/cas_freshwater_v1.xlsx")
sp_list <- unique(cas$valid_name)   # get unique species names

# Preview the first few species
head(sp_list)
# [1] "Amatitlania kanna"        "Atherinella colombiensis" "Gobiesox juradoensis"    
# [4] "Toxotes microlepis"       "Pampus candidus"          "Amblydoras nheco"


# === Prefer a species-level GBIF key.===
# Strategy:
# 1) name_backbone(): if rank == "SPECIES" use usageKey
# 2) else use speciesKey if present (handles subspecies/varieties mapping to their species)
# 3) if SYNONYM and acceptedUsageKey exists, use acceptedUsageKey
# 4) fallback: name_suggest(rank="SPECIES") and try exact canonical match
# 5) drop invalid keys (NA or 1)

get_taxon_key <- function(sp, min_conf = 80) {
  sp_clean <- str_squish(sp)
  
  bb <- tryCatch(name_backbone(name = sp_clean, verbose = FALSE),
                 error = function(e) NULL)
  if (!is.null(bb) && !is.null(bb$confidence) && bb$confidence < min_conf) {
    bb <- NULL
  }
  
  # 1) species rank direct
  if (!is.null(bb) && identical(bb$rank, "SPECIES") && !is.null(bb$usageKey)) {
    key <- as.integer(bb$usageKey)
    if (!is.na(key) && key != 1L) return(key)
  }
  
  # 2) speciesKey
  if (!is.null(bb) && !is.null(bb$speciesKey)) {
    key <- as.integer(bb$speciesKey)
    if (!is.na(key) && key != 1L) return(key)
  }
  
  # 3) acceptedUsageKey
  acc_key <- NULL
  if (!is.null(bb)) {
    acc_key <- bb$acceptedUsageKey
    if (is.null(acc_key)) acc_key <- bb$acceptedKey
  }
  if (!is.null(acc_key)) {
    key <- as.integer(acc_key)
    if (!is.na(key) && key != 1L) return(key)
  }
  
  # 4) fallback to name_suggest (add safe checks here)
  sug <- tryCatch(name_suggest(q = sp_clean, rank = "SPECIES", limit = 10),
                  error = function(e) NULL)
  if (!is.null(sug) && is.data.frame(sug) && nrow(sug) > 0) {
    # exact canonical match
    if ("canonicalName" %in% names(sug)) {
      hit_idx <- which(tolower(sug$canonicalName) == tolower(sp_clean))[1]
      if (!is.na(hit_idx)) {
        hit <- sug[hit_idx, ]
        if (!is.na(hit$key) && hit$key != 1L) return(as.integer(hit$key))
      }
    }
    # top species suggestion
    hit_idx <- which(sug$rank == "SPECIES")[1]
    if (!is.na(hit_idx)) {
      hit <- sug[hit_idx, ]
      if (!is.na(hit$key) && hit$key != 1L) return(as.integer(hit$key))
    }
  }
  
  return(NA_integer_)
}

# === Get GBIF taxon keys for species list ===
message("⏳ Retrieving taxon keys for all species... ",
        "\nThis may take up to 2 hours depending on network speed and GBIF API response time.")

# Apply get_taxon_key() to each species in sp_list
tk_group <- vapply(sp_list, get_taxon_key, integer(1))

# Remove missing or invalid keys (NA)
tk_group <- tk_group[!is.na(tk_group)]

# Identify species not matched in the first round
unmatched_species <- setdiff(sp_list, names(tk_group))

# Second round: try to get taxon keys for unmatched species
tk_res <- vapply(unmatched_species, get_taxon_key, integer(1))
tk_res <- tk_res[!is.na(tk_res)]

# Combine both sets of taxon keys
tk_total <- c(tk_group, tk_res)

#  Download GBIF occurrence data ===
key <- occ_download(
  pred_in("taxonKey", tk_total),          # select taxon keys
  pred("hasCoordinate", TRUE),            # only records with coordinates
  pred_lte("year", 2024),                 # year <= 2024
  format = "DWCA",                        # Darwin Core Archive format
  user = "xxx",                           # GBIF username
  pwd = "xxx",                            # GBIF password
  email = "xxxxxx"                        # GBIF registered email
)

# Save download key object
saveRDS(key, file = "input/raw/occ/gbif/gbif_occ_freshwater_fish.rds")

# === Load GBIF download key object ===
key_obj <- readRDS("input/raw/occ/gbif/gbif_occ_freshwater_fish.rds")

# === Friendly reminder ===
message("⏳ GBIF is processing a large dataset (~36 million records, ~8 GB). ",
        "\nThe download may take up to 3 hours to complete. ",
        "\nPlease check the status periodically with **occ_download_meta()**.")

# === Retrieve download metadata from GBIF ===
meta <- occ_download_meta(key_obj)

print(occ_download_meta(key_obj))
# <<gbif download metadata>>
# Status: SUCCEEDED
# DOI: 10.15468/dl.dzfcsw
# Format: DWCA
# Download key: 0060713-250920141307145
# Created: 2025-10-08T07:22:08.671+00:00
# Modified: 2025-10-08T10:31:28.966+00:00
# Download link: https://api.gbif.org/v1/occurrence/download/request/0060713-250920141307145.zip
# Total records: 36302470

# === Summarize key metadata fields in a tidy data frame ===
download_info <- data.frame(
  n_species      = length(tk_total),         # number of target species
  Status         = meta$status,              # download job status (e.g. RUNNING, SUCCEEDED)
  DOI            = meta$doi,                 # GBIF DOI (generated after completion)
  Format         = meta$request$format,      # file format (e.g. DWCA)
  Key            = meta$key,                 # unique download key
  Created        = meta$created,             # creation timestamp
  Modified       = meta$modified,            # last modified timestamp
  Total_records  = meta$totalRecords,        # total number of records
  downloadLink   = meta$downloadLink,        # direct download link (ZIP file)
  stringsAsFactors = FALSE
)

cat(
  "✅ GBIF Download Summary\n",
  "Species count : ", download_info$n_species, "\n",
  "Status        : ", download_info$Status, "\n",
  "DOI           : ", download_info$DOI, "\n",
  "Format        : ", download_info$Format, "\n",
  "Key           : ", download_info$Key, "\n",
  "Created       : ", download_info$Created, "\n",
  "Modified      : ", download_info$Modified, "\n",
  "Total records : ", format(download_info$Total_records, big.mark = ","), "\n",
  "Download URL  :\n", download_info$downloadLink, "\n"
)
