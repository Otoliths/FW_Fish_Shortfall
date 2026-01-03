library(robis)
library(dplyr)
library(readxl)
library(progressr)
library(purrr)
library(data.table)

handlers(global = TRUE)
handlers("progress")

# -------- Robust OBIS downloader --------
get_obis <- function(scientificname, save_path, ..., verbose = TRUE, sleep_sec = 5) {
  # 1) ensure output dir exists
  if (!dir.exists(save_path)) dir.create(save_path, recursive = TRUE, showWarnings = FALSE)
  
  # 2) sanitize filename: replace spaces and illegal chars
  safe_name <- gsub("\\s+", "_", scientificname)
  safe_name <- gsub("[^A-Za-z0-9_\\-\\.]", "_", safe_name)
  file_name <- file.path(save_path, paste0(safe_name, ".rds"))
  
  # 3) skip if already downloaded
  if (file.exists(file_name)) {
    if (verbose) cat("Already exists (skip):", scientificname, "\n")
    return(invisible(file_name))
  }
  
  # 4) fetch with safety
  if (verbose) cat("Downloading:", scientificname, "\n")
  res <- tryCatch(
    robis::occurrence(scientificname = scientificname, mof = TRUE, ...),
    error = function(e) {
      message("Error fetching ", scientificname, ": ", conditionMessage(e))
      return(NULL)
    }
  )
  
  # 5) handle empty result
  if (is.null(res) || nrow(res) == 0) {
    if (verbose) cat("No records returned:", scientificname, "\n")
    return(invisible(NULL))
  }
  
  # 6) save and sleep (be gentle to API)
  saveRDS(res, file = file_name)
  if (sleep_sec > 0) Sys.sleep(sleep_sec)
  invisible(file_name)
}

# -------- Load species list --------
cas <- readxl::read_excel("input/raw/cas_freshwater_v1.xlsx")

sp_list <- cas %>%
  transmute(valid_name = as.character(valid_name)) %>%
  filter(!is.na(valid_name), nzchar(valid_name)) %>%
  distinct() %>%
  pull(valid_name)

# -------- Iterate with a simple progress bar --------
save_path <- "input/raw/occ/obis"
n <- length(sp_list)
pb <- utils::txtProgressBar(min = 0, max = n, style = 3)

for (i in seq_len(n)) {
  try(
    get_obis(
      scientificname = sp_list[i],
      save_path = save_path,
      enddate = "2024-12-31",
      verbose = (i %% 50 == 0) # reduce console noise; print every ~50 species
    ),
    silent = TRUE
  )
  utils::setTxtProgressBar(pb, i)
}
close(pb)

cat("\nDone. Saved RDS files to:", normalizePath(save_path, mustWork = FALSE), "\n")

# --------------merge data-----------------------
cols_keep <- c("order","family","genus","species",
               "decimalLongitude","decimalLatitude","year",
               "taxonRank","occurrenceStatus","basisOfRecord")

read_select <- function(f) {
  x <- readRDS(f)
  keep <- intersect(names(x), cols_keep)
  if (length(keep) == 0) return(NULL)
  miss <- setdiff(cols_keep, names(x))
  if (length(miss)) x[miss] <- NA
  as_tibble(x[, cols_keep, drop = FALSE])
}

obis_files <- list.files(save_path, full.names = TRUE, pattern = "\\.rds$")

progressr::with_progress({
  p <- progressor(along = obis_files)
  r1 <- map(obis_files, ~{
    res <- read_select(.x)
    p(message = basename(.x))  
    res
  })
})

res <- data.table::rbindlist(r1, use.names = TRUE, fill = TRUE) %>%
  as_tibble()

#arrow::write_parquet(res, "input/raw/occ/obis/obis.parquet")