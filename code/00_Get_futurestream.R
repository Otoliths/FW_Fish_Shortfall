# Bosmans, J., Wanders, N., Bierkens, M.F.P. et al. FutureStreams, a global dataset
# of future streamflow and water temperature. Sci Data 9, 307 (2022). https://doi.org/10.1038/s41597-022-01410-6

#Bosmans, J. et al. FutureStreams [Data set]. Yoda. https://doi.org/10.24416/UU01-T7TVTQ (2022).


# ====================== User config ======================
out_base   <- normalizePath(file.path(getwd(),"input/raw"), mustWork = FALSE)

# ==================== URL specification ==================
base_host <- "https://geo.public.data.uu.nl"
vault     <- "vault-futurestreams/research-futurestreams%5B1633685642%5D/original"

variables <- c("discharge", "waterTemp")
models    <- c("gfdl", "hadgem", "ipsl", "miroc", "noresm")
ranges_hist   <- list(c("1976-01-07","1985-12-30"),
                      c("1986-01-07","1995-12-30"),
                      c("1996-01-07","2005-12-30"))
ranges_rcp2p6 <- list(c("2006-01-07","2019-12-30"),
                      c("2020-01-07","2029-12-30"))
scenarios <- list(hist = ranges_hist, rcp2p6 = ranges_rcp2p6)

# ======================== Helpers ========================
# Ensure a directory exists and return its normalized absolute path
ensure_dir <- function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE)
  normalizePath(p, mustWork = FALSE)
}

# Compose a single file URL following the dataset naming convention
# Pattern: .../{var}/{scen}/{model}/{var}_weekAvg_output_{model}_{scen}_{start}_to_{end}.nc
build_url <- function(var, scen, model, start, end) {
  dir_path <- sprintf("%s/%s/%s/%s/%s", base_host, vault, var, scen, model)
  file_nm  <- sprintf("%s_weekAvg_output_%s_%s_%s_to_%s.nc", var, model, scen, start, end)
  paste0(dir_path, "/", file_nm)
}

# Build a download manifest (one row per file)
build_table <- function() {
  rows <- list()
  for (var in variables) {
    for (scen in names(scenarios)) {
      for (model in models) {
        for (rg in scenarios[[scen]]) {
          url <- build_url(var, scen, model, rg[1], rg[2])
          # Save by variable: out_base/<variable>/<filename>.nc
          dest_dir  <- file.path(out_base, var)
          dest_file <- file.path(dest_dir, basename(url))
          rows[[length(rows) + 1L]] <- data.frame(
            variable = var,
            scenario = scen,
            model    = model,
            start    = rg[1],
            end      = rg[2],
            url      = url,
            dest_dir = dest_dir,
            dest_file = dest_file,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  do.call(rbind, rows)
}

# Locate wget: prefer system PATH; fallback to ./exec/wget(.exe)
find_wget <- function() {
  w <- Sys.which("wget")
  if (nzchar(w)) return(w)
  cand <- file.path(getwd(), "code/exec", if (.Platform$OS.type == "windows") "wget.exe" else "wget")
  if (file.exists(cand)) return(cand)
  ""
}

# Download a set of files sequentially with wget (resumable, robust)
download_with_wget <- function(df, rate_limit = NULL, force = FALSE) {
  # df must contain columns: url, dest_dir, dest_file, variable
  wget_bin <- find_wget()
  if (!nzchar(wget_bin)) {
    stop("wget not found. Put wget(.exe) under ./exec/ or add it to your system PATH.")
  }
  
  n <- nrow(df)
  if (n == 0L) {
    message("Nothing to download: manifest has zero rows.")
    return(invisible(list(ok = 0L, skip = 0L, fail = 0L)))
  }
  
  ok <- 0L; skip <- 0L; fail <- 0L
  t0 <- Sys.time()
  
  for (i in seq_len(n)) {
    u     <- df$url[i]
    ddir  <- df$dest_dir[i]
    dfile <- df$dest_file[i]
    
    ensure_dir(ddir)
    
    # Skip existing non-empty files unless force = TRUE
    if (!force && file.exists(dfile) && isTRUE(file.info(dfile)$size > 0)) {
      skip <- skip + 1L
      message(sprintf("[%d/%d] Skip (exists): %s", i, n, basename(dfile)))
      next
    }
    
    # Build wget args: resumable (-c), retries, timeouts, output name
    args <- c(
      "-c",
      "--tries=10",
      "--timeout=60",
      "--waitretry=20",
      shQuote(u),
      "-O", shQuote(dfile)
    )
    if (!is.null(rate_limit)) {
      # e.g., "2m" for 2 MB/s; "500k" for 500 KB/s
      args <- c("--limit-rate", rate_limit, args)
    }
    
    message(sprintf("→ [%d/%d] Start: %s\n   Dir : %s", i, n, basename(dfile), ddir))
    code <- suppressWarnings(system2(wget_bin, args = args))
    
    # Validate success by exit code and non-zero file size
    if (identical(code, 0L) && file.exists(dfile) && isTRUE(file.info(dfile)$size > 0)) {
      ok <- ok + 1L
      message(sprintf("✓ Done: %s", basename(dfile)))
    } else {
      fail <- fail + 1L
      warning(sprintf("× Failed: %s (exit=%s)", basename(dfile), as.character(code)))
    }
  }
  
  dt <- difftime(Sys.time(), t0, units = "mins")
  message(sprintf("\n[wget summary] ok=%d, skip=%d, fail=%d, total=%d, elapsed=%.1f min\n",
                  ok, skip, fail, n, as.numeric(dt)))
  invisible(list(ok = ok, skip = skip, fail = fail, total = n))
}

# ========================= Main ==========================
# Ensure base output directory exists
ensure_dir(out_base)

# Build manifest
url_tab <- build_table()

download_with_wget(url_tab, rate_limit = NULL, force = FALSE)

# Show a quick tree preview (shallow)
cat("\n== Output preview ==\n")
print(utils::head(list.files(out_base, recursive = TRUE), 20))
cat("\nDone. You can re-run safely for resume; existing files are skipped by default.\n")