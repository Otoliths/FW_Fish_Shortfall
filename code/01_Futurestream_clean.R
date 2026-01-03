#https://github.com/JoyceBosmans/FutureStreams/blob/main/Figures/Plot_derived_vars.R

# =====================================================================
# Annual Water Temperature/discharge Aggregation with terra::app()
# - Computes per-year, per-pixel means from weekly NetCDF stacks
# - Supports single-model annual means and multi-model ensemble means
# - File naming assumed: waterTemp_weekAvg_output_{model}_{scen}_{start}_to_{end}.nc
#   e.g., waterTemp_weekAvg_output_gfdl_hist_1976-01-07_to_1985-12-30.nc
# =====================================================================

# =====================================================================
# compute_model_annual_means (with unified threshold logic for WT / Q)
# - Reads weekly NetCDF/TIF stacks for a model (helper .read_model_weekly)
# - Masks weekly outliers IN-PLACE according to 'var' and 'threshold':
#     * if var == "waterTemp": weekly values > threshold -> NA
#     * if var == "discharge": weekly values < threshold -> NA
#   (If threshold is NULL, defaults: waterTemp -> 350 K, discharge -> 10 m3/s)
# - Computes per-pixel annual mean from the masked weekly layers (using terra::app)
# - Writes one single-layer GeoTIFF per year immediately to disk (or returns rasters)
# - English comments throughout
# =====================================================================

library(terra)
library(lubridate)
library(stringr)

# ---------------------------------------------------------------------
# Helper: read & concatenate weekly chunks for a model into a SpatRaster time-series
# - Assumes filenames like: {var}_weekAvg_output_{model}_{scen}_{YYYY-MM-DD}_to_{YYYY-MM-DD}.nc
# - Returns list(r = SpatRaster (layers = weeks), dates = Date vector for each layer)
# ---------------------------------------------------------------------


# ---- Helper: derive weekly dates either from time(r) or from filename ----
.get_dates <- function(r, f) {
  # Prefer native time stamps if present and aligned
  tt <- time(r)
  if (!is.null(tt) && length(tt) >= nlyr(r)) {
    return(as.Date(tt)[seq_len(nlyr(r))])
  }
  # Otherwise, parse start/end from filename and build a weekly sequence
  m <- str_match(basename(f),
                 ".*_(\\d{4}-\\d{2}-\\d{2})_to_(\\d{4}-\\d{2}-\\d{2})\\.nc$")
  if (any(is.na(m))) stop("Cannot parse start/end dates from filename: ", f)
  d0 <- as.Date(m[2]); d1 <- as.Date(m[3])
  ds <- seq(d0, d1, by = "7 days")
  # Align length to number of layers if needed
  if (length(ds) < nlyr(r)) {
    warning("Date length < layer count; truncating to date length: ", basename(f))
    return(ds)
  }
  ds[seq_len(nlyr(r))]
}

# ---- Helper: list model files (sorted by start date) ----
.list_model_files <- function(in_dir, var, model, scen){
  patt  <- paste0("^", var, "_weekAvg_output_", model, "_", scen,
                  "_\\d{4}-\\d{2}-\\d{2}_to_\\d{4}-\\d{2}-\\d{2}\\.nc$")
  flist <- list.files(in_dir, full.names = TRUE, pattern = patt)
  if (length(flist) == 0) stop("No files found: model=", model, ", scen=", scen)
  starts <- str_match(basename(flist), ".*_(\\d{4}-\\d{2}-\\d{2})_to_")[,2]
  flist[order(ymd(starts))]
}

# ---- Helper: read & concatenate all chunks for a model into a weekly raster time-series ----
.read_model_weekly <- function(in_dir, var, model, scen) {
  patt  <- paste0("^", var, "_weekAvg_output_", model, "_", scen,
                  "_\\d{4}-\\d{2}-\\d{2}_to_\\d{4}-\\d{2}-\\d{2}\\.nc$")
  flist <- list.files(in_dir, full.names = TRUE, pattern = patt)
  if (length(flist) == 0) stop("No files found for model=", model, ", scen=", scen)
  starts <- str_match(basename(flist), ".*_(\\d{4}-\\d{2}-\\d{2})_to_")[,2]
  flist <- flist[order(ymd(starts))]
  
  rs_list <- list(); ds_list <- list()
  for (f in flist) {
    r  <- rast(f)
    # try to obtain dates from the raster time; otherwise parse from filename
    tt <- time(r)
    if (!is.null(tt) && length(tt) >= nlyr(r)) {
      ds <- as.Date(tt)[seq_len(nlyr(r))]
    } else {
      m <- str_match(basename(f),
                     ".*_(\\d{4}-\\d{2}-\\d{2})_to_(\\d{4}-\\d{2}-\\d{2})\\.nc$")
      if (any(is.na(m))) stop("Cannot parse dates from filename: ", f)
      d0 <- as.Date(m[2]); d1 <- as.Date(m[3])
      ds <- seq(d0, d1, by = "7 days")
      if (length(ds) < nlyr(r)) {
        warning("Date length < layer count; truncating for: ", basename(f))
        ds <- ds[seq_len(length(ds))]
        r  <- r[[seq_len(length(ds))]]
      } else if (length(ds) > nlyr(r)) {
        ds <- ds[seq_len(nlyr(r))]
      }
    }
    rs_list[[length(rs_list) + 1]] <- r
    ds_list[[length(ds_list) + 1]] <- ds
  }
  r_all <- do.call(c, rs_list)
  d_all <- do.call(c, ds_list)
  time(r_all) <- d_all
  list(r = r_all, dates = d_all)
}

# =====================================================================
# 1) Single-model annual means (per-year, per-pixel mean via app())
# ---------------------------------------------------------------------
# Returns: SpatRaster with one layer per year (layer names: Y{year})
# Optionally writes a multi-layer GeoTIFF (one layer per year).
# =====================================================================

compute_model_annual_means <- function(in_dir   = "waterTemp",
                                       out_dir  = NULL,
                                       var      = "waterTemp",        # "waterTemp" or "discharge"
                                       model,
                                       scen     = "hist",             # "hist" or "rcp2p6"
                                       threshold = NULL,             # unified threshold (NULL -> use defaults)
                                       cores    = 32,
                                       overwrite = TRUE,
                                       return_rasters = FALSE) {
  # Input validation
  stopifnot(!missing(model))
  if (!is.null(out_dir)) dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Set default threshold if NULL depending on var
  if (is.null(threshold)) {
    if (var == "waterTemp") {
      threshold <- 350   # K
    } else if (var == "discharge") {
      threshold <- 10    # m3/s
    } else {
      stop("Unknown var. Supported: 'waterTemp' or 'discharge'.")
    }
  }
  
  # Read weekly stacks and dates for this model (assumes .read_model_weekly exists)
  mw <- .read_model_weekly(in_dir, var, model, scen)
  yrs <- sort(unique(year(mw$dates)))
  
  # Prepare results containers
  written_files <- character(0)
  if (return_rasters) {
    raster_layers <- vector("list", length(yrs))
    rl_idx <- 1
  } else {
    raster_layers <- NULL
  }
  
  # Loop years: compute annual mean and write immediately
  for (yy in yrs) {
    # select layer indices for this year
    idx <- which(year(mw$dates) == yy)
    if (length(idx) == 0) next
    
    # Construct per-year output filename
    fname_year <- paste0(var, "_annual_mean_", model, "_", scen, "_", yy, ".tif")
    out_path <- if (is.null(out_dir)) fname_year else file.path(out_dir, fname_year)
    
    # Skip if exists and overwrite == FALSE
    if (file.exists(out_path) && !overwrite) {
      message("[skip] Year ", yy, " exists and overwrite = FALSE: ", out_path)
      written_files <- c(written_files, normalizePath(out_path))
      if (return_rasters) {
        # optionally, load existing file into memory for return
        raster_layers[[rl_idx]] <- rast(out_path)
        rl_idx <- rl_idx + 1
      }
      next
    }
    
    # Subset weekly layers for the year (SpatRaster with layers = weeks)
    ryy <- mw$r[[idx]]
    
    # -----------------------------
    # Mask unrealistic weekly values BEFORE annual aggregation
    # - For waterTemp: set weekly values > threshold -> NA
    # - For discharge: set weekly values < threshold -> NA
    # - This removes unrealistic weekly outliers prior to computing the yearly mean
    # -----------------------------
    if (var == "waterTemp") {
      # English: set weekly temperature values greater than threshold to NA
      ryy[ryy > threshold] <- NA
    } else if (var == "discharge") {
      # English: set weekly discharge values smaller than threshold to NA
      ryy[ryy < threshold] <- NA
    } else {
      stop("Unsupported var: ", var)
    }
    
    # Pixel-wise annual mean using app()
    # English comment: compute pixel-wise mean across the selected (masked) weekly layers
    ann <- app(ryy, fun = mean, na.rm = TRUE, cores = cores)
    
    # Set a meaningful layer name
    names(ann) <- paste0("Y", yy)
    
    # Write the single-year raster to disk immediately
    # English comment: write out the computed single-year raster to disk
    writeRaster(ann, filename = out_path, overwrite = overwrite)
    
    # Record the file path
    written_files <- c(written_files, normalizePath(out_path))
    
    # If requested, keep the raster layer in memory for returning later
    if (return_rasters) {
      raster_layers[[rl_idx]] <- ann
      rl_idx <- rl_idx + 1
    }
    
    message("Wrote year ", yy, " -> ", out_path)
    # Free memory for ann (R will garbage collect as needed)
    rm(ann); gc()
  }
  
  # Optionally combine kept raster layers into a multi-layer SpatRaster
  if (return_rasters && length(raster_layers) > 0) {
    # Some entries may be NULL if files were skipped/exists; keep non-null
    raster_layers <- raster_layers[!vapply(raster_layers, is.null, is.null)]
    if (length(raster_layers) > 0) {
      combined <- rast(raster_layers)
      names(combined) <- paste0("Y", yrs[seq_len(nlyr(combined))])
      return(list(files = written_files, raster = combined))
    } else {
      return(list(files = written_files, raster = NULL))
    }
  }
  
  # Default return: vector of file paths written (or existing)
  return(invisible(written_files))
}


# =====================================================================
# 2) Multi-model ensemble annual means (common years across models)
# ---------------------------------------------------------------------
# For each common year, averages the single-model annual layers using app()
# Returns: SpatRaster with one layer per common year (layer names: Y{year})
# Optionally writes a multi-layer GeoTIFF.
# =====================================================================
compute_ensemble_annual_means <- function(in_dir    = "waterTemp_annual_mean",
                                          out_dir   = "waterTemp_annual_mean_ensemble",
                                          var       = "waterTemp",
                                          models    = c("gfdl","hadgem","ipsl","miroc","noresm"),
                                          scen      = "hist",
                                          cores     = 32,
                                          overwrite = TRUE) {
  # validate
  if (!dir.exists(in_dir)) stop("in_dir does not exist: ", in_dir)
  if (!is.null(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("Looking for files in: ", in_dir)
  # collect files per model (pattern strict)
  per_model_files <- list()
  for (m in models) {
    patt <- paste0("^", var, "_annual_mean_", m, "_", scen, "_[0-9]{4}\\.tif$")
    f <- list.files(in_dir, pattern = patt, full.names = TRUE, ignore.case = FALSE)
    if (length(f) == 0) {
      stop("Missing files for required model '", m, "'. Expected pattern: ", patt,
           "\nFunction requires all specified models to be present.")
    }
    per_model_files[[m]] <- sort(f)
  }
  
  # map filepaths to year names
  per_model_year_map <- lapply(per_model_files, function(fp) {
    yrs <- sub(".*_([0-9]{4})\\.tif$", "\\1", basename(fp))
    names(fp) <- yrs
    return(fp)
  })
  
  # compute intersection of years across all specified models
  years_list <- lapply(per_model_year_map, names)
  common_years <- Reduce(intersect, years_list)
  if (length(common_years) == 0) stop("No common years across all specified models.")
  common_years <- sort(as.integer(common_years))
  
  message("Models: ", paste(models, collapse = ", "))
  message("Common years (n=", length(common_years), "): ", paste(common_years, collapse = ", "))
  
  # models tag for filename
  models_tag <- paste(models, collapse = "-")
  
  # compute per-year ensemble and write single-layer tif
  mm_layers_mem <- vector("list", length(common_years))
  idx <- 1
  for (yy in common_years) {
    yy_chr <- as.character(yy)
    # get ordered file for each model
    files_for_year <- vapply(models, function(m) {
      mm <- per_model_year_map[[m]]
      if (yy_chr %in% names(mm)) mm[[yy_chr]] else NA_character_
    }, character(1), USE.NAMES = FALSE)
    
    # sanity: should be all present because we used intersection
    if (any(is.na(files_for_year))) {
      warning("Year ", yy, " skipped because not all models present (unexpected).")
      next
    }
    
    # read & stack rasters
    ras_list <- lapply(files_for_year, rast)
    X <- rast(ras_list)      # layers = models
    
    # compute pixel-wise mean across models
    mm_y <- app(X, fun = mean, na.rm = TRUE, cores = cores)
    names(mm_y) <- paste0("Y", yy)
    
    # write per-year file with required naming convention
    if (!is.null(out_dir)) {
      outf <- file.path(out_dir, sprintf("%s_annual_mean_ensemble_%s_%d.tif",
                                         var, scen, yy))
      writeRaster(mm_y, outf, overwrite = overwrite)
      message("Wrote: ", outf, "  (models=", length(files_for_year), ")")
    }
    
    mm_layers_mem[[idx]] <- mm_y
    idx <- idx + 1
  }
  
  # combine produced layers (drop trailing NULLs if any)
  mm_layers_mem <- mm_layers_mem[!vapply(mm_layers_mem, is.null, logical(1))]
  if (length(mm_layers_mem) == 0) return(invisible(NULL))
  
  r_mm <- rast(mm_layers_mem)
  return(r_mm)
}


# =====================================================================
# Example usage
# =====================================================================
# Save each year's single-layer tif under "waterTemp_annual_mean", do not return combined raster
# =====================================================================
compute_model_annual_means(in_dir = "input/raw/waterTemp",
                          out_dir = "input/processed/waterTemp_annual_mean",
                          var = "waterTemp",
                          model = "gfdl", #c("gfdl", "hadgem", "ipsl", "miroc", "noresm")
                          scen = "hist", # "hist" or "rcp2p6"
                          threshold = 350,    # will default to 350 for waterTemp
                          cores = 40,
                          overwrite = FALSE,
                          return_rasters = FALSE)

compute_model_annual_means(in_dir = "input/raw/discharge",
                          out_dir = "input/processed/discharge_annual_mean",
                          var = "discharge",
                          model = "gfdl", #c("gfdl", "hadgem", "ipsl", "miroc", "noresm")
                          scen = "rcp2p6", # "hist" or "rcp2p6"
                          threshold = 10,    # will default to 10 for discharge
                          cores = 20,
                          overwrite = FALSE,
                          return_rasters = FALSE)
# =====================================================================
# Save each year's single-layer tif under "waterTemp_annual_mean_ensemble", do not return combined raster
# =====================================================================
compute_ensemble_annual_means(
                          in_dir  = "input/processed/waterTemp_annual_mean",
                          out_dir = "input/processed/waterTemp_annual_mean_ensemble",
                          var     = "waterTemp",
                          models  = c("gfdl", "hadgem", "ipsl", "miroc", "noresm"),
                          scen    = "rcp2p6",         # "hist" or "rcp2p6"
                          cores   = 8,
                          overwrite = TRUE
)

compute_ensemble_annual_means(
                          in_dir  = "input/processed/discharge_annual_mean",
                          out_dir = "input/processed/discharge_annual_mean_ensemble",
                          var     = "discharge",
                          models  = c("gfdl", "hadgem", "ipsl", "miroc", "noresm"),
                          scen    = "rcp2p6",         # "hist" or "rcp2p6"
                          cores   = 8,
                          overwrite = TRUE
)


################################################################################
wt <- list.files("input/processed/waterTemp_annual_mean_ensemble",full.names = T)
wt <- rast(wt[1:50])      # layers = models

# compute pixel-wise mean across models
wt <- app(wt, fun = mean, na.rm = TRUE, cores = 20)
names(wt) <- "Y1758-1975"

# write 1758-1975 file with required naming convention
writeRaster(wt, "input/processed/waterTemp_annual_mean_ensemble/waterTemp_annual_mean_ensemble_1758_1975.tif")

#-------------------------------------------------------------------------------
discharge <- list.files("input/processed/discharge_annual_mean_ensemble",full.names = T)
discharge <- rast(discharge[1:50])      # layers = models

# compute pixel-wise mean across models
discharge <- app(discharge, fun = mean, na.rm = TRUE, cores = 20)
names(discharge) <- "Y1758-1975"

# write 1758-1975 file with required naming convention
writeRaster(discharge, "input/processeddischarge_annual_mean_ensemble/discharge_annual_mean_ensemble_1758_1975.tif")

