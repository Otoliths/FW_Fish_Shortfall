get_nelson_et_al_fix <- function(ranges = "50k_50mio") {
  if (any(!ranges %in% mapme.biodiversity:::.nelson_df$range)) {
    index <- which(!ranges %in% mapme.biodiversity:::.nelson_df$range)
    ranges <- ranges[-index]
    if (length(ranges) == 0) {
      stop("No supoorted ranges have been specified.")
    }
  }
  function(x, name = "nelson_et_al", type = "raster", outdir = mapme_options()[["outdir"]], 
           verbose = mapme_options()[["verbose"]]) {
    fps <- .get_traveltime_url_fix(ranges, paste0("traveltime-", 
                                                 ranges, ".tif"), verbose = verbose)
    make_footprints(fps, filenames = fps[["filename"]], what = "raster", 
                    co = c("-co", "COMPRESS=LZW", "-ot", "UInt16", "-a_nodata", 
                           "65535"))
  }
  
}


.get_traveltime_url_fix <- function(range, filenames, verbose = TRUE) {
  urls <- unlist(lapply(range, function(x) {
    index <- mapme.biodiversity:::.nelson_df$index[mapme.biodiversity:::.nelson_df$range == x]
    paste0("https://figshare.com/ndownloader/files/", 
           index)
  }))
  bbox <- c(xmin = -180, ymin = -60, xmax = 180, ymax = 85)
  fps <- st_as_sfc(st_bbox(bbox, crs = "EPSG:4326"))
  fps <- st_as_sf(rep(fps, length(urls)))
  fps[["source"]] <- urls
  fps[["filename"]] <- filenames
  fps
}

.has_internet <- function () 
{
  if (Sys.getenv("mapme_check_connection", unset = "TRUE") == 
      "FALSE") {
    return(TRUE)
  }
  rsp <- httr2::req_perform(httr2::request("www.baidu.com"))
  has_internet <- !httr2::resp_is_error(rsp)
  if (!has_internet) {
    message("There seems to be no internet connection. Cannot download resources.")
  }
  has_internet
}


get_resources <- function (x, ...) 
{
  x <- mapme.biodiversity:::.check_portfolio(x)
  if (!.has_internet()) {
    return(invisible(x))
  }
  funs <- purrr::map(list(...), function(fun) mapme.biodiversity:::.check_resource_fun(fun))
  for (fun in funs) mapme.biodiversity:::.get_single_resource(x, fun)
  invisible(x)
}