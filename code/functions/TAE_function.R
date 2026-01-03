#' Calculate the Taxonomic Effort Index (TAE) from Eschmeyer's Catalog of Fishes
#'
#' These functions provide a workflow to retrieve reference metadata from the
#' [Eschmeyer's Catalog of Fishes](https://researcharchive.calacademy.org/research/ichthyology/catalog/fishcatmain.asp),
#' standardize author information, and compute the **Taxonomic Effort Index (TAE)** following
#' a publication- and author-weighted formulation. The TAE index integrates
#' temporal coverage and authorship information to quantify cumulative taxonomic effort.
#'
#' The workflow consists of four functions:
#' \itemize{
#'   \item \code{cas_ref()}: Retrieve a single reference by ID.
#'   \item \code{get_cas_ref()}: Retrieve multiple references.
#'   \item \code{clean_refs()}: Standardize author name strings.
#'   \item \code{calculate_TAE()}: Compute the final TAE index.
#' }
#'
#' @name TAE-workflow
#' @author Dr. Liuyong Ding <ly_ding@126.com>
NULL


#' Retrieve a Single Reference Record
#'
#' Query Eschmeyer's Catalog of Fishes by reference ID and parse key metadata,
#' including authorship, publication year, and full citation text.
#'
#' @inheritParams get_cas_ref
#'
#' @return A [tibble][tibble::tibble-package] containing:
#' \itemize{
#'   \item \code{ref_id}: Reference ID.
#'   \item \code{ref_authorship}: Authorship string.
#'   \item \code{ref_year}: Publication year.
#'   \item \code{reference}: Full reference text.
#'   \item \code{ref_url}: URL to the original reference.
#' }
#'
#' @importFrom httr GET status_code
#' @importFrom xml2 read_html xml_find_all xml_text xml_find_first xml_attr
#' @importFrom dplyr tibble select
#' @importFrom stringr str_extract
#' @importFrom cli cli_alert_info
#'
#' @seealso [Eschmeyer's Catalog of Fishes](https://researcharchive.calacademy.org/research/ichthyology/catalog/fishcatmain.asp)
#' @family TAE-workflow
#'
#' @export
cas_ref <- function(query, quiet) {
  if (!quiet) cli::cli_alert_info("Searching for ref. {query} on CAS")
  baseurl <- "https://researcharchive.calacademy.org/research/ichthyology/catalog/getref.asp?"
  stopifnot(is.numeric(query))
  
  url <- paste0(baseurl, "id=", query)
  res <- httr::GET(url)
  html <- res %>% xml2::read_html()
  
  if (httr::status_code(res) == 200) {
    table <- dplyr::tibble(
      ref_id = query,
      all_data = xml2::xml_find_all(html, "//dt[@class='result']") %>%
        xml2::xml_text(),
      ref_authorship = stringr::str_extract(all_data, ".+(?=\\s\\s[0-9]{4}(\\s|\\-))"),
      ref_year = stringr::str_extract(all_data, "[0-9]{4}"),
      reference = xml2::xml_find_all(html, "//dd[@class='result']") %>%
        xml2::xml_text(),
      ref_url = xml2::xml_find_all(html, "//dt[@class='result']") %>%
        xml2::xml_find_first('a') %>%
        xml2::xml_attr('href')
    ) %>%
      dplyr::select(-all_data)
    
    return(table)
  } else {
    cat("Error request - the query parameter is not valid")
    browseURL(baseurl)
  }
}


#' Retrieve Multiple Reference Records
#'
#' Fetch multiple reference records by their IDs, with optional CLI progress display.
#'
#' @param query Numeric. One or more reference IDs to query.
#' @param quiet Logical. If \code{TRUE}, suppress CLI messages. Default is \code{FALSE}.
#'
#' @return A [tibble][tibble::tibble-package] with one row per reference.
#'
#' @importFrom purrr map list_rbind
#' @family TAE-workflow
#'
#' @examples
#' \dontrun{
#' # Single reference
#' r1 <- get_cas_ref(4883)
#'
#' # Multiple references
#' r2 <- get_cas_ref(c(4883, 35823))
#' }
#'
#' @export
get_cas_ref <- function(query, quiet = FALSE) {
  
  if (length(query) > 1) {
    res <- purrr::map(
      query,
      .f = function(x) cas_ref(x, quiet),
      .progress = list(
        format = "Query {cli::pb_current+1}/{cli::pb_total} | {cli::pb_bar} { {paste0(round({{cli::pb_current}/{cli::pb_total}}*100),'%')} } [{cli::pb_elapsed}]",
        show_after = 0,
        clear = TRUE
      )
    ) %>%
      purrr::list_rbind()
  } else {
    res <- cas_ref(query, quiet)
  }
  
  return(res)
  cat("Error request - the query parameter is not valid")
  browseURL(baseurl)
}


#' Clean and Standardize Authorship Strings
#'
#' Standardize author name strings by trimming whitespace, swapping first/last names,
#' and normalizing separators. This ensures consistent author counts for TAE calculation.
#'
#' @param df A data frame containing a column with authorship strings.
#' @param col_name Unquoted column name containing authorship strings.
#'
#' @return The input data frame with standardized author strings in the target column.
#'
#' @importFrom dplyr mutate enquo
#' @family TAE-workflow
#'
#' @examples
#' df <- tibble::tibble(ref_authorship = c("Smith,J., and Brown,K.", "Li,Y."))
#' clean_refs(df, ref_authorship)
#'
#' @export
clean_refs <- function(df, col_name) {
  
  process_string <- function(input_string) {
    input_string <- trimws(input_string)
    
    if (grepl("\\.\\,", input_string)) {
      split_by_dot_comma <- strsplit(input_string, "\\.\\,")[[1]]
      first_author <- trimws(split_by_dot_comma[1])
      first_author <- sub("^(\\S+),(.*)$", "\\2. \\1", first_author)
      remaining_authors <- split_by_dot_comma[2]
      remaining_authors <- gsub(" and", ",", remaining_authors)
      remaining_authors <- trimws(remaining_authors)
      final_authors <- paste(first_author, remaining_authors, sep = ", ") %>% trimws()
    } else {
      final_authors <- sub("^(\\S+),\\s*(.*)$", "\\2 \\1", input_string)
    }
    
    return(final_authors)
  }
  
  col_name_quo <- enquo(col_name)
  df <- df %>%
    mutate(!!col_name_quo := sapply(!!col_name_quo, process_string))
  
  return(df)
}


#' Compute the Taxonomic Effort Index (TAE)
#'
#' Calculate the TAE value from reference metadata, following a formulation that
#' combines temporal coverage (\eqn{\sqrt{D}}), author numbers (\eqn{\log(1 + A_i)}),
#' and temporal decay (\eqn{\exp(-\lambda \times \mathrm{age})}).
#'
#' The final formula is:
#' \deqn{\mathrm{TAE} = \sqrt{D} \times \frac{1}{N} \sum_{i=1}^N \log(1 + A_i) \exp(-\lambda (t_{\max} - t_i))}
#' where:
#' \itemize{
#'   \item \eqn{D = t_{\max} - t_{\min} + 1}.
#'   \item \eqn{\lambda = \ln(2)/20}, corresponding to a 20-year half-life.
#'   \item \eqn{A_i} is the number of unique authors for reference \eqn{i}.
#' }
#'
#' @param df A data frame containing \code{ref_year} (numeric) and \code{ref_authorship}.
#' @param lambda Numeric. Temporal decay parameter. Default is \eqn{\ln(2)/20}.
#' @param author_sep Character. Separator used to split authors. Default is \code{","}.
#'
#' @return A single numeric value representing the TAE index.
#'
#' @importFrom stats na.omit
#' @family TAE-workflow
#'
#' @examples
#' \dontrun{
#' refs <- get_cas_ref(c(4883, 35823))
#' tae_value <- calculate_TAE(refs)
#' }
#'
#' @export
calculate_TAE <- function(df,
                          lambda = log(2) / 20,
                          author_sep = ",") {
  
  if (is.null(df) || nrow(df) == 0) return(NA_real_)
  
  df2 <- df[!is.na(df$ref_year), , drop = FALSE]
  if (nrow(df2) == 0) return(NA_real_)
  df2$ref_year <- as.integer(df2$ref_year)
  
  if (exists("clean_refs")) {
    df2 <- clean_refs(df2, "ref_authorship")
  }
  
  split_auth <- strsplit(df2$ref_authorship %||% "", paste0("\\s*", author_sep, "\\s*"), perl = TRUE)
  Ai <- vapply(split_auth, function(x) {
    x <- trimws(x)
    x <- x[nzchar(x)]
    length(unique(x))
  }, integer(1))
  
  Ai[is.na(Ai) | Ai < 0] <- 0L
  
  t_max <- max(df2$ref_year, na.rm = TRUE)
  t_min <- min(df2$ref_year, na.rm = TRUE)
  D <- (t_max - t_min + 1L)
  
  N <- length(Ai)
  if (N == 0) return(NA_real_)
  
  age <- pmax(0, t_max - df2$ref_year)
  w_i <- log1p(Ai) * exp(-lambda * age)
  
  TAE <- sqrt(D) * mean(w_i)
  return(TAE)
}


# small helper for `%||%`
`%||%` <- function(x, y) if (is.null(x)) y else x