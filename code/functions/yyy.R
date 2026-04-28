library(rgbif)
library(dplyr)
library(purrr)
library(tibble)

gbif_preserved_specimens <- function(
    species_name,
    include_synonyms = FALSE,     
    basisOfRecord = "PRESERVED_SPECIMEN",
    verbose = FALSE
) {
  # 1) backbone：
  bb <- tryCatch(name_backbone(name = species_name), error = function(e) NULL)
  if (is.null(bb) || is.null(bb$usageKey)) {
    return(list(
      species_input = species_name,
      accepted_name = NA_character_,
      usageKey = NA_integer_,
      n_preserved_accepted = NA_integer_,
      n_preserved_sum_keys = NA_integer_,
      key_counts = NULL,
      synonyms = NULL,
      status = "backbone_failed"
    ))
  }
  
  acc_key  <- suppressWarnings(as.integer(bb$usageKey))
  acc_name <- bb$scientificName %||% bb$canonicalName %||% species_name
  
  # 2) accepted taxonKey
  n_acc <- tryCatch(
    occ_count(taxonKey = acc_key, basisOfRecord = basisOfRecord),
    error = function(e) NA_integer_
  )
  
  out <- list(
    species_input = species_name,
    accepted_name = acc_name,
    usageKey = acc_key,
    n_preserved_accepted = n_acc,
    n_preserved_sum_keys = NA_integer_,
    key_counts = NULL,
    synonyms = NULL,
    status = "ok"
  )
  
  if (!include_synonyms) return(out)
  
  # 3) synonyms
  nu <- tryCatch(name_usage(key = acc_key, data = "synonyms"), error = function(e) NULL)
  syn <- NULL
  if (!is.null(nu) && !is.null(nu$data) && nrow(nu$data) > 0) {
    syn <- nu$data %>%
      transmute(
        synonym_name = scientificName,
        status = taxonomicStatus,
        rank = rank,
        key = suppressWarnings(as.integer(key))
      ) %>%
      filter(!is.na(key))
  } else {
    syn <- tibble(synonym_name = character(), status = character(), rank = character(), key = integer())
  }
  
  # 4) accepted + distinct synonym keys
  count_key <- function(k) {
    tryCatch(
      occ_count(taxonKey = k, basisOfRecord = basisOfRecord),
      error = function(e) NA_integer_
    )
  }
  
  key_tbl <- tibble(
    type = c("accepted", rep("synonym", nrow(syn))),
    taxonKey = c(acc_key, syn$key)
  ) %>%
    distinct(taxonKey, .keep_all = TRUE)
  
  key_counts <- key_tbl %>%
    mutate(n = map_int(taxonKey, count_key)) %>%
    mutate(n = if_else(is.na(n), 0L, n))
  
  out$synonyms <- syn
  out$key_counts <- key_counts
  out$n_preserved_sum_keys <- sum(key_counts$n, na.rm = TRUE)
  
  if (verbose) {
    message("Species: ", species_name,
            " | acceptedKey: ", acc_key,
            " | N(accepted): ", n_acc,
            " | N(sum keys): ", out$n_preserved_sum_keys,
            " | synonyms: ", nrow(syn))
  }
  
  out
}

`%||%` <- function(x, y) if (!is.null(x) && length(x) > 0 && !is.na(x)) x else y


gbif_preserved_batch <- function(
    species_vec,
    include_synonyms = FALSE,
    basisOfRecord = "PRESERVED_SPECIMEN"
) {
  res <- purrr::map(species_vec, gbif_preserved_specimens,
                    include_synonyms = include_synonyms,
                    basisOfRecord = basisOfRecord)
  
  main <- purrr::map_dfr(res, \(x) tibble(
    species_input = x$species_input,
    accepted_name = x$accepted_name,
    usageKey = x$usageKey,
    n_preserved_accepted = x$n_preserved_accepted,
    n_preserved_sum_keys = x$n_preserved_sum_keys,
    status = x$status
  ))
  
  if (!include_synonyms) return(list(main = main))
  
  
  syn_tbl <- purrr::map_dfr(res, \(x) tibble(
    species_input = x$species_input,
    synonyms = list(x$synonyms)
  ))
  key_tbl <- purrr::map_dfr(res, \(x) tibble(
    species_input = x$species_input,
    key_counts = list(x$key_counts)
  ))
  
  list(main = main, synonyms = syn_tbl, key_counts = key_tbl)
}
