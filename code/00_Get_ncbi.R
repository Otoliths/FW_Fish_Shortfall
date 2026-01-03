library(rentrez)
library(restez)
library(dplyr)
library(tidyverse)
library(stringr)
library(purrr)
library(cli)

fetch_genbank_data <- function(query, gene,ncbi_key = "",path) {
  if (is.null(ncbi_key)) {
    cli::cli_alert_warning("No NCBI API key provided. Please consider applying for one using the `rentrez` package.")
    cli::cli_alert_info("Refer to the `rentrez` documentation for instructions on obtaining an API key.")
  } else {
    set_entrez_key(ncbi_key)
  }
  
  cli::cli_h3("Searching for {.strong   {query}} in NCBI Nuccore...")
  term <- paste(query, "[ORGN] AND", gene, "[GENE]")
  
  # Initialize an empty list to store search results
  search_res <- list()
  
  # Initialize the progress bar
  cli::cli_progress_bar("Searching terms", total = length(term))
  
  temp <- rentrez::entrez_search(db = "nuccore", term = paste(query, "[ORGN]"), retmax = 1, api_key = ncbi_key)
  if(temp$count > 1){
    # Use a for loop for the search process with error handling and progress bar updates
    for (i in seq_along(term)) {
      tryCatch({
        search_res[[i]] <- rentrez::entrez_search(db = "nuccore", term = term[i], retmax = 20000, api_key = ncbi_key)
      }, error = function(e) {
        cli::cli_alert_warning("Failed to search for term '{term[i]}': {conditionMessage(e)}")
        NULL  # Return NULL in case of an error
      })
      
      # Update the progress bar after each search
      cli::cli_progress_update()
    }
  }
  
  # Finish the progress bar
  cli::cli_progress_done()
  
  # Remove any NULL results (failed searches) from the list
  search_res <- search_res[!sapply(search_res, is.null)]
  
  # Remove the zero counts
  search.res.nz <- search_res[which(lapply(search_res, function(x) x$count) > 0)]
  
  # Get IDs and remove duplicates
  search.ids <- unique(unlist(lapply(search.res.nz, function(x) x$ids)))
  # get IDs and remove dups
  cli::cli_alert_info("Retrieved {.strong {length(search.ids)}} IDs.")
  
  # Create directory for output
  dn <- file.path(path,"gb",query)
  if (!dir.exists(dn)) {
    dir.create(dn, recursive = TRUE)
  }
  
  # Check if the files already exist, and only fetch if they don't
  cli::cli_alert_info("Fetching records for {.strong {length(search.ids)}} IDs...")
  
  # Initialize a progress bar with total steps equal to the length of search.ids
  cli::cli_progress_bar("Fetching records", total = length(search.ids))
  
  for (i in seq_along(search.ids)) {
    # Check if the file for the current ID exists
    if (!file.exists(paste0(dn, "/gb_", search.ids[i], ".rds"))) {
      attempts <- 0
      success <- FALSE
      
      # Retry mechanism
      while (attempts < 3 && !success) {
        attempts <- attempts + 1
        
        # Update the progress bar for the current ID being fetched
        #cli::cli_progress_step("Fetching ID {search.ids[i]} (Attempt {attempts})")
        
        # Fetch the record from NCBI
        tryCatch({
          res <- rentrez::entrez_fetch(db = "nuccore", id = search.ids[i], rettype = "gb")
          
          # Save the fetched record as an RDS file
          saveRDS(res, paste0(dn, "/gb_", search.ids[i], ".rds"))
          
          success <- TRUE  # Mark success if fetching is successful
        }, error = function(e) {
          cli::cli_alert_warning("Failed to fetch ID {search.ids[i]} (Attempt {attempts}): {conditionMessage(e)}")
          Sys.sleep(1)  # Wait for a second before retrying
        })
      }
      
      # Update the progress bar without a message if fetching was successful
      if (success) {
        cli::cli_progress_update()
      }
    } else {
      #cli::cli_alert_info("Record for ID {search.ids[i]} already exists, skipping...")
    }
  }
  
  
  cli::cli_progress_done()  # Mark the progress as done
  
  
  cli::cli_alert_info("Loaded {.strong {length(list.files(dn, pattern = '.rds'))}} files.")
  
  ids <- list.files(dn, pattern = ".rds", full.names = TRUE)
  
  GenBank_res <- list()
  
  cli::cli_h3("Processing fetched records...")
  
  # Initialize a progress bar with total steps equal to the length of ids
  cli::cli_progress_bar("Processing records", total = length(ids))
  
  for (i in seq_along(ids)) {
    # Read the RDS file
    record <- readRDS(ids[i])
    
    # Extract relevant information and store it in the GenBank_res list
    GenBank_res[[i]] <- list(
      id = stringr::str_extract(ids[i], "[0-9]+"),
      accession = restez:::extract_accession(record = record),
      organism = restez:::extract_organism(record = record),
      locus = restez:::extract_locus(record = record),
      features = restez:::extract_features(record = record)
    )
    
    # Update the progress bar without a message
    cli::cli_progress_update()
  }
  cli::cli_alert_info("Reading done...")
  
  cli::cli_h3("Cleaning up results...")
  GenBank_res_clean <- do.call("rbind", GenBank_res) %>% 
    tibble::as_tibble() %>% 
    mutate(
      id = map_chr(id, ~ .x[[1]]),
      accession = map_chr(accession, ~ .x[[1]]),
      organism = map_chr(organism, ~ .x[[1]]),
      locus = as.character(locus),
      accession = str_extract(locus, 'accession\\s*=\\s*\"\\S+\"') %>% str_replace_all('accession\\s*=\\s*\"|\"', ''),
      length = str_extract(locus, 'length\\s*=\\s*\"\\S+\"') %>% str_replace_all('length\\s*=\\s*\"|\"', ''),
      mol = str_extract(locus, 'mol\\s*=\\s*\"\\S+\"') %>% str_replace_all('mol\\s*=\\s*\"|\"', ''),
      type = str_extract(locus, 'type\\s*=\\s*\"\\S+\"') %>% str_replace_all('type\\s*=\\s*\"|\"', ''),
      domain = str_extract(locus, 'domain\\s*=\\s*\"\\S+\"') %>% str_replace_all('domain\\s*=\\s*\"|\"', ''),
      date = str_extract(locus, 'date\\s*=\\s*\"\\S+\"') %>% str_replace_all('date\\s*=\\s*\"|\"', '')
    ) %>% 
    select(-locus) %>% 
    mutate(
      features_data = map(features, function(feature_list) {
        tibble(
          organism = map_chr(feature_list, ~ .x$organism[1] %||% NA_character_),
          specimen_voucher = map_chr(feature_list, ~ .x$specimen_voucher[1] %||% NA_character_),
          geo_loc_name = map_chr(feature_list, ~ .x$geo_loc_name[1] %||% NA_character_),
          lat_lon = map_chr(feature_list, ~ .x$lat_lon[1] %||% NA_character_),
          collection_date = map_chr(feature_list, ~ .x$collection_date[1] %||% NA_character_),
          collected_by = map_chr(feature_list, ~ .x$collected_by[1] %||% NA_character_),
          identified_by = map_chr(feature_list, ~ .x$identified_by[1] %||% NA_character_),
          gene = map_chr(feature_list, ~ .x$gene[1] %||% NA_character_),
          product = map_chr(feature_list, ~ .x$product[1] %||% NA_character_)
        )
      })
    ) %>% 
    unnest(features_data, names_sep = "_") %>% 
    select(-features)
  
  # Save the cleaned results
  dt <- file.path(path,"data")
  if (!dir.exists(dt)) {
    dir.create(dt, recursive = TRUE)
  }
  
  saveRDS(GenBank_res_clean, paste0(dt,"/", gsub(" ", "_", query), ".rds"))
  
  cli::cli_alert_info("Data fetching and processing completed...")
  return(GenBank_res_clean)
}

# Example function call

#By default, the NCBI limits users to making only 3 requests per second (and rentrez enforces that limit). Users who register for an “API key” are able to make up to ten requests per second. Getting one of these keys is simple, you just need to register for “my ncbi” account then click on a button in the account settings page.
ncbi_key <- "************************************"
set_entrez_key(ncbi_key)

# make a query for genbank
# Satoh, T. P., Miya, M., Mabuchi, K., & Nishida, M. (2016). Structure and variation of the mitochondrial genome of fishes. BMC genomics, 17, 1-20. 
# Rabosky, D. L., Chang, J., Cowman, P. F., Sallan, L., Friedman, M., Kaschner, K., ... & Alfaro, M. E. (2018). An inverse latitudinal gradient in speciation rate for marine fishes. Nature, 559(7714), 392-395.
# Protein-coding gene
# 58 gens (two main categories: 37mitochondrial genes and 21nuclear genes)
Protein_coding <- c("ND1","ND2","COI","COII","ATP8","ATP6","COIII","ND3","ND4L","ND4","ND5","Cytb","ND6")
tRNA <- c("Arg","Asp","Gly","His","Ile","Leu","Lys","Met","Phe","Ser","Thr","Trp","Val","Ala","Asn","Cys","Gln","Glu","Pro","Sec","Tyr","Pyl")
rRNA <- c("12S","16S","rRNA")
nuclear <- c("4c4","enc1","ficd","glyt","hoxc6a","kiaa1239","myh6","panx2","plagl2","ptr","rag1","rag2","rhodopsin",
             "ripk4","sh3px3","sidkey","sreb2","svep1","tbr1","vcpip","zic1")
gene.syns <- c(Protein_coding,tRNA,rRNA,nuclear,
               "CO1","COX1","COXI","cytochrome",
               "subunit","COB","CYB",
               "mitochondrial","mitochondrion")
# entrez_db_searchable("nucleotide")
# result <- fetch_genbank_data(query = "Schizopyge nukiangensis", gene = gene.syns, ncbi_key, path = "ncbi")

# print(result)

################################################################################

data <- openxlsx::read.xlsx("input/raw/cas_freshwater_v1.xlsx")

for (i in data$valid_name) {
  if(!dir.exists(file.path("input/raw/ncbi/gb",i))){
    tryCatch({
      fetch_genbank_data(query = i, gene = gene.syns, ncbi_key = ncbi_key, path = "input/raw/ncbi/gb")
    }, error = function(e) {
      message("An error occurred: ", e$message)
    })
  }
}