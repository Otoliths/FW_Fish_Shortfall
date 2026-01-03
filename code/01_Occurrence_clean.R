library(data.table)
library(CoordinateCleaner)
library(dplyr)
library(purrr)
library(tidyr)
library(stringr)
library(sf)
library(stringi)
library(readxl)
library(progressr)

handlers(global = TRUE)
handlers("progress")
options(sf_use_s2 = FALSE)

source("code/functions/clean_species_occ.R")

names_only <- names(fread(
  "input/raw/occ/gbif/0060713-250920141307145/occurrence.txt",
  nrows = 0,
  sep = "\t",
  quote = ""
))
head(names_only)

cols_to_keep <- c("species", "decimalLongitude", #"lifeStage","habitat","sex",
                  "decimalLatitude", "year", "basisOfRecord")

gbif <- fread(
  "input/raw/occ/gbif/0060713-250920141307145/occurrence.txt",
  select = cols_to_keep,
  quote = "",
  encoding = "UTF-8",
  showProgress = TRUE
)
gbif$database_name <- "GBIF"

nrow(gbif) #36302470
length(unique(gbif$species)) #14387

#########
obis <- arrow::read_parquet("input/raw/occ/obis/obis.parquet")
obis <- obis[,cols_to_keep]
obis$year <- as.numeric(obis$year)
obis$database_name <- "OBIS"
obis <- obis %>%
  filter(!is.na(decimalLongitude)) %>%
  filter(!is.na(decimalLatitude))
nrow(obis) #36309270
length(unique(obis$species)) #16842

iucn <- read.csv("input/raw/occ/iucn/iucn_points_data.csv")
names(iucn)[19] <- "decimalLongitude"
names(iucn)[20] <- "decimalLatitude"
names(iucn)[3] <- "species"
names(iucn)[17] <- "basisOfRecord"
iucn <- iucn %>%
  filter(!is.na(decimalLongitude)) %>%
  filter(!is.na(decimalLatitude))
iucn <- iucn[,c("species","decimalLongitude","decimalLatitude","year","basisOfRecord")]
iucn$database_name <- "IUCN"
iucn$year <- as.numeric(iucn$year)
nrow(iucn) #113596
length(unique(iucn$species)) #1861

cas <- read.csv("input/raw/occ/cas/cas.csv") %>%
  filter(!is.na(decimalLongitude)) %>%
  filter(!is.na(decimalLatitude))
cas <- cas[,c("species","decimalLongitude","decimalLatitude","year")]
cas$database_name <- "CAS"
cas$basisOfRecord <- "Type_location"
nrow(cas) #5787
length(unique(cas$species)) #5787

occ <- rbind(gbif,obis,iucn,cas)
#arrow::write_parquet(occ, "input/processed/occ.parquet")
occ <- arrow::read_parquet("input/processed/occ.parquet")
table(occ$basisOfRecord)
################################################################################

inland <- readRDS("input/raw/basin/basin_sf_v1.rds")
df <- read_excel("input/raw/cas_freshwater_v1.xlsx") %>%
  dplyr::select(valid_name,basin) %>%
  mutate(basin = str_split(basin, ";")) %>%        
  mutate(basin = lapply(basin, unique)) %>%        
  mutate(basin = lapply(basin, sort))   %>%
  unnest(basin)
sp_list <- unique(df$valid_name)


# Process a single species
#clean_species_occ(sp = "Garra mcclellandi", data = df, polygon = inland, occ,out_dir = "input/processed/occ")

# Batch process all species (with progress bar)
progressr::with_progress({
  p <- progressor(along = sp_list)
  
  # run one species at a time, show species name as the tick message
  invisible(map(sp_list, ~{
    suppressMessages(    # keep console clean; progress bar shows status
      clean_species_occ(.x, data = df, polygon = inland, occ,out_dir = "input/processed/occ")
    )
    p(message = .x)
  }))
})


