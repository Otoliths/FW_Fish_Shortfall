# We extracted the number of preserved specimens per species deposited in biological collections worldwide using the
# function occ_count in the rgbif R package. For that, we created search queries containing valid species name plus their
# unique synonyms (i.e. invalid names that can be traced back to a single valid name) and set the basisOfRecord argument to the
# preserved specimen

source("code/functions/yyy.r")
df <- openxlsx::read.xlsx("input/raw/cas_freshwater_v1.xlsx")

batch <- gbif_preserved_batch(df$valid_name,include_synonyms = TRUE)

batch_df <- batch[["main"]]
specimen_df <- batch_df %>%
  filter(status == "ok") %>%
  select(
    valid_name = species_input,
    preserved_specimen = n_preserved_sum_keys
  ) %>%
  mutate(
    preserved_specimen = ifelse(is.na(preserved_specimen), 0, preserved_specimen)
  )

saveRDS(specimen_df,"input/data_prep/preserved_specimens.rds")
