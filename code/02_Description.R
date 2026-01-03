library(dplyr)
library(purrr)
library(tidyr)
library(stringr)


df <- openxlsx::read.xlsx("input/raw/cas_freshwater_v1.xlsx")

df_cleaned <- df %>%
  mutate(basin = str_split(basin, ";")) %>%        
  mutate(basin = lapply(basin, unique)) %>%        
  mutate(basin = lapply(basin, sort))   %>%
  unnest(basin) %>%
  dplyr::select(valid_name,year_description,basin)

biogeographic_list <- read.csv("input/raw/biogeographic_list.csv")
df_cleaned <- df_cleaned %>% left_join(biogeographic_list, by = "basin")
df_cleaned <- df_cleaned[,c("basin_id","basin","valid_name","year_description")]

head(df_cleaned)

out_dir <- "input/processed/basin_sp"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

df_cleaned %>%
  group_by(basin_id) %>%
  group_walk(~ {
    b <- .y$basin_id[[1]]
    out <- dplyr::bind_cols(.y, .x)   
    file_path <- file.path(out_dir, paste0(b, ".csv"))
    write.csv(out, file_path, row.names = FALSE)
    message("Saved: ", file_path)
  })
