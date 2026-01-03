library(dplyr)
data <- openxlsx::read.xlsx("input/raw/cas_freshwater_v1.xlsx")
biogeographic_list <- read.csv("input/raw/biogeographic_list.csv")
linnaean <- read.csv("output/tables/basin_linnaean_shortfall.csv")
sp <- data %>%
  mutate(basin = str_split(basin, ";")) %>%        # 1. 拆分
  mutate(basin = lapply(basin, unique)) %>%        # 2. 去重复
  mutate(basin = lapply(basin, sort)) %>%
  unnest(basin) %>%
  dplyr::select(valid_name,year_description,basin)


sp_recent_2000 <- sp %>%
  filter(year_description >= 2000, year_description <= 2024) %>%
  left_join(biogeographic_list[,1:2],by = "basin")



lambda_df <- sp_recent_2000 %>%
  group_by(basin_id) %>%
  summarise(
    n_desc = n(),                     
    lambda = n_desc / 25              
  )
range(lambda_df$lambda)

lambda_full <- biogeographic_list %>%
  left_join(lambda_df, by = "basin_id") %>%
  mutate(
    n_desc = replace_na(n_desc, 0),
    lambda = replace_na(lambda, 0)
  )

lambda_q <- quantile(lambda_full$lambda[lambda_full$lambda > 0], 0.05)

#λq​=0.04 species per year.
#In the slowest 5% of basins where at least one new species has been described since 2000, 
#the long-term rate corresponds to describing one species every ~25 years.

df <- lambda_full %>%
  left_join(linnaean[, c("basin_id", "SRdesc")], by = "basin_id") %>%
  rename(Ux = SRdesc)

df <- df %>%
  mutate(
    lambda_adj = pmax(lambda, lambda_q),  
    cost_years = Ux / lambda_adj          
  )

range(df$cost_years,na.rm = T)
hist(log(df$cost_years+1))
saveRDS(df,"input/processed/basin_cost_taxonomy.rds")


#################################################################################
library(dplyr)
data <- readRDS("input/processed/country_sp.rds")
country_list <- read.csv("input/raw/country_list.csv")
linnaean <- read.csv("output/tables/country_linnaean_shortfall.csv")
sp <- data %>%
  mutate(iso3 = str_split(iso3, ";")) %>%        
  mutate(iso3 = lapply(iso3, unique)) %>%        
  mutate(iso3 = lapply(iso3, sort)) %>%
  unnest(iso3) %>%
  dplyr::select(valid_name,year_description,iso3)


sp_recent_2000 <- sp %>%
  filter(year_description >= 2000, year_description <= 2024) %>%
  left_join(country_list[,1:2],by = "iso3")



lambda_df <- sp_recent_2000 %>%
  group_by(iso3) %>%
  summarise(
    n_desc = n(),                     
    lambda = n_desc / 25              
  )
range(lambda_df$lambda)

lambda_full <- country_list %>%
  left_join(lambda_df, by = "iso3") %>%
  mutate(
    n_desc = replace_na(n_desc, 0),
    lambda = replace_na(lambda, 0)
  )

lambda_q <- quantile(lambda_full$lambda[lambda_full$lambda > 0], 0.05)

#λq​=0.04 species per year.
#In the slowest 5% of basins where at least one new species has been described since 2000, 
#the long-term rate corresponds to describing one species every ~25 years.

df <- lambda_full %>%
  left_join(linnaean[, c("iso3", "SRdesc")], by = "iso3") %>%
  rename(Ux = SRdesc)

df <- df %>%
  mutate(
    lambda_adj = pmax(lambda, lambda_q),  
    cost_years = Ux / lambda_adj          
  )

range(df$cost_years,na.rm = T)
hist(log(df$cost_years+1))
saveRDS(df,"input/processed/country_cost_taxonomy.rds")
