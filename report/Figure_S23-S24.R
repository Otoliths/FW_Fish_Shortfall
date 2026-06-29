# Supplementary Figure S23-S24
# Principal component analysis (PCA) of Sustainable Development Goal (SDG)

#https://dashboards.sdgindex.org/downloads
#Localizing the SDG Index with Machine Learning and Satellite Imagery
#library(data.table) #fread
rm(list = ls())
library(dplyr)
library(sf)
library(readxl)
library(stringr)
options(sf_use_s2 = FALSE)   # Disable spherical geometry to avoid topology issues
sdg <- read.csv("https://hub.arcgis.com/api/v3/datasets/eae2def2bf3c4c01b9ff358fcd1a5a13_0/downloads/data?format=csv&spatialRefId=3857&where=1%3D1")
names(sdg)
merit <- c("iso3","Name","Overall_Score","Overall_Rank","Region",
           paste0("Goal_",seq(1:17),"_Score"))
sdg <- sdg[,merit]
group <- read.csv("input/raw/sdg/Metadata_Country_API_NY.GDP.PCAP.CD_DS2_en_csv_v2_134819.csv")
sdg <- sdg %>% left_join(group[,c(1,3)], by = c("iso3"="Country.Code"))
# write.csv(sdg, "input/processed//sdg_2024.csv",row.names = F)
country_sdg <- sdg %>% select(iso3,starts_with("Goal_"))

names(country_sdg) <- gsub("Goal","SDG",names(country_sdg))
names(country_sdg) <- gsub("_Score","",names(country_sdg))
country_list <- read.csv("input/raw/country_list.csv")
country_sdg <- country_sdg %>% left_join(country_list[,c(1,4,5),],by = "iso3") 


sdg_cols <- paste0("SDG_", 1:17)

sdg_final <- country_sdg %>%
  mutate(
    income_group = if_else(is.na(income_group), "Unknown", income_group),
    continent    = if_else(is.na(continent), "Unknown", continent)
  ) %>%
  group_by(income_group, continent) %>%
  mutate(across(all_of(sdg_cols),
                ~ if_else(is.na(.x), median(.x, na.rm = TRUE), .x))) %>%
  ungroup() %>%
  group_by(income_group) %>%
  mutate(across(all_of(sdg_cols),
                ~ if_else(is.na(.x), median(.x, na.rm = TRUE), .x))) %>%
  ungroup() %>%
  group_by(continent) %>%
  mutate(across(all_of(sdg_cols),
                ~ if_else(is.na(.x), median(.x, na.rm = TRUE), .x))) %>%
  ungroup() %>%
  mutate(across(all_of(sdg_cols),
                ~ if_else(is.na(.x), median(.x, na.rm = TRUE), .x))) %>%
  mutate(
    income_group = na_if(income_group, "Unknown"),
    continent    = na_if(continent, "Unknown")
  )


country_iso <- readxl::read_excel("input/raw/country_iso.xls")
basin <- readRDS("input/raw/basin/basin_sf_v1.rds") %>%
  st_drop_geometry() %>%
  mutate(country = str_split(country, ";")) %>%     # split multiple basins
  mutate(country = lapply(country, unique)) %>%     # remove duplicate basin names
  mutate(country = lapply(country, sort)) %>%       # sort alphabetically for consistency
  unnest(country)  %>% 
  left_join(country_iso[,1:2], by = "country")


basin_sdg <- basin %>% select(basin_id,iso3) %>%
  left_join(sdg_final[,1:18],by = "iso3") %>%
  group_by(basin_id) %>%
  summarise(across(starts_with("SDG_"), ~ mean(. , na.rm = TRUE), .names = "mean_{.col}")) %>%
  na.omit()

names(basin_sdg) <- gsub("mean_","",names(basin_sdg))
#saveRDS(basin_sdg,"input/processed/sdg_basin.rds")


################################################################################
basin_sdg <- readRDS("input/processed/sdg_basin.rds")
biogeographic_list <- read.csv("input/raw/biogeographic_list.csv")
basin_sdg <- basin_sdg %>% left_join(biogeographic_list[,c(1,3),],by = "basin_id")


sdg.pca <- prcomp(basin_sdg[,2:18], center = TRUE,scale. = TRUE) 
sdg.pca$PC1 <- sdg.pca$x[, 1] 
sdg.pca$PC2 <- sdg.pca$x[, 2] 
sdg.pca$PC3 <- sdg.pca$x[, 3] 

sdg.pca$rotation


library(ggplot2)
library(ggfortify)
autoplot(sdg.pca, data = basin_sdg, 
         colour = 'biogeographic_realm', 
         loadings = TRUE, 
         loadings.label = TRUE, 
         loadings.colour = "black", 
         loadings.label.colour = "black", 
         loadings.size =2, 
         loadings.label.repel = TRUE,
         label.repel = TRUE,
         label = TRUE, 
         max.overlaps = 5,
         label.label = "basin_id",
         label.colour='biogeographic_realm', 
         label.alpha=0.5, 
         label.position=position_jitter(width=0.012,height=0.012), 
         shape=19, alpha=0.7, size=2, label.size = 2, lwd=2) + 
  scale_colour_manual(values = c(
    "Afrotropic" = "#D55E00",
    "Australasia" = "#0072B2",
    "Indomalayan" = "#009E73",
    "Nearctic" = "#F0E442",
    "Neotropic" = "#E69F00",
    "Oceania" = "#56B4E9",
    "Palearctic" = "#CC79A7"
  )) + 
  theme_minimal()+
  #scale_x_continuous(limits = c(-0.12,0.12))+
  guides(color = guide_legend(nrow = 1))+
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.position = "top",
        legend.key.size = unit(0.2,"cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7))
ggsave("figures/supplement/Figure_S23.png",dpi = 300, units = "cm", width = 18, height = 16)

autoplot(sdg.pca, data = basin_sdg, 
         x = 1, y = 3,
         colour = 'biogeographic_realm', 
         loadings = TRUE, 
         loadings.label = TRUE, 
         loadings.colour = "black", 
         loadings.label.colour = "black", 
         loadings.size =2, 
         loadings.label.repel = TRUE,
         label.repel = TRUE,
         label = TRUE, 
         max.overlaps = 5,
         label.label = "basin_id",
         label.colour='biogeographic_realm', 
         label.alpha=0.5, 
         label.position=position_jitter(width=0.012,height=0.012), 
         shape=19, alpha=0.7, size=2, label.size = 2, lwd=2) + 
  scale_colour_manual(values = c(
    "Afrotropic" = "#D55E00",
    "Australasia" = "#0072B2",
    "Indomalayan" = "#009E73",
    "Nearctic" = "#F0E442",
    "Neotropic" = "#E69F00",
    "Oceania" = "#56B4E9",
    "Palearctic" = "#CC79A7"
  )) + 
  theme_minimal()+
  #scale_x_continuous(limits = c(-0.12,0.12))+
  guides(color = guide_legend(nrow = 1))+
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.position = "top",
        legend.key.size = unit(0.2,"cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7))


# We adopt a multi-objective prioritization framework integrating biodiversity irreplaceability, affordability, and feasibility. 
# Biodiversity irreplaceability is captured by the magnitude of Linnaean, Wallacean, and Darwinian shortfalls. 
# Affordability reflects the integrated cost of taxonomic description, field sampling, and genetic sequencing. 
# Feasibility is represented by two orthogonal socio-institutional gradients: 
#   PC1 (human pressure and socioeconomic development) and PC2 (governance capacity and environmental stewardship).
# Feasibility is represented by two orthogonal gradients: Socioeconomics and Governance.

# Knowledge Gap (Shortfall/NSCI)
# Investment (Cost)
# Socioeconomics (PC1)
# Governance (PC2)

################################################################################
library(dplyr)
library(sf)
library(readxl)
options(sf_use_s2 = FALSE)   # Disable spherical geometry to avoid topology issues
sdg <- read.csv("https://hub.arcgis.com/api/v3/datasets/eae2def2bf3c4c01b9ff358fcd1a5a13_0/downloads/data?format=csv&spatialRefId=3857&where=1%3D1")
names(sdg)
merit <- c("iso3","Name","Overall_Score","Overall_Rank","Region",
           paste0("Goal_",seq(1:17),"_Score"))
sdg <- sdg[,merit]
group <- read.csv("input/raw/sdg/Metadata_Country_API_NY.GDP.PCAP.CD_DS2_en_csv_v2_134819.csv")
sdg <- sdg %>% left_join(group[,c(1,3)], by = c("iso3"="Country.Code"))
# write.csv(sdg, "input/processed//sdg_2024.csv",row.names = F)

country_sdg <- sdg %>% select(iso3,starts_with("Goal_"))

names(country_sdg) <- gsub("Goal","SDG",names(country_sdg))
names(country_sdg) <- gsub("_Score","",names(country_sdg))

# Handling missing SDG indicators
# 
# Missing values in national SDG indicators (SDG 1–17) were conservatively imputed 
# using a rule-based hierarchical strategy. For each SDG indicator, missing values 
# were first replaced by the median of countries within the same income group and 
# continent. If no valid observations were available at this level, we sequentially
# fell back to medians calculated within the same income group, within the same continent,
# and finally across all countries. This approach avoids model-based assumptions while 
# preserving large-scale socioeconomic structure.
country_list <- read.csv("input/raw/country_list.csv")
country_sdg <- country_sdg %>% left_join(country_list[,c(1,4,5),],by = "iso3") 


sdg_cols <- paste0("SDG_", 1:17)

country_sdg_final <- country_sdg %>%
  mutate(
    income_group = if_else(is.na(income_group), "Unknown", income_group),
    continent    = if_else(is.na(continent), "Unknown", continent)
  ) %>%
  group_by(income_group, continent) %>%
  mutate(across(all_of(sdg_cols),
                ~ if_else(is.na(.x), median(.x, na.rm = TRUE), .x))) %>%
  ungroup() %>%
  group_by(income_group) %>%
  mutate(across(all_of(sdg_cols),
                ~ if_else(is.na(.x), median(.x, na.rm = TRUE), .x))) %>%
  ungroup() %>%
  group_by(continent) %>%
  mutate(across(all_of(sdg_cols),
                ~ if_else(is.na(.x), median(.x, na.rm = TRUE), .x))) %>%
  ungroup() %>%
  mutate(across(all_of(sdg_cols),
                ~ if_else(is.na(.x), median(.x, na.rm = TRUE), .x))) %>%
  mutate(
    income_group = na_if(income_group, "Unknown"),
    continent    = na_if(continent, "Unknown")
  )

#saveRDS(country_sdg_final,"input/processed/sdg_country.rds")

################################################################################
country_sdg <- readRDS("input/processed/sdg_country.rds")

sdg.pca <- prcomp(country_sdg[,2:18], center = TRUE,scale. = TRUE) 
sdg.pca$PC1 <- sdg.pca$x[, 1] 
sdg.pca$PC2 <- sdg.pca$x[, 2] 
sdg.pca$PC3 <- sdg.pca$x[, 3] 

sdg.pca$rotation


library(ggplot2)
library(ggfortify)
autoplot(sdg.pca, data = country_sdg, 
         colour = 'continent', 
         loadings = TRUE, 
         loadings.label = TRUE, 
         loadings.colour = "black", 
         loadings.label.colour = "black", 
         loadings.size =2, 
         loadings.label.repel = TRUE,
         label.repel = TRUE,
         label = TRUE, 
         max.overlaps = 5,
         label.label = "iso3",
         label.colour='continent', 
         label.alpha=0.5, 
         label.position=position_jitter(width=0.012,height=0.012), 
         shape=19, alpha=0.7, size=2, label.size = 2, lwd=2) + 
  scale_colour_manual(values = c(
    "Africa"   = "#D55E00",
    "Oceania" = "#0072B2",
    "Asia"  = "#009E73",
    "North America"     = "#F0E442",
    "South America"    = "#E69F00",
    "Europe"   = "#CC79A7"
  ),na.translate = FALSE) + 
  theme_minimal()+
  #scale_x_continuous(limits = c(-0.12,0.12))+
  guides(color = guide_legend(nrow = 1))+
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.position = "top",
        legend.key.size = unit(0.2,"cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7))
ggsave("figures/supplement/Figure_S24.png",dpi = 300, units = "cm", width = 18, height = 16)

autoplot(sdg.pca, data = country_sdg, 
         x = 1, y = 3,
         colour = 'continent', 
         loadings = TRUE, 
         loadings.label = TRUE, 
         loadings.colour = "black", 
         loadings.label.colour = "black", 
         loadings.size =2, 
         loadings.label.repel = TRUE,
         label.repel = TRUE,
         label = TRUE, 
         max.overlaps = 5,
         label.label = "basin_id",
         label.colour='continent', 
         label.alpha=0.5, 
         label.position=position_jitter(width=0.012,height=0.012), 
         shape=19, alpha=0.7, size=2, label.size = 2, lwd=2) + 
  scale_colour_manual(values = c(
    "Africa"   = "#D55E00",
    "Oceania" = "#0072B2",
    "Asia"  = "#009E73",
    "North America"     = "#F0E442",
    "South America"    = "#E69F00",
    "Europe"   = "#CC79A7"
  ),na.translate = FALSE) + 
  theme_minimal()+
  #scale_x_continuous(limits = c(-0.12,0.12))+
  guides(color = guide_legend(nrow = 1))+
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 9),
        legend.position = "top",
        legend.key.size = unit(0.2,"cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7))
