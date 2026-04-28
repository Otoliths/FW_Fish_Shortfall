# ------------------------------------------------------------------------------
# Supplementary Figure S17 Potential predictor variables
# ------------------------------------------------------------------------------

rm(list = ls())
library(dplyr)
library(ggplot2)
library(magrittr)
library(data.table)
library(arrow)
library(patchwork)
library(ggh4x)

body_size <- readRDS("input/data_prep/body_size_final_single.rds")####
preserved_specimen <- readRDS("input/data_prep/preserved_specimens.rds")####

range_size <- read_parquet("input/data_prep/basin_range_size.parquet")
rarity <- read_parquet("input/data_prep/basin_rarity.parquet")
latitude <- read_parquet("input/data_prep/basin_latitude_avg.parquet")
watershed_area <- readRDS("input/data_prep/watershed_area.rds")

discharge <- read_parquet("input/data_prep/basin_discharge_avg_year.parquet")
watertemp <- read_parquet("input/data_prep/basin_watertemp_avg_year.parquet")

elevation <- read_parquet("input/data_prep/basin_elevation.parquet")
population_density <- read_parquet("input/data_prep/basin_population_density.parquet")

taxonomic_effort <- readRDS("input/data_prep/taxonomic_effort_year.rds")###
taxonomic_activity <- readRDS("input/data_prep/taxonomic_activity_year.rds")###
sampling_effort <- read_parquet("input/data_prep/basin_sampling_effort.parquet")
sequencing_effort <- read_parquet("input/data_prep/basin_sequencing_effort.parquet")



p01 <- ggplot(data = body_size, aes(log_length_max)) +
  #ggrastr::rasterize(geom_rug(col = "#6D4D80", alpha = 0.5), dpi = 300)+
  geom_density(fill = "#6C8F5D80", colour = "#6C8F5D80", alpha = 0.5)+
  xlab(expression(log[10] * "(Body size (cm))"))+
  ylab("Density") +
  theme_classic()+
  annotate("text",x = Inf, y = Inf,label = "Biology",hjust = 1.05, colour = "#6C8F5D80",
           vjust = 1.2,size = 2.5,fontface = "bold")

p02 <- ggplot(data = preserved_specimen, aes(preserved_specimen+1)) +
  #ggrastr::rasterize(geom_rug(col = "#6D4D80", alpha = 0.5), dpi = 300)+
  geom_density(fill = "#6C8F5D80", colour = "#6C8F5D80", alpha = 0.5)+
  scale_x_log10(guide = "axis_logticks",
                breaks = c(1,100,100000,100000000),                
                labels = expression(10^0, 10^2, 10^5, 10^8)
  )+
  xlab("N. preserved specimen")+
  ylab("Density") +
  theme_classic()+
  annotate("text",x = Inf, y = Inf,label = "Biology",hjust = 1.05, colour = "#6C8F5D80",
           vjust = 1.2,size = 2.5,fontface = "bold")


p03 <- ggplot(data = range_size, aes(log10(range_size+1))) +
  #ggrastr::rasterize(geom_rug(col = "#6D4D80", alpha = 0.5), dpi = 300)+
  geom_density(fill = "#6D4D8080", colour = "#6D4D8080", alpha = 0.5)+
  xlab(expression(log[10] * "(Range size)")) + 
  ylab("Density") +
  theme_classic()+
  annotate("text",x = Inf, y = Inf,label = "Geography",hjust = 1.05, colour = "#6D4D8080",
           vjust = 1.2,size = 2.5,fontface = "bold")

p04 <- rarity %>%
  group_by(basin_id, sampling_year) %>%
  summarise(
    rarity_longterm = mean(relative_rarity_cummean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = sampling_year, y = rarity_longterm, group = basin_id)) +
  ggrastr::rasterize(geom_line(alpha = 0.1, linewidth = 0.2, na.rm = TRUE,colour = "#6D4D8080"),dpi=300) +
  xlab("Year") +
  ylab("Range rarity") +
  scale_x_continuous(breaks = c(1800,1900,2020),limits = c(1758,2024))+
  theme_classic()+
  annotate("text",x = Inf, y = Inf,label = "Geography",hjust = 1.05, colour = "#6D4D8080",
           vjust = 1.2,size = 2.5,fontface = "bold")

p05 <- ggplot(data = latitude, aes(latitude)) +
  #ggrastr::rasterize(geom_rug(col = "#6D4D80", alpha = 0.5), dpi = 300)+
  geom_density(fill = "#6D4D8080", colour = "#6D4D8080", alpha = 0.5)+
  xlab("Latitude") +
  ylab("Density") +
  theme_classic()+
  annotate("text",x = Inf, y = Inf,label = "Geography",hjust = 1.05, colour = "#6D4D8080",
           vjust = 1.2,size = 2.5,fontface = "bold")

p06 <- ggplot(data = watershed_area, aes(log10(area_km2))) +
  #ggrastr::rasterize(geom_rug(col = "#6D4D80", alpha = 0.5), dpi = 300)+
  geom_density(fill = "#6D4D8080", colour = "#6D4D8080", alpha = 0.5)+
  xlab(expression(log[10] * "(Watershed area (km"^2*"))")) +
  ylab("Density") +
  theme_classic()+
  annotate("text",x = Inf, y = Inf,label = "Geography",hjust = 1.05, colour = "#6D4D8080",
           vjust = 1.2,size = 2.5,fontface = "bold")


p07 <- ggplot(discharge, aes(x = year, y = discharge, group = basin_id)) +
  ggrastr::rasterize(geom_line(alpha = 0.02, linewidth = 0.2, na.rm = TRUE,colour = "#A3473E80"),dpi = 300) +
  xlab("Year") +
  ylab(expression("Streamflow (m"^3*"/s)")) +
  theme_classic() +
  scale_y_log10(
    breaks = c(0.001, 1, 10, 100, 1000),
    labels = function(x) {
      ifelse(
        x %in% c(0, 1, 10, 100, 1000),
        parse(text = paste0("10^", log10(x))),
        as.character(x)
      )
    }
  )+
  scale_x_continuous(limits = c(1976,2024), breaks = c(1980,1990,2000,2010,2020))+
  annotate("text",x = Inf, y = Inf,label = "Environment",hjust = 1.1, colour = "#A3473E80",
           vjust = 1.2,size = 2.5,fontface = "bold")


p08 <- ggplot(data = watertemp, aes(x = year, y = watertemp-273.15, group = basin_id)) +
  ggrastr::rasterize(geom_line(alpha = 0.02, linewidth = 0.2, na.rm = TRUE,colour = "#A3473E80"),dpi=300) +
  xlab("Year") +
  ylab(expression("Water temperature ("*degree*C*")")) +
  theme_classic()+
  scale_x_continuous(limits = c(1976,2024), breaks = c(1980,1990,2000,2010,2020))+
  annotate("text",x = Inf, y = Inf,label = "Environment",hjust = 1.1, colour = "#A3473E80",
           vjust = 1.2,size = 2.5,fontface = "bold")


p09 <- ggplot(data = elevation, aes(elevation)) +
  #ggrastr::rasterize(geom_rug(col = "#6D4D80", alpha = 0.5), dpi = 300)+
  geom_density(fill = "#4E6E8E80", colour = "#4E6E8E80", alpha = 0.5)+
  xlab("Elevation (m)") +
  ylab("Density") +
  theme_classic()+
  annotate("text",x = Inf, y = Inf,label = "Access",hjust = 1.05, colour = "#4E6E8E80",
           vjust = 1.2,size = 2.5,fontface = "bold")

p10 <- population_density %>%
  group_by(basin_id, year) %>%
  summarise(
    population_density = mean(population_density, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = year, y = log10(population_density+1), group = basin_id)) +
  ggrastr::rasterize(geom_line(alpha = 0.1, linewidth = 0.2, na.rm = TRUE,colour = "#4E6E8E80"),dpi=300) +
  xlab("Year") +
  ylab(expression(log[10] * "(Human density)")) +
  scale_x_continuous(breaks = c(1800,1900,2020),limits = c(1758,2024))+
  theme_classic()+
  annotate("text",x = Inf, y = Inf,label = "Access",hjust = 1.05, colour = "#4E6E8E80",
           vjust = 1.2,size = 2.5,fontface = "bold")

p11 <- ggplot(data = taxonomic_effort, aes(x= year, y = TAE,group = as.factor(valid_name))) +
  ggrastr::rasterize(geom_line(alpha = 0.05, linewidth = 0.05, na.rm = TRUE,colour = "#8A6A3F80"),dpi = 300) +
  xlab("Year") + 
  ylab("Taxonomic effort") +
  theme_classic()+
  scale_x_continuous(breaks = c(1800,1900,2020),limits = c(1758,2024))+
  annotate("text",x = Inf, y = Inf,label = "Activity",hjust = 1.05, colour = "#8A6A3F80",
           vjust = 1.2,size = 2.5,fontface = "bold")


p12 <- ggplot(data = taxonomic_activity, aes(x= year, y = TAA,group = as.factor(valid_name))) +
  ggrastr::rasterize(geom_line(alpha = 0.05, linewidth = 0.05, na.rm = TRUE,colour = "#8A6A3F80"),dpi = 300) +
  xlab("Year") + 
  ylab("Taxonomic activity") +
  theme_classic()+
  scale_x_continuous(breaks = c(1800,1900,2020),limits = c(1758,2024))+
  annotate("text",x = Inf, y = Inf,label = "Activity",hjust = 1.05, colour = "#8A6A3F80",
           vjust = 1.2,size = 2.5,fontface = "bold")

p13 <- sampling_effort %>%
  group_by(basin_id, sampling_year) %>%
  summarise(
    sampling_incidence = sum(sampling_effort, na.rm = TRUE),  # 
    .groups = "drop"
  ) %>%
  arrange(basin_id, sampling_year) %>%
  group_by(basin_id) %>%
  mutate(
    sampling_cum = cumsum(sampling_incidence)
  ) %>%
  ungroup() %>%
  ggplot(aes(x = sampling_year, y = sampling_cum, group = basin_id)) +
  ggrastr::rasterize(geom_line(alpha = 0.1, linewidth = 0.2, na.rm = TRUE,colour = "#8A6A3F80"),dpi=300) +
  xlab("Year") +
  ylab("Sampling effort") +
  scale_x_continuous(breaks = c(1800,1900,2020),limits = c(1758,2024))+
  scale_y_log10(
    breaks = c(0.001, 1, 10, 100, 1000,10000),
    labels = function(x) {
      ifelse(
        x %in% c(0, 1, 10, 100, 1000,10000),
        parse(text = paste0("10^", log10(x))),
        as.character(x)
      )
    }
  )+
  theme_classic()+
  annotate("text",x = Inf, y = Inf,label = "Activity",hjust = 1.05, colour = "#8A6A3F80",
           vjust = 1.2,size = 2.5,fontface = "bold")


p14 <- sequencing_effort %>%
  group_by(basin_id, sequencing_year) %>%
  summarise(
    sequencing_incidence = sum(sequencing_effort, na.rm = TRUE),  # 
    .groups = "drop"
  ) %>%
  arrange(basin_id, sequencing_year) %>%
  group_by(basin_id) %>%
  mutate(
    sequencing_cum = cumsum(sequencing_incidence)
  ) %>%
  ungroup() %>%
  ggplot(aes(x = sequencing_year, y = sequencing_cum, group = basin_id)) +
  ggrastr::rasterize(geom_line(alpha = 0.1, linewidth = 0.2, na.rm = TRUE,colour = "#8A6A3F80"),dpi=300) +
  xlab("Year") +
  ylab("Sequencing effort") +
  scale_x_continuous(breaks = c(1990,2000,2010,2020),limits = c(1985,2024))+
  scale_y_log10(
    breaks = c(0.001, 1, 10, 100, 1000,10000),
    labels = function(x) {
      ifelse(
        x %in% c(0, 1, 10, 100, 1000,10000),
        parse(text = paste0("10^", log10(x))),
        as.character(x)
      )
    }
  )+
  theme_classic()+
  annotate("text",x = Inf, y = Inf,label = "Activity",hjust = 1.05, colour = "#8A6A3F80",
           vjust = 1.2,size = 2.5,fontface = "bold")


library(cowplot)
plots <- list(
  p01, p02, p03, p04,
  p05, p06, p07, p08,
  p09, p10, p11, p12,
  p13, p14
)
base_theme <- theme(plot.margin = margin(0.4,0.4,0.4,0.4),
                    axis.title = element_text(size = 8),
                    axis.text = element_text(size = 6))

plots <- lapply(plots, `+`, base_theme)
plot_grid(
  plotlist   = plots,
  nrow       = 4,
  labels     = LETTERS,
  label_fontface = "bold",
  label_colour   = "black",
  label_size     = 9,
  label_x        = 0.02,
  label_y        = 0.98,
  hjust          = 0,
  vjust          = 1
)



ggsave("figures/supplement/Figure_S17.png",dpi = 300, units = "cm", width = 21, height = 20)

