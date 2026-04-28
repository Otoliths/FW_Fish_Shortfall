# ------------------------------------------------------------------------------
# Supplementary Figure S11 
# National-scale patterns of biogeographical ignorance in freshwater fish biodiversity
# ------------------------------------------------------------------------------
library(dplyr)
#library(ggpubr)
library(cowplot)
library(tidyr)
library(Hmisc)
library(rnaturalearth)
library(ggplot2)
library(biscale)
#library(stringr)
library(ggtext)
library(brms)         # Bayesian multilevel models interface
library(arrow)
library(rstan)
library(tidybayes)
library(sf)
library(viridis)
library(scales)
library(legendry) #https://github.com/teunbrand/legendry/issues/80
sf_use_s2(FALSE)

#source("reports/zzz.r")

linnaean <- read.csv("output/tables/country_linnaean_shortfall.csv")
wallacean <- read.csv("output/tables/country_wallacean_shortfall.csv")
darwinian <- read.csv("output/tables/country_darwinian_shortfall.csv")

df_country <- linnaean[,c(1,10)] %>% 
  left_join(wallacean[,c(1,6)], by = "iso3") %>%
  left_join(darwinian[,c(1,6)], by = "iso3")

country <- read.csv("input/raw/country_list.csv")
df <- df_country %>% left_join(country,by = "iso3")
names(df)[2:4] <- c("LS","WS","DS")
head(df)
df$continent <- as.factor(df$continent)
df <- df %>% na.omit()

library(dplyr)
classify_within_continent_quantile <- function(data,
                                               value_col,
                                               realm_col = "continent",
                                               q = 0.25) {
  stopifnot(value_col %in% names(data), realm_col %in% names(data))
  stopifnot(q > 0 & q < 0.5)
  
  data %>%
    group_by(.data[[realm_col]]) %>%
    mutate(
      hi = quantile(.data[[value_col]], 1 - q, na.rm = TRUE, type = 7),
      lo = quantile(.data[[value_col]], q,     na.rm = TRUE, type = 7),
      state = case_when(
        .data[[value_col]] >= hi ~ "high",
        .data[[value_col]] <= lo ~ "low",
        TRUE ~ NA_character_
      )
    ) %>%
    ungroup() %>%
    pull(state)
}
# 
# library(purrr)
# library(tidyr)
# qs <- seq(0.01, 0.30, by = 0.01)
# 
# scenario_by_q <- map_dfr(qs, function(qi) {
#   
#   df_tmp <- df %>%
#     mutate(
#       LS_state = classify_within_continent_quantile(., "LS", q = qi),
#       WS_state = classify_within_continent_quantile(., "WS", q = qi),
#       DS_state = classify_within_continent_quantile(., "DS", q = qi),
#       scenario = case_when(
#         LS_state == "high" & WS_state == "high" & DS_state == "high" ~ "LS↑ WS↑ DS↑",
#         LS_state == "high" & WS_state == "high" & DS_state == "low"  ~ "LS↑ WS↑ DS↓",
#         LS_state == "high" & WS_state == "low"  & DS_state == "high" ~ "LS↑ WS↓ DS↑",
#         LS_state == "high" & WS_state == "low"  & DS_state == "low"  ~ "LS↑ WS↓ DS↓",
#         LS_state == "low"  & WS_state == "high" & DS_state == "high" ~ "LS↓ WS↑ DS↑",
#         LS_state == "low"  & WS_state == "high" & DS_state == "low"  ~ "LS↓ WS↑ DS↓",
#         LS_state == "low"  & WS_state == "low"  & DS_state == "high" ~ "LS↓ WS↓ DS↑",
#         TRUE ~ NA_character_
#       )
#     )
#   
#   df_tmp %>%
#     count(scenario) %>%
#     mutate(q = qi)
# })

# scenario_by_q %>%
#   filter(!is.na(scenario)) %>%
#   arrange(q, desc(n))
# 
# library(ggplot2)
# 
# ggplot(
#   scenario_by_q %>% filter(!is.na(scenario)),
#   aes(x = q, y = n, colour = scenario)
# ) +
#   geom_line() +
#   labs(
#     x = "Quantile threshold (q)",
#     y = "Number of countries",
#     colour = "Scenario"
#   ) +
#   theme_minimal()

#Quantile thresholds for scenario classification were selected based on analytical scale. 
#For drainage basins, upper and lower quartiles (q = 0.25) were used to ensure broad spatial coverage. 
#For country-level analyses, a slightly more conservative threshold (q = 0.30) was adopted to account for 
#smaller sample sizes and to avoid over-expansion of dominant scenarios. Sensitivity analyses across alternative thresholds yielded consistent patterns.

df_q30 <- df %>%
  mutate(
    LS_state = classify_within_continent_quantile(., "LS", q = 0.3),
    WS_state = classify_within_continent_quantile(., "WS", q = 0.3),
    DS_state = classify_within_continent_quantile(., "DS", q = 0.3)
  )
df_q30 <- df_q30 %>%
  mutate(
    scenario = case_when(
      LS_state == "high" & WS_state == "high" & DS_state == "high" ~ "LS↑ WS↑ DS↑",
      LS_state == "high" & WS_state == "high" & DS_state == "low"  ~ "LS↑ WS↑ DS↓",
      LS_state == "high" & WS_state == "low"  & DS_state == "high" ~ "LS↑ WS↓ DS↑",
      LS_state == "high" & WS_state == "low"  & DS_state == "low"  ~ "LS↑ WS↓ DS↓",
      LS_state == "low"  & WS_state == "high" & DS_state == "high" ~ "LS↓ WS↑ DS↑",
      LS_state == "low"  & WS_state == "high" & DS_state == "low"  ~ "LS↓ WS↑ DS↓",
      LS_state == "low"  & WS_state == "low"  & DS_state == "high" ~ "LS↓ WS↓ DS↑",
      TRUE ~ NA_character_
    )
  )


#write.csv(df_q30,"output/tables/country_shortfalls.csv",row.names = F)
country_shortfalls <- read.csv("output/tables/country_shortfalls.csv" )

#world_map <- ne_countries(scale = 50, type = "countries", returnclass = "sf")
moll_proj <- st_crs("+proj=moll")
lat_points <- data.frame(lon = 0, lat = c(-60, 90)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
projected_points <- st_transform(lat_points, crs = moll_proj)
y_limits <- st_coordinates(projected_points)[, "Y"]


hotspots <- st_read("input/raw/hotspots_2016_1/hotspots_2016_1_fix.shp")
inland <- readRDS("input/raw/country_20251212.rds")
scenario_sf <- inland %>% 
  left_join(country_shortfalls[,c(1,12)], by = "iso3") %>%
  st_wrap_dateline(options = c("WRAPDATELINE=YES"))

hotspots <- st_transform(hotspots, st_crs(scenario_sf))

scenario_sf <- scenario_sf %>%
  mutate(
    is_hotspot = ifelse(lengths(st_intersects(geometry, hotspots %>% filter(Type != "outer limit"))) > 0,TRUE,FALSE)
  )

overlap_tab <- scenario_sf %>%
  st_drop_geometry() %>%
  filter(!is.na(scenario)) %>%  
  group_by(is_hotspot, scenario) %>%
  summarise(
    n_country = n(),
    .groups = "drop"
  ) %>%
  group_by(is_hotspot) %>%
  mutate(
    prop = n_country/ sum(n_country)
  )

overlap_tab
# A tibble: 6 × 4
# Groups:   is_hotspot [2]
# is_hotspot scenario    n_country   prop
# <lgl>      <chr>           <int>  <dbl>
# 1 FALSE      LS↑ WS↑ DS↑         4 0.667 
# 2 FALSE      LS↓ WS↑ DS↑         2 0.333 
# 3 TRUE       LS↑ WS↑ DS↑        26 0.897 
# 4 TRUE       LS↑ WS↑ DS↓         1 0.0345
# 5 TRUE       LS↑ WS↓ DS↓         1 0.0345
# 6 TRUE       LS↓ WS↑ DS↑         1 0.0345


#overlap_tab <- rbind(overlap_tab,remain)
scenario_sf <- scenario_sf %>%
  mutate(
    is_shortfall = ifelse(is.na(scenario),FALSE,TRUE)
  )
scenario_sf$is_hotspot <- as.factor(scenario_sf$is_hotspot)
scenario_sf$is_shortfall <- as.factor(scenario_sf$is_shortfall)
scenario_sf <- bi_class(scenario_sf,x=is_shortfall,y= is_hotspot,dim = 2)

custom_pal <- c(
  "1-1" = "#d3d3d3", # low x, low y
  "2-1" = "#4A3A3B", # high x, low y
  "1-2" = "#F16C31", # low x, high y
  "2-2" = "#7C0809" # high x, high y
)

leg <- bi_pal(pal = "PurpleOr", dim = 2)
grid_dat <- leg$layers[[1]][["data"]]

pts <- scenario_sf %>% st_drop_geometry() %>%
  left_join(
    leg$layers[[1]][["data"]] %>% select(bi_class, x, y),
    by = "bi_class"
  )
count_df <- pts %>% 
  count(bi_class, name = "n_cells") %>%
  right_join(
    grid_dat %>% select(bi_class, x, y),
    by = "bi_class"
  ) %>%
  replace_na(list(n_cells = 0)) %>%
  mutate(n_cells_ratio = n_cells/sum(n_cells))

grid_dat$bi_fill <- custom_pal

bi_legend <- ggplot() +
  geom_tile(data = grid_dat,aes(x = x, y = y, fill = bi_fill),colour = "white") +
  scale_fill_identity() +
  geom_text(data = count_df,aes(x = x, y = y, label = percent(n_cells_ratio, accuracy = 1)),
            colour   = "white",fontface = "bold", size = 2.2)+
  scale_x_continuous(breaks = 1:2,labels = c("N", "Y"),expand = c(0.015, 0.015)) +
  scale_y_continuous(breaks = 1:2,labels = c("N", "Y"),expand = c(0.015, 0.015)) +
  labs(
    x = expression(Shortfall), 
    y = expression(Hotspot)
  ) +
  coord_equal(clip = "off") +
  # annotate("text",x = Inf, y = 4.8, hjust = 1.2, vjust = 0,size  = 2.2,label = deparse(label_expr),parse = TRUE)+
  theme_void(base_size = 8) +
  theme(
    axis.title  = element_text(size = 8, colour = "black", face = "bold"),
    axis.text   = element_text(size = 6, colour = "black"),
    axis.title.y.left = element_text(angle = 90),
    axis.text.y.left = element_text(angle = 90),
    plot.margin = margin(0,0,0,0),
    panel.ontop = TRUE
  )

map <- scenario_sf %>% 
  #filter(scenario == "LS↑ WS↑ DS↑") %>%
  ggplot() +
  #ggrastr::rasterise(geom_sf(data = world_map, fill = "lightgrey", colour = NA),dpi = 300)+
  ggrastr::rasterise(geom_sf(aes(fill = bi_class),colour = "white", linewidth = 0.03, show.legend = FALSE),dpi = 300)+
  bi_scale_fill(pal = custom_pal, dim = 2) +
  bi_theme()+
  coord_sf(crs = "+proj=moll +lon_0=0", ylim = y_limits, expand = FALSE)+
  theme(plot.margin = margin(0,0,0,0))

p2 <- ggdraw() +
  draw_plot(map, 0, 0, 1, 1) +
  draw_plot(bi_legend, 0.09, 0.16, 0.34, 0.34)


#################################################################################
count_scenario_1 <- pts %>%
  select(is_hotspot,scenario) %>%
  mutate(scenario = ifelse(scenario == "LS↑ WS↑ DS↑","LS↑ WS↑ DS↑",NA),
         is_shortfall = ifelse(is.na(scenario),FALSE,TRUE)) %>%
  mutate(is_shortfall = factor(is_shortfall),
         is_hotspot = factor(is_hotspot)) %>%
  bi_class(x=is_shortfall,y= is_hotspot,dim = 2) %>%
  group_by(scenario,bi_class) %>%
  summarise(n_cells = n(),.groups = "drop") %>%
  mutate(scenario = 'LS* "" %up% "" *WS* "" %up% "" *DS* "" %up% ""')

count_scenario_2 <- pts %>%
  select(is_hotspot,scenario) %>%
  mutate(scenario = ifelse(scenario == "LS↑ WS↑ DS↓","LS↑ WS↑ DS↓",NA),
         is_shortfall = ifelse(is.na(scenario),FALSE,TRUE)) %>%
  mutate(is_shortfall = factor(is_shortfall),
         is_hotspot = factor(is_hotspot)) %>%
  bi_class(x=is_shortfall,y= is_hotspot,dim = 2) %>%
  group_by(scenario,bi_class) %>%
  summarise(n_cells = n(),.groups = "drop") %>%
  rbind(data.frame(scenario = NA,bi_class = "2-1",n_cells = 0)) %>%
  mutate(scenario = 'LS* "" %up% "" *WS* "" %up% "" *DS* "" %down% ""')

count_scenario_4 <- pts %>%
  select(is_hotspot,scenario) %>%
  mutate(scenario = ifelse(scenario == "LS↑ WS↓ DS↓","LS↑ WS↓ DS↓",NA),
         is_shortfall = ifelse(is.na(scenario),FALSE,TRUE)) %>%
  mutate(is_shortfall = factor(is_shortfall),
         is_hotspot = factor(is_hotspot)) %>%
  bi_class(x=is_shortfall,y= is_hotspot,dim = 2) %>%
  group_by(scenario,bi_class) %>%
  summarise(n_cells = n(),.groups = "drop") %>%
  rbind(data.frame(scenario = NA,bi_class = "2-1",n_cells = 0)) %>%
  mutate(scenario = 'LS* "" %up% "" *WS* "" %down% "" *DS* "" %down% ""')   

count_scenario_5 <- pts %>%
  select(is_hotspot,scenario) %>%
  mutate(scenario = ifelse(scenario == "LS↓ WS↑ DS↑","LS↓ WS↑ DS↑",NA),
         is_shortfall = ifelse(is.na(scenario),FALSE,TRUE)) %>%
  mutate(is_shortfall = factor(is_shortfall),
         is_hotspot = factor(is_hotspot)) %>%
  bi_class(x=is_shortfall,y= is_hotspot,dim = 2) %>%
  group_by(scenario,bi_class) %>%
  summarise(n_cells = n(),.groups = "drop")  %>%
  mutate(scenario = 'LS* "" %down% "" *WS* "" %up% "" *DS* "" %up% ""')    



count_scenario <- rbind(count_scenario_1,count_scenario_2,count_scenario_4,count_scenario_5)
rm(count_scenario_1,count_scenario_2,count_scenario_4,count_scenario_5)

count_scenario <- count_scenario %>%
  left_join(grid_dat,by = "bi_class")


count_scenario$scenario <- factor(count_scenario$scenario,
                                  levels = c('LS* "" %up% "" *WS* "" %up% "" *DS* "" %up% ""',
                                             'LS* "" %up% "" *WS* "" %up% "" *DS* "" %down% ""',
                                             'LS* "" %up% "" *WS* "" %down% "" *DS* "" %down% ""',
                                             'LS* "" %down% "" *WS* "" %up% "" *DS* "" %up% ""'
                                  ))


p1 <- ggplot(data = count_scenario) +
  geom_tile(aes(x = x, y = y, fill = bi_fill),colour = "white") +
  scale_fill_identity() +
  geom_text(aes(x = x, y = y, label = format(n_cells, big.mark=",")),
            colour   = "white",fontface = "bold", size = 2,hjust = 0.5, vjust = 0.5)+
  scale_x_continuous(breaks = 1:2,labels = c("N", "Y"),expand = c(0.015, 0.015)) +
  scale_y_continuous(breaks = 1:2,labels = c("N", "Y"),expand = c(0.015, 0.015)) +
  labs(
    x = "Combinations of biodiversity shortfalls", 
    y = expression(Hotspot)
  ) +
  coord_equal() +
  facet_wrap(~ scenario,nrow = 1, labeller = label_parsed)+
  theme_void()+
  theme(
    axis.title  = element_text(size = 8, colour = "black"),
    axis.text   = element_text(size = 6, colour = "black"),
    axis.title.y.left = element_text(angle = 90),
    axis.text.y.left = element_text(angle = 90),
    strip.text = element_text(size = 6),
    panel.spacing = unit(0, "lines"),
    plot.margin = margin(0,0,0,0),
    panel.ontop = TRUE)

################################################################################
linnaean <- read.csv("output/tables/country_linnaean_shortfall.csv")
wallacean <- read.csv("output/tables/country_wallacean_shortfall.csv")
darwinian <- read.csv("output/tables/country_darwinian_shortfall.csv")

df_country <- linnaean[,c(1,10)] %>% 
  left_join(wallacean[,c(1,6)], by = "iso3") %>%
  left_join(darwinian[,c(1,6)], by = "iso3")

min_max_normalize_safe <- function(x, epsilon = 1e-4) {
  min_val <- min(x,na.rm = T)
  max_val <- max(x,na.rm = T)
  range_val <- max_val - min_val + epsilon
  return((x - min_val) / range_val)
}


df_country$SRdesc <- log10(df_country$SRdesc+1)
df_country$SRnoloc <- log10(df_country$SRnoloc+1)
df_country$SRnoseq <- log10(df_country$SRnoseq+1)

df_country$SRdesc <- min_max_normalize_safe(df_country$SRdesc)
df_country$SRnoloc <- min_max_normalize_safe(df_country$SRnoloc)
df_country$SRnoseq <- min_max_normalize_safe(df_country$SRnoseq)


df_nsc <- df_country %>%
  mutate(
    SL = SRdesc,
    SW = SRnoloc,
    SD = SRnoseq
  ) %>%
  filter(!is.na(SL)) %>%
  mutate(
    SW_adj = (1 - SL) * SW, 
    SD_adj = (1 - SL) * SD,
    NSCI = SL + SW_adj + SD_adj
  )
# normalize 0–1
df_nsc <- df_nsc %>%
  mutate(
    NSCI_norm = (NSCI - min(NSCI, na.rm = TRUE)) /
      (max(NSCI, na.rm = TRUE) - min(NSCI, na.rm = TRUE))
  )


saveRDS(df_nsc,"input/processed/country_shortfall.rds")
# classification by percentiles
# q_high <- quantile(df_nsc$NSCI_rank_norm, 0.75, na.rm = TRUE)
# q_low  <- quantile(df_nsc$NSCI_rank_norm, 0.25, na.rm = TRUE)
# 
# df_nsc <- df_nsc %>%
#   mutate(
#     NSCI_class = case_when(
#       NSCI_rank_norm >= q_high ~ "Composite↑",
#       NSCI_rank_norm <= q_low  ~ "Composite↓",
#       TRUE ~ NA_character_
#     )
#   )

dt <- inland %>% left_join(df_nsc[,c(1,11)], by = "iso3") %>%
  st_wrap_dateline(options = c("WRAPDATELINE=YES"))


p3 <- ggplot(data = dt) +
  #ggrastr::rasterise(geom_sf(data = world_map, fill = "lightgrey", colour = NA),dpi = 300)+
  ggrastr::rasterise(geom_sf(aes(fill = NSCI_norm),colour = "white", linewidth = 0.03),dpi = 300)+
  scale_fill_viridis_c(
    name = "NSCI",
    na.value = "lightgrey",
    option = "C",
    direction = -1,
    begin = 0.1,
    end = 0.99,
    breaks = seq(0, 1, 0.5),
    labels = scales::number_format(accuracy = 0.1),
    guide = legendry::compose_sandwich(
      middle = gizmo_histogram(just = 1),
      text   = "axis_base"
    )
  )+
  bi_theme()+
  coord_sf(crs = "+proj=moll +lon_0=0", ylim = y_limits, expand = FALSE)+
  theme(
    legend.position = c(0.06, 0.13),
    legend.justification = c(0, 0),
    legend.background = element_blank(),
    legend.key.height = unit(0.5,"cm"),
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8,margin = margin(l = 4)),
    legend.axis.line  = element_line(linewidth = 0.2, color = "black")
  )


ggdraw() +
  draw_plot(p1, x = 0.18,   y = 2/3, width = 0.6, height = 0.13) +
  draw_label("A",x = 0.17,y = 2/3 + 0.14,
             fontface = "bold",size = 12,hjust = 0,vjust = 1) +
  draw_plot(p2, x = 0,   y = 1/3, width = 1.0, height = 0.33) +  
  draw_label("B",x = 0.17, y = 1/3 + 0.33 - 0.02,
             fontface = "bold",size = 12,hjust = 0,vjust = 1) +
  draw_plot(p3, x = 0,   y = 0,   width = 1.0, height = 0.33)+
  draw_label("C",x = 0.17,y = 0.31,
             fontface = "bold",size = 12,hjust = 0,vjust = 1)
ggsave("figures/supplement/Figure_S11.png",width = 20, height = 20, units = "cm",dpi = 300)
