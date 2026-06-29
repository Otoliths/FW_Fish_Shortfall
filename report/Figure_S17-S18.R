# Supplementary Figure S17 Temporal decline in documentation lags following taxonomic description
# Supplementary Figure S18 Global patterns and delays in first georeferenced freshwater fish records following taxonomic description

rm(list = ls())
library(ggplot2)
library(dplyr)
library(scales)
library(ggpmisc)
library(patchwork)
basin_wallacean <- readRDS("output/stan_survdata_basin_wallacean_single.rds")

taxon_geo_lag_clean <- basin_wallacean %>%
  select(valid_name,time,event,year_description) %>%
  filter(time >= 0 & event == 1) %>%
  group_by(valid_name, year_description) %>%  
  summarise(
    mean_lag_time = mean(time, na.rm = TRUE),  
    geo_records_count = n(),  
    .groups = "drop"  
  )

loess_model <- loess(mean_lag_time ~ year_description, data = taxon_geo_lag_clean)
y <- taxon_geo_lag_clean$mean_lag_time
y_hat <- predict(loess_model)  
total_var <- var(y, na.rm = TRUE)  
resid_var <- var(y - y_hat, na.rm = TRUE)  
loess_r2 <- round(1 - (resid_var / total_var), 3)
anno_text <- paste0(
  "n = ", format(loess_model$n, big.mark = ","), "\n",
  "Loess R² = ", loess_r2
)


p1 <- ggplot(taxon_geo_lag_clean, 
             aes(x = year_description, 
                 y = mean_lag_time)) +  
  geom_point(alpha = 0.3, size = 0.8, color = "#2166AC")+
  geom_smooth(method = "loess",
              se = TRUE,
              color = "#B2182B",   
              fill = "#F4A582", 
              linetype = "solid",
              linewidth = 1,        
              span = 0.75) +        
  annotate("text", x = 1770, y = Inf, label = "(A)", size = 5, fontface = "bold",hjust = 0, vjust = 2) +
  annotate("text",
           x = max(taxon_geo_lag_clean$year_description) - 100,  
           y = 250,                                
           label = anno_text,
           hjust = 0, vjust = 1,
           size = 3.5,              
           color = "black") +
  labs(
    x = "Taxonomic description year",
    y = "Mean lag time (years) to first georeferenced occurrence"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(colour = "black",size = 10),
    axis.text = element_text(colour = "black",size = 9),
  ) +
  scale_x_continuous(breaks = seq(1760, 2020, 50),expand = c(0.01,0.01)) +
  scale_y_continuous(breaks = seq(0, max(taxon_geo_lag_clean$mean_lag_time) + 10, 50),expand = c(0.01,0.01))



basin_darwinian <- readRDS("output/stan_survdata_basin_darwinian_single.rds")

taxon_seq_lag_clean <- basin_darwinian %>%
  select(valid_name,time,event,year_description) %>%
  filter(time >= 0 & event == 1) %>%
  group_by(valid_name, year_description) %>%  
  summarise(
    mean_lag_time = mean(time, na.rm = TRUE),  
    geo_records_count = n(),  
    .groups = "drop"  
  )

lm_model <- lm(mean_lag_time ~ year_description, data = taxon_seq_lag_clean)
lm_sum <- summary(lm_model)

n_sample_lm <- nrow(taxon_seq_lag_clean)
lm_r2 <- round(lm_sum$r.squared, 3)



p2 <- ggplot(taxon_seq_lag_clean, 
             aes(x = year_description, 
                 y = mean_lag_time)) +  
  geom_point(alpha = 0.2, color = "#006837", size = 0.8) +
  geom_smooth(method = "lm", 
              se = TRUE, 
              color = "#B2182B", 
              fill = "#F4A582",
              linewidth = 1) +
  annotate("text", x = 1770, y = Inf, label = "(B)", size = 5, fontface = "bold",hjust = 0, vjust = 2) +
  annotate("text",
           x = max(taxon_seq_lag_clean$year_description) - 100,  
           y = 250,                                
           label = paste0(
             "n = ", format(n_sample_lm, big.mark = ","), "\n",
             "Linear R² = ", lm_r2
           ),
           hjust = 0, vjust = 1,
           size = 3.5,              
           color = "black") +
  labs(
    x = "Taxonomic description year",
    y = "Mean lag time (years) to first sequenced record"
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(colour = "black",size = 10),
    axis.text = element_text(colour = "black",size = 9),
  ) +
  scale_x_continuous(breaks = seq(1760, 2020, 50),expand = c(0.01,0.01)) +
  scale_y_continuous(breaks = seq(0, max(taxon_geo_lag_clean$mean_lag_time) + 10, 50),expand = c(0.01,0.01))



p1 | p2
ggsave(
  filename = "figures/supplement/Figure_S17.png",
  width    = 20,
  height   = 12,
  units    = "cm",
  dpi      = 300
)


# =========================
# 0) Packages (no duplicates)
# =========================
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(sf)
library(rnaturalearth)
library(hexbin)
library(gganimate)
library(transformr)   # gganimate needs it for sf tweening
library(gifski)

# =========================
# 1) Read & clean data
# =========================
rr <- readRDS("input/processed/basin_sp_geo.rds")

rr2 <- rr %>%
  mutate(
    year_description = suppressWarnings(as.integer(year_description)),
    year_geolocation = suppressWarnings(as.integer(year_geolocation)),
    decimalLongitude = suppressWarnings(as.numeric(decimalLongitude)),
    decimalLatitude  = suppressWarnings(as.numeric(decimalLatitude))
  ) %>%
  filter(
    !is.na(valid_name),
    !is.na(year_description),
    !is.na(year_geolocation),
    !is.na(decimalLongitude),
    !is.na(decimalLatitude)
  )

dat <- rr2
if (nrow(dat) == 0) stop("No data after filtering.")

# =========================
# 2) Common time axis
#    Use geolocation years as the animation clock
# =========================
t_min <- min(dat$year_geolocation, na.rm = TRUE)
t_max <- max(dat$year_geolocation, na.rm = TRUE)
time_seq <- seq(t_min, t_max, by = 1L)

# =========================
# 3) Base map & projection
# =========================
world <- ne_countries(scale = "medium", returnclass = "sf") %>%
  filter(name != "Antarctica")

crs_world <- "+proj=robin"
world_proj <- st_transform(world, crs_world)

# =========================
# 4) Map points (first georeferenced occurrences)
# =========================
map_df <- dat %>%
  transmute(
    valid_name,
    year_description,
    year_geolocation,
    lon = decimalLongitude,
    lat = decimalLatitude
  )

# cumulative expansion for map (points persist after appearance year)
map_cum <- map_df %>%
  tidyr::expand_grid(frame_year = time_seq) %>%
  filter(frame_year >= year_geolocation) %>%
  mutate(frame_year = as.integer(frame_year))

map_cum_sf <- st_as_sf(map_cum, coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs_world) %>%
  cbind(st_coordinates(.)) %>%
  as.data.frame()

p_map_cum <- ggplot() +
  geom_sf(
    data = world_proj,
    fill = "#121212",
    colour = "#2E2E2E",
    linewidth = 0.2
  ) +
  stat_binhex(
    data = map_cum_sf,
    aes(x = X, y = Y, fill = after_stat(count)),
    bins = 200,
    alpha = 0.9
  ) +
  scale_fill_viridis_c(option = "C", trans = "sqrt", name = "Density") +
  coord_sf(expand = FALSE) +
  labs(
    title = "First georeferenced occurrences (cumulative)",
    subtitle = "Hex-binned density | Year: {frame_time}",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.background  = element_rect(fill = "#0E0E0E", colour = NA),
    panel.background = element_rect(fill = "#0E0E0E", colour = NA),
    panel.grid       = element_blank(),
    plot.title       = element_text(hjust = 0.5, colour = "grey90", face = "bold"),
    plot.subtitle    = element_text(hjust = 0.5, colour = "grey75"),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    legend.position  = c(0.1, 0.05),
    legend.justification = c(0, 0),
    legend.key.size  = unit(0.5, "cm"),
    legend.title     = element_text(colour = "grey85", size = 10),
    legend.text      = element_text(colour = "grey80", size = 9),
    legend.background = element_rect(fill = "#0E0E0E", colour = NA),
    legend.key       = element_rect(fill = "#0E0E0E", colour = NA)
  ) +
  transition_time(frame_year) +
  ease_aes("linear")

# =========================
# 5) Latitudinal cumulative pattern
# =========================
lat_step <- 0.5

lat_year_first <- map_df %>%
  distinct(valid_name, year_geolocation, lat) %>%
  mutate(lat_bin = floor(lat / lat_step) * lat_step) %>%
  count(year_geolocation, lat_bin, name = "new_spp") %>%
  rename(frame_year = year_geolocation) %>%
  mutate(frame_year = as.integer(frame_year))

lat_bins_all <- sort(unique(floor(map_df$lat / lat_step) * lat_step))

lat_new_full <- expand_grid(frame_year = time_seq, lat_bin = lat_bins_all) %>%
  left_join(lat_year_first, by = c("frame_year", "lat_bin")) %>%
  mutate(new_spp = replace_na(new_spp, 0L)) %>%
  arrange(lat_bin, frame_year) %>%
  group_by(lat_bin) %>%
  mutate(
    cum_spp  = cumsum(new_spp),
    prev_cum = cum_spp - new_spp
  ) %>%
  ungroup()

p_lat_stack <- ggplot(lat_new_full) +
  geom_col(
    aes(x = lat_bin, y = prev_cum, fill = prev_cum),
    alpha = 0.25,
    width = lat_step
  ) +
  geom_col(
    aes(x = lat_bin, y = new_spp, fill = new_spp),
    width = lat_step
  ) +
  scale_fill_viridis_c(option = "C", trans = "sqrt", guide = "none") +
  labs(
    title = "Latitudinal pattern (cumulative)",
    subtitle = paste0("Latitude bins = ", lat_step, " | Year: {frame_time}"),
    x = "Latitude (binned)",
    y = "First georeferenced occurrence (cumulative)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.background  = element_rect(fill = "#0E0E0E", colour = NA),
    panel.background = element_rect(fill = "#0E0E0E", colour = NA),
    panel.grid       = element_blank(),
    plot.title       = element_text(hjust = 0.5, colour = "grey90", face = "bold"),
    plot.subtitle    = element_text(hjust = 0.5, colour = "grey75"),
    axis.title.x     = element_text(colour = "grey80"),
    axis.title.y     = element_text(colour = "grey80"),
    axis.text        = element_text(colour = "grey70"),
    axis.ticks       = element_line(colour = "grey50")
  ) +
  transition_time(frame_year) +
  view_follow(fixed_x = TRUE, fixed_y = FALSE) +
  ease_aes("linear")

# =========================
# 6) Lag cumulative density (desc -> first georef)
# =========================
lag_df <- map_df %>%
  distinct(valid_name, year_description, year_geolocation) %>%
  mutate(lag = year_geolocation - year_description) %>%
  filter(is.finite(lag), lag >= 0)

lag_cum <- lag_df %>%
  expand_grid(frame_year = time_seq) %>%
  filter(frame_year >= year_geolocation) %>%
  mutate(frame_year = as.integer(frame_year))

p_lag_dark <- ggplot() +
  stat_binhex(
    data = lag_cum,
    aes(x = year_description, y = lag, fill = after_stat(count)),
    bins = 200,
    alpha = 0.9
  ) +
  scale_fill_viridis_c(option = "C", trans = "sqrt", name = "Density") +
  labs(
    title = "Lag to first georeferenced occurrence (cumulative)",
    subtitle = "Hex-binned density | Year: {frame_time}",
    x = "Taxonomic description year",
    y = "Lag time (years) to first georeferenced occurrence"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.background  = element_rect(fill = "#0E0E0E", colour = NA),
    panel.background = element_rect(fill = "#0E0E0E", colour = NA),
    panel.grid       = element_blank(),
    plot.title       = element_text(hjust = 0.5, colour = "grey90", face = "bold"),
    plot.subtitle    = element_text(hjust = 0.5, colour = "grey75"),
    axis.title.x     = element_text(colour = "grey80"),
    axis.title.y     = element_text(colour = "grey80"),
    axis.text        = element_text(colour = "grey70"),
    axis.ticks       = element_line(colour = "grey50"),
    legend.position  = c(0.8, 0.7),
    legend.justification = c(0, 0),
    legend.key.size  = unit(0.5, "cm"),
    legend.title     = element_text(colour = "grey85", size = 10),
    legend.text      = element_text(colour = "grey80", size = 9),
    legend.background = element_rect(fill = "#0E0E0E", colour = NA),
    legend.key       = element_rect(fill = "#0E0E0E", colour = NA)
  ) +
  transition_time(frame_year) +
  view_follow(fixed_x = TRUE, fixed_y = FALSE) +
  ease_aes("linear")

# =========================
# 7) Render & save (only once each)
# =========================
anim_map <- animate(
  p_map_cum,
  nframes = length(time_seq),
  fps = 12,
  width = 1200,
  height = 560,
  renderer = gifski_renderer()
)

anim_lat <- animate(
  p_lat_stack,
  nframes = length(time_seq),
  fps = 12,
  width = 600,
  height = 560,
  renderer = gifski_renderer()
)

anim_lag <- animate(
  p_lag_dark,
  nframes = length(time_seq),
  fps = 12,
  width = 600,
  height = 560,
  renderer = gifski_renderer()
)

anim_save("figures/supplement/Figure_S18_map.gif", animation = anim_map)
anim_save("figures/supplement/Figure_S18_lat.gif", animation = anim_lat)
anim_save("figures/supplement/Figure_S18_lag.gif", animation = anim_lag)


p1 <- ggplot() +
  geom_sf(
    data = world_proj,
    #fill = "#121212",
    #colour = "#2E2E2E",
    linewidth = 0.2
  ) +
  stat_binhex(
    data = map_cum_sf %>% filter(frame_year == 2024),
    aes(x = X, y = Y, fill = after_stat(count)),
    bins = 200,
    alpha = 0.9
  ) +
  scale_fill_viridis_c(option = "C", trans = "sqrt", name = "Density",breaks = c(0,100,200,300)) +
  coord_sf(expand = FALSE) +
  labs(
    title = "First georeferenced occurrences (cumulative)",
    subtitle = "Hex-binned density | Year: 2024",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    #plot.background  = element_rect(fill = "#0E0E0E", colour = NA),
    #panel.background = element_rect(fill = "#0E0E0E", colour = NA),
    panel.grid       = element_blank(),
    plot.title       = element_text(hjust = 0.5,  face = "bold", size = 8),
    plot.subtitle    = element_text(hjust = 0.5,  size = 7),
    # plot.title       = element_text(hjust = 0.5, colour = "grey90", face = "bold", size = 8),
    # plot.subtitle    = element_text(hjust = 0.5, colour = "grey75", size = 7),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    legend.position  = c(0.1, 0.05),
    legend.justification = c(0, 0),
    legend.key.size  = unit(0.4, "cm"),
    # legend.title     = element_text(colour = "grey85", size = 8),
    # legend.text      = element_text(colour = "grey80", size = 8),
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 8),
    #legend.background = element_rect(fill = "#0E0E0E", colour = NA),
    #legend.key       = element_rect(fill = "#0E0E0E", colour = NA)
  )


p2 <- ggplot(lat_new_full %>% filter(frame_year == 2024),) +
  geom_col(
    aes(x = lat_bin, y = prev_cum, fill = prev_cum),
    alpha = 0.5,
    width = lat_step
  ) +
  geom_col(
    aes(x = lat_bin, y = new_spp, fill = new_spp),
    width = lat_step
  ) +
  scale_fill_viridis_c(option = "C", trans = "sqrt", guide = "none") +
  labs(
    title = "Latitudinal pattern (cumulative)",
    subtitle = paste0("Latitude bins = ", lat_step, " | Year: 2024"),
    x = "Latitude (binned)",
    y = "First georeferenced occurrence (cumulative)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    # plot.background  = element_rect(fill = "#0E0E0E", colour = NA),
    # panel.background = element_rect(fill = "#0E0E0E", colour = NA),
    #panel.grid       = element_blank(),
    plot.title       = element_text(hjust = 0.5, colour = "black", face = "bold", size = 8),
    plot.subtitle    = element_text(hjust = 0.5, colour = "black", size = 7),
    # axis.title.x     = element_text(colour = "grey80", size = 8),
    # axis.title.y     = element_text(colour = "grey80", size = 8),
    # axis.text        = element_text(colour = "grey70"),
    # axis.ticks       = element_line(colour = "grey50")
    axis.title.x     = element_text(colour = "black", size = 8),
    axis.title.y     = element_text(colour = "black", size = 8),
    axis.text        = element_text(colour = "black"),
    axis.ticks       = element_line(colour = "black")
  )

p3 <- ggplot() +
  stat_binhex(
    data = lag_cum %>% filter(frame_year == 2024),
    aes(x = year_description, y = lag, fill = after_stat(count)),
    bins = 200,
    alpha = 0.9
  ) +
  scale_fill_viridis_c(option = "C", trans = "sqrt", name = "Density",breaks = c(0,20,40,60)) +
  labs(
    title = "Lag to first georeferenced occurrence (cumulative)",
    subtitle = "Hex-binned density | Year: 2024",
    x = "Taxonomic description year",
    y = "Lag time (years) to first georeferenced occurrence"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    # plot.background  = element_rect(fill = "#0E0E0E", colour = NA),
    # panel.background = element_rect(fill = "#0E0E0E", colour = NA),
    #panel.grid       = element_blank(),
    plot.title       = element_text(hjust = 0.5, colour = "black", face = "bold", size = 8),
    plot.subtitle    = element_text(hjust = 0.5, colour = "black", size = 7),
    axis.title.x     = element_text(colour = "black", size = 8),
    axis.title.y     = element_text(colour = "black", size = 8),
    axis.text        = element_text(colour = "black"),
    axis.ticks       = element_line(colour = "black"),
    legend.position  = c(0.8, 0.6),
    legend.justification = c(0, 0),
    legend.key.size  = unit(0.4, "cm"),
    # legend.title     = element_text(colour = "grey85", size = 8),
    # legend.text      = element_text(colour = "grey80", size = 8),
    legend.title     = element_text(size = 8),
    legend.text      = element_text(size = 8),
    # legend.background = element_rect(fill = "#0E0E0E", colour = NA),
    # legend.key       = element_rect(fill = "#0E0E0E", colour = NA)
  )

library(cowplot)
cowplot::ggdraw() +
  cowplot::draw_plot(p1, x = 0.00, y = 0.50, width = 1.00, height = 0.50) +
  cowplot::draw_plot(p2, x = 0.00, y = 0.00, width = 0.50, height = 0.50) +
  cowplot::draw_plot(p3, x = 0.50, y = 0.00, width = 0.50, height = 0.50) +
  cowplot::draw_plot_label(
    label = c("A", "B", "C"),
    x     = c(0.01, 0.01, 0.51),  
    y     = c(0.99, 0.49, 0.49),  
    hjust = 0, vjust = 1,
    size  = 14,
    colour = "black",
    fontface = "bold"
  ) 
ggsave("figures/supplement/Figure_S18.png",dpi = 300, units = "cm", width = 21, height = 18)