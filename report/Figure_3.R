# ==============================================================================
# Fig.3 Shortfall scenarios vs hotspots + NSCI map (here-based I/O, de-duplicated)
# ==============================================================================

# ---- Libraries ----
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(cowplot)
library(sf)
library(rnaturalearth)
library(biscale)
library(ggrastr)
library(viridis)
library(scales)
library(legendry)
library(here)

sf::sf_use_s2(FALSE)
options(warn = -1)

# ---- Paths ----
path_linnaean  <- here("output", "tables", "basin_linnaean_shortfall.csv")
path_wallacean <- here("output", "tables", "basin_wallacean_shortfall.csv")
path_darwinian <- here("output", "tables", "basin_darwinian_shortfall.csv")

path_basin_sf  <- here("input", "raw", "basin", "basin_sf_v1.rds")
path_biorealm  <- here("input", "raw", "biogeographic_list.csv")
path_hotspots  <- here("input", "raw", "hotspots_2016_1", "hotspots_2016_1_fix.shp")

path_out_tbl   <- here("output", "tables", "basin_shortfalls.csv")
path_out_rds   <- here("input", "processed", "basin_shortfall.rds")
path_out_fig   <- here("figures", "main", "Figure_3.png")

# ---- Helpers ----
mm01 <- function(x, eps = 1e-4) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (diff(rng) + eps)
}

log_mm01 <- function(x) mm01(log10(x + 1))

within_realm_state <- function(data, value_col, realm_col = "biogeographic_realm", q = 0.25) {
  stopifnot(value_col %in% names(data), realm_col %in% names(data))
  stopifnot(q > 0 && q < 0.5)
  
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

scenario_from_states <- function(ls_state, ws_state, ds_state) {
  case_when(
    ls_state == "high" & ws_state == "high" & ds_state == "high" ~ "LS↑ WS↑ DS↑",
    ls_state == "high" & ws_state == "high" & ds_state == "low"  ~ "LS↑ WS↑ DS↓",
    ls_state == "high" & ws_state == "low"  & ds_state == "high" ~ "LS↑ WS↓ DS↑",
    ls_state == "high" & ws_state == "low"  & ds_state == "low"  ~ "LS↑ WS↓ DS↓",
    ls_state == "low"  & ws_state == "high" & ds_state == "high" ~ "LS↓ WS↑ DS↑",
    ls_state == "low"  & ws_state == "high" & ds_state == "low"  ~ "LS↓ WS↑ DS↓",
    ls_state == "low"  & ws_state == "low"  & ds_state == "high" ~ "LS↓ WS↓ DS↑",
    TRUE ~ NA_character_
  )
}

label_expr_for_scenario <- function(s) {
  dplyr::recode(
    s,
    "LS↑ WS↑ DS↑" = 'LS*""%up%""*WS*""%up%""*DS*""%up%""',
    "LS↑ WS↑ DS↓" = 'LS*""%up%""*WS*""%up%""*DS*""%down%""',
    "LS↑ WS↓ DS↑" = 'LS*""%up%""*WS*""%down%""*DS*""%up%""',
    "LS↑ WS↓ DS↓" = 'LS*""%up%""*WS*""%down%""*DS*""%down%""',
    "LS↓ WS↑ DS↑" = 'LS*""%down%""*WS*""%up%""*DS*""%up%""',
    "LS↓ WS↑ DS↓" = 'LS*""%down%""*WS*""%up%""*DS*""%down%""',
    "LS↓ WS↓ DS↑" = 'LS*""%down%""*WS*""%down%""*DS*""%up%""',
    .default = NA_character_
  )
}

# ---- Read inputs ----
linnaean  <- read.csv(path_linnaean)
wallacean <- read.csv(path_wallacean)
darwinian <- read.csv(path_darwinian)

biorealm  <- read.csv(path_biorealm)
hotspots  <- st_read(path_hotspots, quiet = TRUE)
inland_sf <- readRDS(path_basin_sf)

# ---- Join + standardize (for scenario classification) ----
basin_vals <- linnaean %>%
  select(basin_id, LS = SRdesc) %>%
  left_join(wallacean %>% select(basin_id, WS = SRnoloc), by = "basin_id") %>%
  left_join(darwinian %>% select(basin_id, DS = SRnoseq), by = "basin_id")

df <- basin_vals %>%
  left_join(biorealm, by = "basin_id") %>%
  mutate(
    biogeographic_realm = factor(biogeographic_realm)
  ) %>%
  drop_na(LS, WS, DS, biogeographic_realm)

# ---- Sensitivity (q grid) ----
qs <- seq(0.01, 0.30, by = 0.01)

scenario_by_q <- map_dfr(qs, function(qi) {
  tmp <- df %>%
    mutate(
      LS_state = within_realm_state(., "LS", q = qi),
      WS_state = within_realm_state(., "WS", q = qi),
      DS_state = within_realm_state(., "DS", q = qi),
      scenario = scenario_from_states(LS_state, WS_state, DS_state),
      q = qi
    )
  
  tmp %>%
    count(q, scenario, name = "n") %>%
    filter(!is.na(scenario))
})

p_sens <- ggplot(scenario_by_q, aes(x = q, y = n, colour = scenario)) +
  geom_line() +
  labs(x = "Quantile threshold (q)", y = "Number of basins", colour = "Scenario") +
  theme_minimal()

# ---- Scenario table (q = 0.3) ----
df_q30 <- df %>%
  mutate(
    LS_state = within_realm_state(., "LS", q = 0.3),
    WS_state = within_realm_state(., "WS", q = 0.3),
    DS_state = within_realm_state(., "DS", q = 0.3),
    scenario = scenario_from_states(LS_state, WS_state, DS_state)
  )

# If you already have basin_shortfalls_fix.csv, keep using it.
# Otherwise you can write it from df_q30:
# dir.create(here("output", "tables"), showWarnings = FALSE, recursive = TRUE)
# write.csv(df_q30, path_out_tbl, row.names = FALSE)
basin_shortfalls <- read.csv(path_out_tbl)

# ---- Basemap bounds (Mollweide y-limits) ----
world_map <- ne_countries(scale = 50, type = "countries", returnclass = "sf")

moll_crs <- st_crs("+proj=moll")
lat_sf <- data.frame(lon = 0, lat = c(-60, 90)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs = moll_crs)
y_limits <- st_coordinates(lat_sf)[, "Y"]

# ---- Scenario sf + hotspot overlay ----
scenario_sf <- inland_sf %>%
  left_join(basin_shortfalls %>% select(basin_id, scenario), by = "basin_id") %>%
  st_wrap_dateline(options = "WRAPDATELINE=YES")

hotspots <- st_transform(hotspots, st_crs(scenario_sf))

scenario_sf <- scenario_sf %>%
  mutate(
    is_hotspot = lengths(st_intersects(geometry, hotspots %>% filter(Type != "outer limit"))) > 0,
    is_shortfall = !is.na(scenario)
  )

# ---- Bivariate class: Shortfall (x) × Hotspot (y) ----
scenario_sf <- scenario_sf %>%
  mutate(
    is_hotspot  = factor(is_hotspot,  levels = c(FALSE, TRUE)),
    is_shortfall= factor(is_shortfall,levels = c(FALSE, TRUE))
  ) %>%
  bi_class(x = is_shortfall, y = is_hotspot, dim = 2)


# ---- In-plot bivariate legend with proportions ----
custom_pal <- c(
  "1-1" = "#d3d3d3", # no shortfall, no hotspot
  "2-1" = "#4A3A3B", # shortfall, no hotspot
  "1-2" = "#F16C31", # no shortfall, hotspot
  "2-2" = "#7C0809"  # shortfall, hotspot
)

pal_obj  <- bi_pal(pal = "PurpleOr", dim = 2)
grid_df  <- pal_obj$layers[[1]][["data"]]
grid_df$bi_fill <- custom_pal

pts <- scenario_sf %>%
  st_drop_geometry() %>%
  left_join(grid_df %>% select(bi_class, x, y), by = "bi_class")

# ---- Panel A: scenario × hotspot counts (facetted) ----
scenario_levels <- c(
  "LS↑ WS↑ DS↑","LS↑ WS↑ DS↓","LS↑ WS↓ DS↑","LS↑ WS↓ DS↓",
  "LS↓ WS↑ DS↑","LS↓ WS↑ DS↓","LS↓ WS↓ DS↑"
)

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
  mutate(scenario = 'LS* "" %up% "" *WS* "" %up% "" *DS* "" %down% ""')

count_scenario_3 <- pts %>%
  select(is_hotspot,scenario) %>%
  mutate(scenario = ifelse(scenario == "LS↑ WS↓ DS↑","LS↑ WS↓ DS↑",NA),
         is_shortfall = ifelse(is.na(scenario),FALSE,TRUE)) %>%
  mutate(is_shortfall = factor(is_shortfall),
         is_hotspot = factor(is_hotspot)) %>%
  bi_class(x=is_shortfall,y= is_hotspot,dim = 2) %>%
  group_by(scenario,bi_class) %>%
  summarise(n_cells = n(),.groups = "drop") %>%
  mutate(scenario = 'LS* "" %up% "" *WS* "" %down% "" *DS* "" %up% ""')  

count_scenario_4 <- pts %>%
  select(is_hotspot,scenario) %>%
  mutate(scenario = ifelse(scenario == "LS↑ WS↓ DS↓","LS↑ WS↓ DS↓",NA),
         is_shortfall = ifelse(is.na(scenario),FALSE,TRUE)) %>%
  mutate(is_shortfall = factor(is_shortfall),
         is_hotspot = factor(is_hotspot)) %>%
  bi_class(x=is_shortfall,y= is_hotspot,dim = 2) %>%
  group_by(scenario,bi_class) %>%
  summarise(n_cells = n(),.groups = "drop") %>%
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

count_scenario_6 <- pts %>%
  select(is_hotspot,scenario) %>%
  mutate(scenario = ifelse(scenario == "LS↓ WS↑ DS↓","LS↓ WS↑ DS↓",NA),
         is_shortfall = ifelse(is.na(scenario),FALSE,TRUE)) %>%
  mutate(is_shortfall = factor(is_shortfall),
         is_hotspot = factor(is_hotspot)) %>%
  bi_class(x=is_shortfall,y= is_hotspot,dim = 2) %>%
  group_by(scenario,bi_class) %>%
  summarise(n_cells = n(),.groups = "drop") %>%
  rbind(data.frame(scenario = NA,bi_class = "2-1",n_cells = 0)) %>%
  mutate(scenario = 'LS* "" %down% "" *WS* "" %up% "" *DS* "" %down% ""')     

count_scenario_7 <- pts %>%
  select(is_hotspot,scenario) %>%
  mutate(scenario = ifelse(scenario == "LS↓ WS↓ DS↑","LS↓ WS↓ DS↑",NA),
         is_shortfall = ifelse(is.na(scenario),FALSE,TRUE)) %>%
  mutate(is_shortfall = factor(is_shortfall),
         is_hotspot = factor(is_hotspot)) %>%
  bi_class(x=is_shortfall,y= is_hotspot,dim = 2) %>%
  group_by(scenario,bi_class) %>%
  summarise(n_cells = n(),.groups = "drop")  %>%
  mutate(scenario = 'LS* "" %down% "" *WS* "" %down% "" *DS* "" %up% ""')      

count_scenario <- rbind(count_scenario_1,count_scenario_2,count_scenario_3,count_scenario_4,
                        count_scenario_5,count_scenario_6,count_scenario_7)
rm(count_scenario_1,count_scenario_2,count_scenario_3,count_scenario_4,
   count_scenario_5,count_scenario_6,count_scenario_7)

count_scenario <- count_scenario %>%
  left_join(grid_df,by = "bi_class")

count_scenario$scenario <- factor(count_scenario$scenario,
                                  levels = c('LS* "" %up% "" *WS* "" %up% "" *DS* "" %up% ""',
                                             'LS* "" %up% "" *WS* "" %up% "" *DS* "" %down% ""',
                                             'LS* "" %up% "" *WS* "" %down% "" *DS* "" %up% ""',
                                             'LS* "" %up% "" *WS* "" %down% "" *DS* "" %down% ""',
                                             'LS* "" %down% "" *WS* "" %up% "" *DS* "" %up% ""',
                                             'LS* "" %down% "" *WS* "" %up% "" *DS* "" %down% ""',
                                             'LS* "" %down% "" *WS* "" %down% "" *DS* "" %up% ""')
)

p_counts <- ggplot(count_scenario) +
  geom_tile(aes(x = x, y = y, fill = bi_fill), colour = "white") +
  scale_fill_identity() +
  geom_text(
    aes(x = x, y = y, label = format(n_cells, big.mark = ",")),
    colour = "white", fontface = "bold", size = 2
  ) +
  scale_x_continuous(breaks = 1:2, labels = c("N", "Y"), expand = c(0.015, 0.015)) +
  scale_y_continuous(breaks = 1:2, labels = c("N", "Y"), expand = c(0.015, 0.015)) +
  labs(x = "Combinations of biodiversity shortfalls", y = "Hotspot") +
  coord_equal() +
  facet_wrap(~ scenario, nrow = 1, labeller = label_parsed) +
  theme_void() +
  theme(
    axis.title  = element_text(size = 8, colour = "black"),
    axis.text   = element_text(size = 6, colour = "black"),
    axis.title.y.left = element_text(angle = 90),
    axis.text.y.left  = element_text(angle = 90),
    strip.text = element_text(size = 6),
    panel.spacing = unit(0, "lines"),
    plot.margin = margin(0, 0, 0, 0)
  )

# ---- Panel B: Map: shortfall × hotspot ----

count_df <- pts %>%
  count(bi_class, name = "n_cells") %>%
  right_join(grid_df %>% select(bi_class, x, y), by = "bi_class") %>%
  replace_na(list(n_cells = 0)) %>%
  mutate(pct = n_cells / sum(n_cells))

p_bi_legend <- ggplot() +
  geom_tile(data = grid_df, aes(x = x, y = y, fill = bi_fill), colour = "white") +
  scale_fill_identity() +
  geom_text(
    data = count_df,
    aes(x = x, y = y, label = scales::percent(pct, accuracy = 1)),
    colour = "white", fontface = "bold", size = 2.2
  ) +
  scale_x_continuous(breaks = 1:2, labels = c("N", "Y"), expand = c(0.015, 0.015)) +
  scale_y_continuous(breaks = 1:2, labels = c("N", "Y"), expand = c(0.015, 0.015)) +
  labs(x = "Shortfall", y = "Hotspot") +
  coord_equal(clip = "off") +
  theme_void(base_size = 8) +
  theme(
    axis.title = element_text(size = 8, colour = "black", face = "bold"),
    axis.text  = element_text(size = 6, colour = "black"),
    axis.title.y.left = element_text(angle = 90),
    axis.text.y.left  = element_text(angle = 90),
    plot.margin = margin(0, 0, 0, 0)
  )

p_bimap <- ggplot(scenario_sf) +
  ggrastr::rasterise(geom_sf(data = world_map, fill = "lightgrey", colour = NA), dpi = 300) +
  ggrastr::rasterise(
    geom_sf(aes(fill = bi_class), colour = "white", linewidth = 0.03, show.legend = FALSE),
    dpi = 300
  ) +
  bi_scale_fill(pal = custom_pal, dim = 2) +
  bi_theme() +
  coord_sf(crs = "+proj=moll +lon_0=0", ylim = y_limits, expand = FALSE) +
  theme(plot.margin = margin(0, 0, 0, 0))

p_map_with_legend <- ggdraw() +
  draw_plot(p_bimap, 0, 0, 1, 1) +
  draw_plot(p_bi_legend, 0.09, 0.16, 0.34, 0.34)


# ---- Panel C: Map: NSCI ----
df_nsc <- basin_vals %>%
  mutate(
    SL = log_mm01(LS),
    SW = log_mm01(WS),
    SD = log_mm01(DS),
    SW_adj = (1 - SL) * SW,
    SD_adj = (1 - SL) * SD,
    NSCI   = SL + SW_adj + SD_adj,
    NSCI_norm = mm01(NSCI)
  ) %>%
  select(basin_id, NSCI_norm)

#dir.create(dirname(path_out_rds), showWarnings = FALSE, recursive = TRUE)
#saveRDS(df_nsc, path_out_rds)

dt_nsc <- inland_sf %>%
  left_join(df_nsc, by = "basin_id") %>%
  st_wrap_dateline(options = "WRAPDATELINE=YES")

p_nsc <- ggplot(dt_nsc) +
  ggrastr::rasterise(geom_sf(data = world_map, fill = "lightgrey", colour = NA), dpi = 300) +
  ggrastr::rasterise(geom_sf(aes(fill = NSCI_norm), colour = "white", linewidth = 0.03), dpi = 300) +
  scale_fill_viridis_c(
    name = "NSCI",
    na.value = "lightgrey",
    option = "C",
    alpha = 1,
    begin = 0.1,
    end = 0.99,
    direction = -1,
    breaks = seq(0, 1, 0.5),
    labels = scales::number_format(accuracy = 0.1),
    guide = legendry::compose_sandwich(
      middle = gizmo_histogram(hist.args = list(breaks = 8),just = 1),
      #middle = legendry::gizmo_density(just = 1), #https://teunbrand.github.io/legendry/articles/guide_composition.html
      text   = "axis_base"
    )
  ) +
  bi_theme() +
  coord_sf(crs = "+proj=moll +lon_0=0", ylim = y_limits, expand = FALSE) +
  theme(
    legend.position = c(0.06, 0.13),
    legend.justification = c(0, 0),
    legend.background = element_blank(),
    legend.key.height = unit(0.5, "cm"),
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8, margin = margin(l = 4)),
    legend.axis.line = element_line(linewidth = 0.2, colour = "black")
  )

# ---- Assemble + export ----
p_final <- ggdraw() +
  draw_plot(p_counts, x = 0.18, y = 2/3, width = 0.6, height = 0.13) +
  draw_label("A", x = 0.17, y = 2/3 + 0.14, fontface = "bold", size = 12, hjust = 0, vjust = 1) +
  draw_plot(p_map_with_legend, x = 0, y = 1/3, width = 1.0, height = 0.33) +
  draw_label("B", x = 0.17, y = 1/3 + 0.33 - 0.02, fontface = "bold", size = 12, hjust = 0, vjust = 1) +
  draw_plot(p_nsc, x = 0, y = 0, width = 1.0, height = 0.33) +
  draw_label("C", x = 0.17, y = 0.31, fontface = "bold", size = 12, hjust = 0, vjust = 1)

dir.create(here("figures", "main"), showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = path_out_fig,
  plot = p_final,
  width = 20, height = 20, units = "cm",
  dpi = 300
)

