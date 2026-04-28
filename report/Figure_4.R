# ============================================================
# Fig.4 Prioritization of freshwater fish biodiversity collections based on multi-objective trade-offs
# ============================================================

rm(list = ls())

library(dplyr)
library(ggplot2)
library(GGally)
library(rPref)
library(patchwork)
library(legendry)
library(ggh4x)
library(cowplot)
library(here)

# extra libs used later in your code
library(emoa)     # Pareto dominance (minimization) emoa::nds_rank()
library(ggrepel)  # geom_text_repel
library(sf)
library(rnaturalearth)
library(ggrastr)

sf_use_s2(FALSE)

# -----------------------------
# Helpers
# -----------------------------
read_rds_h  <- function(...) readRDS(here::here(...))
read_csv_h  <- function(...) read.csv(here::here(...), stringsAsFactors = FALSE)

min_max_normalize_safe <- function(x, epsilon = 1e-6) {
  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)
  (x - min_val) / (max_val - min_val + epsilon)
}

# Realm-specific Pareto frontier helper
add_realm_front <- function(df, obj_cols) {
  df %>%
    group_by(biogeographic_realm) %>%
    group_modify(~ {
      obj_mat <- as.matrix(.x[, obj_cols, drop = FALSE])
      fronts  <- emoa::nds_rank(t(obj_mat))
      stopifnot(length(fronts) == nrow(.x))
      .x %>%
        mutate(
          front        = fronts,
          front_rank   = 1 - percent_rank(front),              # 1 = best (front=1), 0 = worst
          front_scaled = rank(-front, ties.method = "min")     # larger = better
        )
    }) %>%
    ungroup()
}

# tag best/worst basins by q-quantiles of front_scaled
tag_strategy <- function(df, q = 0.10) {
  df %>%
    group_by(biogeographic_realm) %>%
    mutate(
      strategy = case_when(
        front_scaled >= quantile(front_scaled, 1 - q, na.rm = TRUE) ~ "best",
        front_scaled <= quantile(front_scaled, q,     na.rm = TRUE) ~ "worst",
        TRUE ~ "not_significant"
      )
    ) %>%
    ungroup()
}


# -----------------------------
# Data I/O (here)
# -----------------------------
cost      <- read_rds_h("input", "processed", "cost_basin_all.rds")
sdg       <- read_rds_h("input", "processed", "sdg_basin.rds")
shortfall <- read_rds_h("input", "processed", "basin_shortfall.rds")

biogeographic_list <- read_csv_h("input", "raw", "biogeographic_list.csv") %>%
  mutate(
    biogeographic_realm = recode(
      biogeographic_realm,
      "Afrotropic"   = "AFR",
      "Australasia"  = "AUS",
      "Indomalayan"  = "IND",
      "Nearctic"     = "NEA",
      "Neotropic"    = "NEO",
      "Palearctic"   = "PAL",
      "Oceania"      = "OCE"
    )
  )

# -----------------------------
# PCA: SDG → Socioeconomics / Governance
# -----------------------------
sdg.pca <- prcomp(sdg[, 2:18], center = TRUE, scale. = TRUE)
sdg.pca$rotation

# PC1(Socioeconomics)# ensure higher value = stronger socioeconomic capacity
# Governance = PC2  # already correct direction
sdg$PC1 <- -sdg.pca$x[, 1]
sdg$PC2 <-  sdg.pca$x[, 2]
names(sdg)[19:20] <- c("Socioeconomics", "Governance")

# -----------------------------
# Merge to analysis table
# -----------------------------
data <- shortfall %>%
  select(basin_id = basin_id, Shortfall = NSCI_norm) %>%
  left_join(cost %>% select(basin_id = basin_id, Investment = cost_scale), by = "basin_id") %>%
  left_join(sdg  %>% select(basin_id, Socioeconomics, Governance), by = "basin_id") %>%
  left_join(biogeographic_list[, c(1, 3)], by = "basin_id")

basin_shortfalls <- read_csv_h("output", "tables", "basin_shortfalls.csv") %>%
  select(basin_id, scenario) %>%
  na.omit() %>%
  pull(basin_id)

# data <- data %>% filter(basin_id %in% basin_shortfalls)
data <- data %>% filter(basin_id %in% basin_shortfalls)

##################################################################################
df_nsc <- shortfall %>% left_join(biogeographic_list[,c(1,3)],by = "basin_id")

basins <- df_nsc %>% select(basin_id,NSCI_norm,biogeographic_realm) %>% na.omit()
# Step 1: Realm-based expectation model
model <- lm(NSCI_norm ~ biogeographic_realm, data = basins)

# Step 2: Realm-adjusted residuals
basins$resid_NSCI <- residuals(model)

# Step 3: Compute μ and σ of residuals
mu_resid <- mean(basins$resid_NSCI, na.rm = TRUE)
sd_resid <- sd(basins$resid_NSCI, na.rm = TRUE)

# Step 4: Cutoff values consistent with the published method
z <- 1.96
upper_cut <- mu_resid + z * sd_resid
lower_cut <- mu_resid - z * sd_resid

# Step 5: Classification
basins$NSCI_high   <- basins$resid_NSCI >= upper_cut
basins$NSCI_low    <- basins$resid_NSCI <= lower_cut
basins$NSCI_normal <- basins$resid_NSCI > lower_cut & basins$resid_NSCI < upper_cut
shortfall_8 <- basins %>% filter(NSCI_high == TRUE) %>% pull(basin_id)


data <- data %>% filter(basin_id %in% unique(basin_shortfalls, shortfall_8))

min_max_normalize_safe <- function(x, epsilon = 1e-6) {
  min_val <- min(x,na.rm = T)
  max_val <- max(x,na.rm = T)
  range_val <- max_val - min_val + epsilon
  return((x - min_val) / range_val)
}

data$Investment <- min_max_normalize_safe(data$Investment)
data$Shortfall <- min_max_normalize_safe(data$Shortfall)
data$Socioeconomics <- min_max_normalize_safe(data$Socioeconomics)
data$Governance <- min_max_normalize_safe(data$Governance)
rm(cost,sdg,sdg.pca,shortfall,biogeographic_list,basin_shortfalls,min_max_normalize_safe)
#-------------------------------------------------------------------------------
# Maps (use here for inputs, keep your styling)
library(sf)
library(rnaturalearth)
library(ggplot2)
sf_use_s2(FALSE)

inland <- read_rds_h("input", "raw", "basin", "basin_sf_v1.rds")
world_map <- rnaturalearth::ne_countries(scale = 50, type = "countries", returnclass = "sf")

moll_proj <- st_crs("+proj=moll")
lat_points <- data.frame(lon = 0, lat = c(-60, 90)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
projected_points <- st_transform(lat_points, crs = moll_proj)
y_limits <- st_coordinates(projected_points)[, "Y"]


################################################################################
# Biological-urgency baseline
df_plot_realm_1 <- data

# Linear model
model <- lm(Shortfall ~ biogeographic_realm, data = df_plot_realm_1)

# Realm-adjusted residuals
df_plot_realm_1$residuals <- residuals(model)

# Summary stats
res_mean <- mean(df_plot_realm_1$residuals, na.rm = TRUE)
res_sd   <- sd(df_plot_realm_1$residuals, na.rm = TRUE)

critical_value_upper <- res_mean + 1.96 * res_sd
critical_value_lower <- res_mean - 1.96 * res_sd

# Color gradient functions
get_green_shade  <- colorRampPalette(c("#BFBFBF", "#00B050"))
get_yellow_shade <- colorRampPalette(c("#FFC000", "#BFBFBF"))

# Function to map residuals to colors
residual_to_fill <- function(x, min_res, max_res) {
  sapply(x, function(v) {
    if (is.na(v)) return(NA_character_)
    
    if (v < 0) {
      # yellow -> grey
      idx <- floor(100 * (v - min_res) / (0 - min_res)) + 1
      idx <- max(1, min(100, idx))
      get_yellow_shade(100)[idx]
    } else {
      # grey -> green
      idx <- floor(100 * (v - 0) / (max_res - 0)) + 1
      idx <- max(1, min(100, idx))
      get_green_shade(100)[idx]
    }
  })
}

min_res <- min(df_plot_realm_1$residuals, na.rm = TRUE)
max_res <- max(df_plot_realm_1$residuals, na.rm = TRUE)

# Add fronts and fill to df_plot_realm_1
df_plot_realm_1 <- df_plot_realm_1 %>%
  mutate(
    fronts = case_when(
      residuals > critical_value_upper ~ "best",
      residuals < critical_value_lower ~ "worst",
      TRUE ~ "not_significant"
    ),
    fill = residual_to_fill(residuals, min_res, max_res)
  )

# Histogram breaks
breaks <- seq(min_res, max_res, length.out = 16)

# Histogram data
hist_data <- hist(df_plot_realm_1$residuals, breaks = breaks, plot = FALSE)

hist_df <- data.frame(
  xmin = head(hist_data$breaks, -1),
  xmax = tail(hist_data$breaks, -1),
  xmid = hist_data$mids,
  count = hist_data$counts,
  density = hist_data$density
)

# Map fill colors to histogram bins
hist_df$fill <- residual_to_fill(hist_df$xmid, min_res, max_res)

# Normal curve
x_vals <- seq(min_res, max_res, length.out = 400)
curve_df <- data.frame(
  x = x_vals,
  y = dnorm(x_vals, mean = res_mean, sd = res_sd) * nrow(df_plot_realm_1) * diff(breaks)[1]
)

# Label positions
ymax_hist <- max(hist_df$count, na.rm = TRUE)
label_y_top <- ymax_hist * 1.08
label_y_mid <- ymax_hist * 0.62

p1 <- ggplot() +
  geom_rect(data = hist_df,aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = count),
            fill = hist_df$fill,color = "grey40",linewidth = 0.2) +
  # normal curve scaled to counts
  geom_line(data = curve_df,aes(x = x, y = y),color = "black",linewidth = 0.2) +
  # critical lines
  geom_vline(xintercept = critical_value_lower,color = "#FFC000",linetype = 2,linewidth = 0.3) +
  geom_vline(xintercept = critical_value_upper,color = "#00B050",linetype = 2,linewidth = 0.3) +
  labs(
    x = expression(Residuals[NSCI]),
    y = "Frequency"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  scale_x_continuous(breaks = c(-0.5,0,0.5))+
  theme_classic() +
  theme(
    axis.title = element_text(size = 6),
    axis.text = element_text(colour = "black",size = 5),
    axis.line = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.ticks.length = unit(0.05, "cm"),
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  )

p1


case1 <- inland %>% left_join(df_plot_realm_1[,c(1,7,9)], by = "basin_id")
pp1 <- ggplot()+
  ggrastr::rasterise(geom_sf(data = world_map %>%
                               st_wrap_dateline(options = c("WRAPDATELINE=YES","DATELINEOFFSET=180")),
                             fill = "#d3d3d3", colour = NA),dpi = 300) +
  ggrastr::rasterise(geom_sf(data =case1 ,aes(fill = fill),colour = "white", linewidth = 0.03),dpi = 300)+
  scale_fill_identity(na.value = "#d3d3d3")+
  theme_void()+
  coord_sf(crs = "+proj=moll +lon_0=0", ylim = y_limits, expand = FALSE)+
  theme(legend.position = "none",
        plot.margin = margin(0,0,0,0))

pp1 <- ggdraw() +
  draw_plot(pp1, 0, 0, 1, 1) +                             
  draw_plot(p1, 0.06, 0.02, 0.2, 0.5) 


################################################################################
# Prepare data: Shortfall + Investment
df2 <- data %>%
  select(basin_id,biogeographic_realm, Shortfall, Investment) %>%
  filter(!is.na(Investment), !is.na(Shortfall)) %>%
  mutate(
    # Convert objectives to minimization form for nds_rank():
    # High Shortfall (good) → -Shortfall (so smaller = better)
    Shortfall_obj  = -Shortfall,
    Investment_obj = Investment # low Investment (good) → keep as original value
  )


df_plot_realm_2 <- add_realm_front(df2, c("Shortfall_obj", "Investment_obj")) %>%
  tag_strategy(q = 0.10)

p2 <- ggplot() +
  geom_point(data = df_plot_realm_2,aes(x = Investment,y = Shortfall,colour = front_rank),size  = 0.1,alpha = 0.5) +
  scale_colour_gradientn(name = "Pareto frontier",colours = c("#FFC000", "#BFBFBF", "#00B050"),
                         limits = c(0, 1),breaks = c(0, 1), labels = c("Worst","Best"),
                         guide = guide_colbar(show = c(TRUE, TRUE), oob = "squish"))+
  labs(x = "Cost", y = expression(Shortfall[NSCI])) +
  theme_classic() +
  scale_x_continuous(breaks  = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0,0.5, 1)) +
  theme(
    axis.text  = element_text(colour = "black", size = 5),
    axis.title = element_text(colour = "black", size = 6),
    legend.title = element_text(size = 6, hjust = 0.5),
    legend.text  = element_text(size = 5),
    axis.line = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.ticks.length = unit(0.05, "cm"),
    legend.key.width = unit(0.3,"cm"),
    legend.key.height = unit(0.15,"cm"),
    legend.background = element_blank(),
    legend.position = "none",
    legend.title.position = "top",
    legend.direction = "horizontal",
    plot.margin = margin(0,0,0,0)
  )

p_legend <- ggplot() +
  geom_point(data = df_plot_realm_2,aes(x = Investment,y = Shortfall,colour = front_rank),size  = 0.1,alpha = 0.5) +
  scale_colour_gradientn(name = "Pareto frontier",colours = c("#FFC000", "#BFBFBF", "#00B050"),
                         limits = c(0, 1),breaks = c(0, 1), labels = c("Worst","Best"),
                         guide = guide_colbar(show = c(TRUE, TRUE), oob = "squish"))+
  labs(x = "Cost", y = expression(Shortfall[NSCI])) +
  theme_void() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    plot.background = element_blank(),
    panel.grid = element_blank(),
    legend.title = element_text(size = 5, hjust = 0.5),
    legend.text  = element_text(size = 4),
    legend.key.width = unit(0.2,"cm"),
    legend.key.height = unit(0.1,"cm"),
    legend.background = element_blank(),
    legend.position = "top",  # 图例居中
    legend.title.position = "top",
    legend.direction = "horizontal",
    plot.margin = margin(0,0,0,0)
  )

case2 <- inland %>% left_join(df_plot_realm_2[,c(2,8)], by = "basin_id")
pp2 <- ggplot()+
  ggrastr::rasterise(geom_sf(data = world_map %>%
                               st_wrap_dateline(options = c("WRAPDATELINE=YES","DATELINEOFFSET=180")),
                             fill = "#d3d3d3", colour = NA),dpi = 300) +
  ggrastr::rasterise(geom_sf(data =case2 ,aes(fill = front_rank),colour = "white", linewidth = 0.03),dpi = 300)+
  scale_fill_gradientn(name = "Pareto frontier",
                       colours = c("#FFC000", "grey90","#00B050"),
                       limits = c(0, 1),                 
                       breaks = c(0, 1), 
                       na.value = "grey70",
                       labels = c("Worst","Best"),
                       guide = guide_colbar(
                         show = c(TRUE, TRUE), 
                         oob = "squish")
  )+
  theme_void()+
  coord_sf(crs = "+proj=moll +lon_0=0", ylim = y_limits, expand = FALSE)+
  theme(legend.position = "none",
        plot.margin = margin(0,0,0,0))

pp2 <- ggdraw() +
  draw_plot(pp2, 0, 0, 1, 1) +                             
  draw_plot(p2, 0.06, 0.02, 0.2, 0.5) +
  draw_plot(get_legend(p_legend), 0.07, 0.63, 0.06, 0.03) +
  draw_label(label = expression(Delta[NSCI] == -0.24),
             x = 0.58, y = 0.08,hjust = 1, vjust = 1, size = 7.2)+
  draw_label(label = expression("E-test, " * italic(p) * " < 0.01"),
             x = 0.8, y = 0.08,hjust = 1, vjust = 1, size = 7.2)+
  draw_label(label = expression(BWR[NSCI] * ' = 5.25'),
             x = 0.3, y = 0.98,hjust = 1, vjust = 1, size = 7.5)


################################################################################
#    Goal: prioritize basins with
#      - high Shortfall   (large knowledge gaps)
#      - high Socioeconomics  (strong capacity to fill those gaps)
# ------------------------------------------------------------
df3 <- data %>%
  select(basin_id, biogeographic_realm,Shortfall, Socioeconomics) %>%       
  filter(!is.na(Shortfall), !is.na(Socioeconomics)) %>%
  mutate(
    # Both objectives need to be maximized,
    # so convert them into minimization for nds_rank():
    Shortfall_obj = -Shortfall, # higher shortfall = better priority
    Socio_obj = -Socioeconomics # higher socioeconomic capacity = better
  )

df_plot_realm_3 <- add_realm_front(df3, c("Shortfall_obj", "Socio_obj")) %>%
  tag_strategy(q = 0.10)

p3 <- ggplot() +
  geom_point(data = df_plot_realm_3,aes(x = Socioeconomics,y = Shortfall,colour = front_rank),size  = 0.1,alpha = 0.5) +
  scale_colour_gradientn(name = "Pareto frontier",colours = c("#FFC000", "#BFBFBF", "#00B050"),
                         limits = c(0, 1),breaks = c(0, 1), labels = c("Worst","Best"),
                         guide = guide_colbar(show = c(TRUE, TRUE), oob = "squish"))+
  #coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
  labs(x = "Income", y = expression(Shortfall[NSCI])) +
  theme_classic() +
  scale_x_continuous(breaks  = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0,0.5, 1)) +
  theme(
    axis.text  = element_text(colour = "black", size = 5),
    axis.title = element_text(colour = "black", size = 6),
    legend.title = element_text(size = 6, hjust = 0.5),
    legend.text  = element_text(size = 5),
    axis.line = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.ticks.length = unit(0.05, "cm"),
    legend.key.width = unit(0.3,"cm"),
    legend.key.height = unit(0.15,"cm"),
    legend.background = element_blank(),
    legend.position = "none",
    legend.title.position = "top",
    legend.direction = "horizontal",
    plot.margin = margin(0,0,0,0)
  )


case3 <- inland %>% left_join(df_plot_realm_3[,c(2,8)], by = "basin_id")
pp3 <- ggplot()+
  ggrastr::rasterise(geom_sf(data = world_map %>%
                               st_wrap_dateline(options = c("WRAPDATELINE=YES","DATELINEOFFSET=180")),
                             fill = "#d3d3d3", colour = NA),dpi = 300)+
  ggrastr::rasterise(geom_sf(data =case3 ,aes(fill = front_rank),colour = "white", linewidth = 0.03),dpi = 300)+
  scale_fill_gradientn(name = "Pareto frontier",
                       colours = c("#FFC000", "grey90","#00B050"),
                       limits = c(0, 1),                 
                       breaks = c(0, 1), 
                       na.value = "grey70",
                       labels = c("Worst","Best"),
                       guide = guide_colbar(
                         show = c(TRUE, TRUE), 
                         oob = "squish")
  )+
  theme_void()+
  coord_sf(crs = "+proj=moll +lon_0=0", ylim = y_limits, expand = FALSE)+
  theme(legend.position = "none",
        plot.margin = margin(0,0,0,0))

pp3 <- ggdraw() +
  draw_plot(pp3, 0, 0, 1, 1) +                             
  draw_plot(p3, 0.06, 0.02, 0.2, 0.5) +
  draw_label(label = expression(Delta[NSCI] == -0.34),
             x = 0.58, y = 0.08,hjust = 1, vjust = 1, size = 7.2)+
  draw_label(label = expression("E-test, " * italic(p) * " < 0.01"),
             x = 0.8, y = 0.08,hjust = 1, vjust = 1, size = 7.2)+
  draw_label(label = expression(BWR[NSCI] * ' = 3.05'),
             x = 0.3, y = 0.98,hjust = 1, vjust = 1, size = 7.5)

################################################################################
#    Goal: prioritize basins with
#      - high Shortfall   (large knowledge gaps)
#      - high Governance  (strong capacity to fill those gaps)
# ------------------------------------------------------------
df4 <- data %>%
  select(basin_id, biogeographic_realm,Governance, Shortfall) %>%
  filter(!is.na(Governance), !is.na(Shortfall)) %>%
  mutate(
    # Convert both objectives to "minimization" form for nds_rank():
    # We want to MAXIMIZE Shortfall → use -Shortfall
    Shortfall_obj  = -Shortfall,
    # We want to MAXIMIZE Governance → use -Governance
    Governance_obj = -Governance
  )

df_plot_realm_4 <- add_realm_front(df4, c("Shortfall_obj", "Governance_obj")) %>%
  tag_strategy(q = 0.10)

p4 <- ggplot() +
  geom_point(data = df_plot_realm_4,aes(x = Governance,y = Shortfall,colour = front_rank),size  = 0.1,alpha = 0.5) +
  scale_colour_gradientn(name = "Pareto frontier",colours = c("#FFC000", "#BFBFBF", "#00B050"),
                         limits = c(0, 1),breaks = c(0, 1), labels = c("Worst","Best"),
                         guide = guide_colbar(show = c(TRUE, TRUE), oob = "squish"))+
  #coord_fixed(xlim = c(0,1), ylim = c(0,1)) +
  labs(x = "Protection", y = expression(Shortfall[NSCI])) +
  theme_classic() +
  scale_x_continuous(breaks  = c(0, 0.5, 1)) +
  scale_y_continuous(breaks = c(0,0.5, 1)) +
  theme(
    axis.text  = element_text(colour = "black", size = 5),
    axis.title = element_text(colour = "black", size = 6),
    legend.title = element_text(size = 6, hjust = 0.5),
    legend.text  = element_text(size = 5),
    axis.line = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.ticks.length = unit(0.05, "cm"),
    legend.key.width = unit(0.3,"cm"),
    legend.key.height = unit(0.15,"cm"),
    legend.background = element_blank(),
    legend.position = "none",
    legend.title.position = "top",
    legend.direction = "horizontal",
    plot.margin = margin(0,0,0,0)
  )

case4 <- inland %>% left_join(df_plot_realm_4[,c(2,8)], by = "basin_id")
pp4 <- ggplot()+
  ggrastr::rasterise(geom_sf(data = world_map %>%
                               st_wrap_dateline(options = c("WRAPDATELINE=YES","DATELINEOFFSET=180")),
                             fill = "#d3d3d3", colour = NA),dpi = 300)+
  ggrastr::rasterise(geom_sf(data =case4 ,aes(fill = front_rank),colour = "white", linewidth = 0.03),dpi = 300)+
  scale_fill_gradientn(name = "Pareto frontier",
                       colours = c("#FFC000", "grey90","#00B050"),
                       limits = c(0, 1),                 
                       breaks = c(0, 1), 
                       na.value = "grey70",
                       labels = c("Worst","Best"),
                       guide = guide_colbar(
                         show = c(TRUE, TRUE), 
                         oob = "squish")
  )+
  theme_void()+
  coord_sf(crs = "+proj=moll +lon_0=0", ylim = y_limits, expand = FALSE)+
  theme(legend.position = "none",
        plot.margin = margin(0,0,0,0))

pp4 <- ggdraw() +
  draw_plot(pp4, 0, 0, 1, 1) +                             
  draw_plot(p4, 0.06, 0.02, 0.2, 0.5)  +
  draw_label(label = expression(Delta[NSCI] == -0.43),
             x = 0.58, y = 0.08,hjust = 1, vjust = 1, size = 7.2)+
  draw_label(label = expression("E-test, " * italic(p) * " < 0.01"),
             x = 0.8, y = 0.08,hjust = 1, vjust = 1, size = 7.2)+
  draw_label(label = expression(BWR[NSCI] * ' = 8.59'),
             x = 0.3, y = 0.98,hjust = 1, vjust = 1, size = 7.5)
################################################################################
# ------------------------------------------------------------
#  Prepare data for 4-objective trade-off:
#    Goal: High Shortfall × Low Investment × High Socioeconomics × High Governance
#
#    Assumptions:
#    - Shortfall:       higher = larger knowledge gaps (priority ↑)
#    - Investment:      higher = higher total cost (priority ↓)
#    - Socioeconomics:  higher = stronger socioeconomic capacity
#    - Governance:      higher = stronger institutional/governance capacity
# ------------------------------------------------------------
df5 <- data %>%
  select(basin_id, biogeographic_realm, Investment, Socioeconomics, Governance, Shortfall) %>%
  na.omit() %>%
  mutate(
    # Convert all objectives to "minimize" form for nds_rank():
    # We want to MAXIMIZE Shortfall → minimize -Shortfall
    Shortfall_obj      = -Shortfall,
    # We want to MINIMIZE Investment → keep as is
    Investment_obj     =  Investment,
    # We want to MAXIMIZE Socioeconomics → minimize -Socioeconomics
    Socioeconomics_obj = -Socioeconomics,
    # We want to MAXIMIZE Governance → minimize -Governance
    Governance_obj     = -Governance
  )


df_plot_realm_5 <- add_realm_front(df5, c("Shortfall_obj", "Investment_obj","Socioeconomics_obj","Governance_obj")) %>%
  tag_strategy(q = 0.10)


p5 <- ggparcoord(df_plot_realm_5,columns = c(7,8,9,10), groupColumn = 12,
                 scale = "uniminmax", 
                 showPoints = F,
                 alphaLines = 0.1)+ 
  scale_colour_gradientn(colours = c("#FFC000", "grey90","#00B050"))+
  #geom_point(size = 1)+
  geom_path(linewidth = 0.01,alpha = 0.5)+
  #annotate(geom = "text",x = 3.6,y = 1.05,label = "BWR = 34", size = 2)+
  theme_minimal()+
  #coord_fixed(ratio = 3) +
  scale_x_discrete(expand = c(0,0),labels = c("NSCI","Cost","Income","Protection"))+
  #scale_x_discrete(expand = c(0,0),labels = c("S","C","I","P"))+
  xlab("")+
  ylab("")+
  scale_y_continuous(breaks = c(0,0.5, 1)) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey50", linewidth = 0.3),
        legend.position = "none",
        axis.text.y = element_text(colour = "black",size = 5),
        axis.text.x = element_text(colour = "black",size = 5, angle = 45, vjust = 1.1, hjust = 0.8),
        axis.line.x  = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05, "cm"),
        plot.margin = margin(0,0,0,0)
  )

case5 <- inland %>% left_join(df_plot_realm_5[,c(2,12)], by = "basin_id")

pp5 <- ggplot()+
  ggrastr::rasterise(geom_sf(data = world_map %>%
                               st_wrap_dateline(options = c("WRAPDATELINE=YES","DATELINEOFFSET=180")),
                             fill = "#d3d3d3", colour = NA),dpi = 300)+
  ggrastr::rasterise(geom_sf(data =case5 ,aes(fill = front_rank),colour = "white", linewidth = 0.03),dpi = 300)+
  scale_fill_gradientn(name = "Pareto frontier",
                       colours = c("#FFC000", "grey90","#00B050"),
                       limits = c(0, 1),                 
                       breaks = c(0, 1), 
                       na.value = "grey70",
                       labels = c("Worst","Best"),
                       guide = guide_colbar(
                         show = c(TRUE, TRUE), 
                         oob = "squish")
  )+
  theme_void()+
  coord_sf(crs = "+proj=moll +lon_0=0", ylim = y_limits, expand = FALSE)+
  theme(legend.position = "none",
        plot.margin = margin(0,0,0,0))
pp5 <- ggdraw() +
  draw_plot(pp5, 0, 0, 1, 1) +                             
  draw_plot(p5, 0.01, 0.01, 0.2, 0.55)  +
  draw_label(label = expression(Delta[NSCI] == -0.28),
             x = 0.58, y = 0.08,hjust = 1, vjust = 1, size = 7.2)+
  draw_label(label = expression("E-test, " * italic(p) * " < 0.01"),
             x = 0.8, y = 0.08,hjust = 1, vjust = 1, size = 7.2)+
  draw_label(label = expression(BWR[NSCI] * ' = 2.18'),
             x = 0.3, y = 0.98,hjust = 1, vjust = 1, size = 7.5)


# ------------------------------------------------------------
#  Prepare data for 6-objective trade-off:
#    Goal: High Shortfall × Low Investment × High Socioeconomics × High Governance
#
#    Assumptions:
#    - Shortfall:       higher = larger knowledge gaps (priority ↑)
#    - Investment-adj:      higher = higher total cost (priority ↓)
#    - Socioeconomics:  higher = stronger socioeconomic capacity
#    - Governance:      higher = stronger institutional/governance capacity
# ------------------------------------------------------------
df6 <- data %>%
  select(basin_id, biogeographic_realm, Investment, Socioeconomics, Governance, Shortfall) %>%
  mutate(Investment_adj = Investment* (1-0.5*(1-Socioeconomics))) %>%
  na.omit() %>%
  mutate(
    # Convert all objectives to "minimize" form for nds_rank():
    # We want to MAXIMIZE Shortfall → minimize -Shortfall
    Shortfall_obj      = -Shortfall,
    # We want to MINIMIZE Investment → keep as is
    Investment_obj     =  Investment_adj,
    # We want to MAXIMIZE Governance → minimize -Governance
    Governance_obj     = -Governance
  )

df_plot_realm_6 <- add_realm_front(df6, c("Shortfall_obj", "Investment_obj","Governance_obj")) %>%
  tag_strategy(q = 0.10)


p6 <- ggparcoord(df_plot_realm_6,columns = c(8,9,10), groupColumn = 12,
                 scale = "uniminmax", 
                 showPoints = F,
                 alphaLines = 0.1)+ 
  scale_colour_gradientn(colours = c("#FFC000", "grey90","#00B050"))+
  #geom_point(size = 1)+
  geom_path(linewidth = 0.01,alpha = 0.5)+
  #annotate(geom = "text",x = 3.6,y = 1.05,label = "BWR = 34", size = 2)+
  theme_minimal()+
  #coord_fixed(ratio = 3) +
  scale_x_discrete(expand = c(0,0),labels = c("NSCI",expression(Cost[adj]),"Protection"))+
  xlab("")+
  ylab("")+
  scale_y_continuous(breaks = c(0,0.5, 1)) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey50", linewidth = 0.3),
        legend.position = "none",
        axis.text.y = element_text(colour = "black",size = 5),
        axis.text.x = element_text(colour = "black",size = 5, angle = 45, vjust = 1.1, hjust = 0.8),
        axis.line.x  = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2),
        axis.ticks.length = unit(0.05, "cm"),
        plot.margin = margin(0,0,0,0)
  )

case6 <- inland %>% left_join(df_plot_realm_6[,c(2,12)], by = "basin_id")
pp6 <- ggplot()+
  ggrastr::rasterise(geom_sf(data = world_map %>%
                               st_wrap_dateline(options = c("WRAPDATELINE=YES","DATELINEOFFSET=180")),
                             fill = "#d3d3d3", colour = NA),dpi = 300)+
  ggrastr::rasterise(geom_sf(data =case6 ,aes(fill = front_rank),colour = "white", linewidth = 0.03),dpi = 300)+
  scale_fill_gradientn(name = "Pareto frontier",
                       colours = c("#FFC000", "grey90","#00B050"),
                       limits = c(0, 1),                 
                       breaks = c(0, 1), 
                       na.value = "grey70",
                       labels = c("Worst","Best"),
                       guide = guide_colbar(
                         show = c(TRUE, TRUE), 
                         oob = "squish")
  )+
  theme_void()+
  coord_sf(crs = "+proj=moll +lon_0=0", ylim = y_limits, expand = FALSE)+
  theme(legend.position = "none",
        plot.margin = margin(0,0,0,0))
pp6 <- ggdraw() +
  draw_plot(pp6, 0, 0, 1, 1) +                             
  draw_plot(p6, 0.01, 0.01, 0.2, 0.55)  +
  draw_label(label = expression(Delta[NSCI] == -0.25),
             x = 0.58, y = 0.08,hjust = 1, vjust = 1, size = 7.2)+
  draw_label(label = expression("E-test, " * italic(p) * " < 0.01"),
             x = 0.8, y = 0.08,hjust = 1, vjust = 1, size = 7.2)+
  draw_label(label = expression(BWR[NSCI] * ' = 3.80'),
             x = 0.3, y = 0.98,hjust = 1, vjust = 1, size = 7.5)

row_pair <- function(left_plot, right_plot, labels) {
  cowplot::plot_grid(
    left_plot  + theme(plot.margin = margin(5, 5, 0, 0, unit = "pt")),
    right_plot + theme(plot.margin = margin(5, 0, 0, 0, unit = "pt")),
    ncol = 2,
    rel_widths = c(0.5, 0.5),
    labels = labels,
    label_size = 9,
    label_fontface = "bold",
    align = "h",
    axis  = "tb"
  )
}

row1 <- row_pair(pp1, pp2, labels = c("A","B"))
row2 <- row_pair(pp3, pp4, labels = c("C","D"))
row3 <- row_pair(pp5, pp6, labels = c("E","F"))

fig <- row1 / row2 / row3
fig

# Save with here()
ggsave(
  filename = here::here("figures", "main", "Figure_4.png"),
  plot = fig,
  dpi = 300,
  units = "cm",
  width = 18,
  height = 15
)


################################################################################
get_denominator_vector <- function(df_subset, scenario) {
  valid_scenarios <- c("scenario_2", "scenario_3", "scenario_4", "scenario_5","scenario_6")
  if (!scenario %in% valid_scenarios) {
    stop(paste("Unsupported scenario:", scenario, "| Valid scenarios:", paste(valid_scenarios, collapse = ", ")))
  }
  if (scenario == "scenario_2") {
    return(df_subset$Investment)
  } else if (scenario == "scenario_3") {
    return(df_subset$Socioeconomics)
  } else if (scenario == "scenario_4") {
    return(df_subset$Governance)
  }else if (scenario == "scenario_5") {
    return(df_subset$Investment * df_subset$Socioeconomics * df_subset$Governance)
  } else if (scenario == "scenario_6") {
    return(df_subset$Investment * df_subset$Governance)
  }
}
calculate_BWR_with_CI <- function(df, scenario) {
  best  <- df %>% filter(strategy == "best")
  worst <- df %>% filter(strategy == "worst")
  if (nrow(best) == 0 | nrow(worst) == 0) {
    result <- tibble(
      scenario = scenario,         
      n_best = 0, 
      n_worst = 0,
      BWR = NA_real_,               
      BWR_CI_95_lower = NA_real_,
      BWR_CI_95_upper = NA_real_
    )
    cat(sprintf("Scenario %s: No valid best/worst group (n_best=%d, n_worst=%d)\n", 
                scenario, result$n_best, result$n_worst))
    return(result)
  }
  
  eps <- 1e-6
  denom_best <- get_denominator_vector(best, scenario)
  E_best_raw <- best$Shortfall / denom_best
  E_best     <- E_best_raw[is.finite(E_best_raw) & !is.na(E_best_raw) & denom_best > eps]
  denom_worst <- get_denominator_vector(worst, scenario)
  E_worst_raw <- worst$Shortfall / denom_worst
  E_worst     <- E_worst_raw[is.finite(E_worst_raw) & !is.na(E_worst_raw) & denom_worst > eps]
  if (length(E_best) == 0 | length(E_worst) == 0) {
    result <- tibble(
      scenario = scenario,
      n_best = length(E_best),
      n_worst = length(E_worst),
      BWR = NA_real_,
      BWR_CI_95_lower = NA_real_,
      BWR_CI_95_upper = NA_real_
    )
    cat(sprintf("Scenario %s: No valid E_best/E_worst (n_best=%d, n_worst=%d)\n", 
                scenario, result$n_best, result$n_worst))
    return(result)
  }
  mean_best  <- mean(E_best, na.rm = TRUE)
  mean_worst <- mean(E_worst, na.rm = TRUE)
  BWR        <- mean_best / mean_worst
  
  se_best  <- sd(E_best, na.rm = TRUE) / sqrt(length(E_best))
  se_worst <- sd(E_worst, na.rm = TRUE) / sqrt(length(E_worst))
  if (mean_best == 0 | mean_worst == 0) {
    result <- tibble(
      scenario = scenario,
      n_best = length(E_best),
      n_worst = length(E_worst),
      BWR = 0,
      BWR_CI_95_lower = 0,
      BWR_CI_95_upper = 0
    )
    cat(sprintf("Scenario %s: BWR=0 (mean_best=%.4f, mean_worst=%.4f)\n", 
                scenario, mean_best, mean_worst))
    return(result)
  }
  se_log_BWR <- sqrt((se_best / mean_best)^2 + (se_worst / mean_worst)^2)
  df_total   <- length(E_best) + length(E_worst) - 2
  t_crit     <- if (df_total > 0) qt(0.975, df = df_total) else 1.96 
  
  log_BWR     <- log(BWR)
  CI_lower    <- exp(log_BWR - t_crit * se_log_BWR)
  CI_upper    <- exp(log_BWR + t_crit * se_log_BWR)
  result <- tibble(
    scenario = scenario,
    n_best = length(E_best),          
    n_worst = length(E_worst),
    BWR = round(BWR, 3),              
    BWR_CI_95_lower = round(CI_lower, 4),  
    BWR_CI_95_upper = round(CI_upper, 3)   
  )
  cat(sprintf(
    "BWR = %.3f (95%% CI: %.3f–%.3f), n_best = %d, n_worst = %d\n",
    result$BWR, result$BWR_CI_95_lower, result$BWR_CI_95_upper,
    result$n_best, result$n_worst
  ))
  return(result)
}

# ---------------------- scenario_2 ---------------------- #
calculate_BWR_with_CI(df_plot_realm_2 %>% mutate(scenario = "scenario_2") , scenario = "scenario_2")
# BWR = 5.252 (95% CI: 4.244–6.499), n_best = 225, n_worst = 176

# ---------------------- scenario_3 ---------------------- #
calculate_BWR_with_CI(df_plot_realm_3 %>% mutate(scenario = "scenario_3") , scenario = "scenario_3")
# BWR = 3.052 (95% CI: 1.609–5.790), n_best = 203, n_worst = 346

# ---------------------- scenario_4 ---------------------- #
calculate_BWR_with_CI(df_plot_realm_4 %>% mutate(scenario = "scenario_4") , scenario = "scenario_4")
# BWR = 8.587 (95% CI: 5.835–12.635), n_best = 183, n_worst = 426

# ---------------------- scenario_5 ---------------------- #
calculate_BWR_with_CI(df_plot_realm_5 %>% mutate(scenario = "scenario_5") , scenario = "scenario_5")
# BWR = 2.180 (95% CI: 1.500–3.171), n_best = 322, n_worst = 174

# ---------------------- scenario_6 ---------------------- #
calculate_BWR_with_CI(df_plot_realm_6 %>% mutate(scenario = "scenario_6") , scenario = "scenario_6")
# BWR = 3.804 (95% CI: 2.866–5.050), n_best = 225, n_worst = 189

