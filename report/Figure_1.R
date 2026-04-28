# ==============================================================================
# Fig.1 Bivariate shortfall maps (LS–WS, LS–DS, WS–DS) with in-plot cell-count legend
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(biscale)
library(classInt)
library(cowplot)
library(here)
library(ggrastr)
set.seed(123)

# ---- I/O ----
linnaean_sf <- read.csv(here("output", "tables", "basin_linnaean_shortfall.csv"))
wallacean_sf <- read.csv(here("output", "tables", "basin_wallacean_shortfall.csv"))
darwinian_sf <- read.csv(here("output", "tables", "basin_darwinian_shortfall.csv"))

basin_sf <- readRDS(here("input", "raw", "basin", "basin_sf_v1.rds"))

# ---- Helpers ----
mm01 <- function(x, eps = 1e-4) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (diff(rng) + eps)
}

log_mm01 <- function(x) mm01(log10(x + 1))

spearman_label <- function(x, y) {
  ct <- suppressWarnings(
    cor.test(x, y, method = "spearman", use = "complete.obs", exact = FALSE)
  )
  rho <- unname(ct$estimate)
  p <- ct$p.value
  if (is.na(p)) {
    return("")
  }

  rho_txt <- sprintf("%.2f", rho)

  if (p < 0.001) {
    sprintf('r * "=" * %s * ", " * italic(p) * "<" * 0.001', rho_txt)
  } else {
    sprintf('r * "=" * %s * ", " * italic(p) * "=" * %s', rho_txt, signif(p, 2))
  }
}


make_legend_counts <- function(df_pts, grid_df, x_lab, y_lab) {
  label_txt <- spearman_label(df_pts$x_raw, df_pts$y_raw)

  counts <- df_pts %>%
    count(bi_class, name = "n_cells") %>%
    right_join(grid_df %>% select(bi_class, x, y), by = "bi_class") %>%
    replace_na(list(n_cells = 0))

  ggplot() +
    geom_tile(data = grid_df, aes(x = x, y = y, fill = bi_fill), colour = "white") +
    scale_fill_identity() +
    geom_text(
      data = counts,
      aes(x = x, y = y, label = n_cells),
      colour = "white",
      fontface = "bold",
      size = 1.8
    ) +
    scale_x_continuous(breaks = 1:4, labels = c("Low", "", "", "High"), expand = c(0.015, 0.015)) +
    scale_y_continuous(breaks = 1:4, labels = c("Low", "", "", "High"), expand = c(0.015, 0.015)) +
    labs(x = x_lab, y = y_lab) +
    coord_equal(clip = "off") +
    annotate(
      "text",
      x = Inf, y = 4.8,
      hjust = 1.2, vjust = 0,
      size = 2.2,
      label = label_txt,
      parse = TRUE
    ) +
    theme_void(base_size = 8) +
    theme(
      axis.title = element_text(size = 8, colour = "black", face = "bold"),
      axis.text = element_text(size = 6, colour = "black"),
      axis.title.y.left = element_text(angle = 90),
      axis.text.y.left = element_text(angle = 90),
      plot.margin = margin(t = 20, r = 3, b = 3, l = 3),
      panel.ontop = TRUE
    )
}

make_bimap <- function(sf_dat, world_sf, y_limits, pal = "PurpleOr", dim = 4) {
  ggplot(sf_dat) +
    ggrastr::rasterise(geom_sf(data = world_sf, fill = "lightgrey", colour = NA), dpi = 300) +
    ggrastr::rasterise(
      geom_sf(aes(fill = bi_class), colour = "white", linewidth = 0.03, show.legend = FALSE),
      dpi = 300
    ) +
    bi_scale_fill(pal = pal, dim = dim) +
    bi_theme() +
    coord_sf(crs = "+proj=moll +lon_0=0", ylim = y_limits, expand = FALSE)
}

make_panel <- function(dt_sf, x, y, x_lab, y_lab, pal = "PurpleOr", dim = 4) {
  stopifnot(nrow(dt_sf) == length(x), nrow(dt_sf) == length(y))

  dt_xy <- dt_sf %>%
    mutate(
      x_raw = x,
      y_raw = y
    )

  # k-means breaks (computed on non-NA, then applied to full vectors)
  x_brks <- classIntervals(stats::na.omit(dt_xy$x_raw), n = dim, style = "kmeans")$brks
  y_brks <- classIntervals(stats::na.omit(dt_xy$y_raw), n = dim, style = "kmeans")$brks

  bivar_sf <- dt_xy %>%
    mutate(
      x_bin = cut(x_raw, breaks = x_brks, include.lowest = TRUE),
      y_bin = cut(y_raw, breaks = y_brks, include.lowest = TRUE)
    ) %>%
    select(basin_id, x_raw, y_raw, x_bin, y_bin) %>%
    drop_na(x_raw, y_raw, x_bin, y_bin) %>%
    bi_class(x = x_bin, y = y_bin, style = "quantile", dim = dim)

  pal_obj <- bi_pal(pal = pal, dim = dim)
  grid_df <- pal_obj$layers[[1]][["data"]]

  pts <- bivar_sf %>%
    st_drop_geometry() %>%
    left_join(grid_df %>% select(bi_class, x, y), by = "bi_class")

  legend_g <- make_legend_counts(pts, grid_df, x_lab = x_lab, y_lab = y_lab)
  map_g <- make_bimap(bivar_sf, world_sf = world_map, y_limits = y_limits, pal = pal, dim = dim)

  ggdraw() +
    draw_plot(map_g, 0, 0, 1, 1) +
    draw_plot(legend_g, -0.1, 0.1, 0.55, 0.55)
}


# ---- Build analysis table ----
shortfalls <- linnaean_sf %>%
  select(basin_id, ls = SRdesc) %>%
  left_join(wallacean_sf %>% select(basin_id, ws = SRnoloc), by = "basin_id") %>%
  left_join(darwinian_sf %>% select(basin_id, ds = SRnoseq), by = "basin_id") %>%
  mutate(
    ls = log_mm01(ls),
    ws = log_mm01(ws),
    ds = log_mm01(ds)
  )

dt_sf <- basin_sf %>%
  left_join(shortfalls, by = "basin_id") %>%
  st_wrap_dateline(options = "WRAPDATELINE=YES")

# ---- Basemap + projection bounds ----
world_map <- ne_countries(scale = 50, type = "countries", returnclass = "sf")

moll_crs <- st_crs("+proj=moll")
lat_sf <- data.frame(lon = 0, lat = c(-60, 90)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs = moll_crs)

y_limits <- st_coordinates(lat_sf)[, "Y"]

# ---- Panels ----
p_ls_ws <- make_panel(
  dt_sf,
  x = dt_sf$ls, y = dt_sf$ws,
  x_lab = expression(LS %->% ""), y_lab = expression(WS %->% "")
)

p_ls_ds <- make_panel(
  dt_sf,
  x = dt_sf$ls, y = dt_sf$ds,
  x_lab = expression(LS %->% ""), y_lab = expression(DS %->% "")
)

p_ws_ds <- make_panel(
  dt_sf,
  x = dt_sf$ws, y = dt_sf$ds,
  x_lab = expression(WS %->% ""), y_lab = expression(DS %->% "")
)

# ---- Stack + export ----
p_combined <- plot_grid(
  p_ls_ws + theme(plot.margin = margin(0, 0, 0, 0)),
  p_ls_ds + theme(plot.margin = margin(0, 0, 0, 0)),
  p_ws_ds + theme(plot.margin = margin(0, 0, 0, 0)),
  ncol = 1,
  align = "v",
  axis = "lr",
  labels = c("A", "B", "C"),
  label_size = 10,
  label_fontface = "bold",
  label_x = 0.10,
  label_y = 0.95
)

ggsave(
  filename = here("figures", "main", "Figure_1.png"),
  plot = p_combined,
  width = 15, height = 20, units = "cm",
  dpi = 300
)
