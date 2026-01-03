# ------------------------------------------------------------------------------
# Supplementary Figure S1
# National-level bivariate patterns of biodiversity knowledge shortfalls in freshwater fishes

library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(tidyr)
library(Hmisc)
library(sf)
library(rnaturalearth)
library(ggplot2)
library(biscale)
library(classInt)
library(stringr)
library(ggtext)


linnaean <- read.csv("output/tables/country_linnaean_shortfall.csv")
wallacean <- read.csv("output/tables/country_wallacean_shortfall.csv")
darwinian <- read.csv("output/tables/country_darwinian_shortfall.csv")

df <- linnaean[,c(1,10)] %>% 
  left_join(wallacean[,c(1,6)], by = "iso3") %>%
  left_join(darwinian[,c(1,6)], by = "iso3")

min_max_normalize_safe <- function(x, epsilon = 1e-4) {
  min_val <- min(x,na.rm = T)
  max_val <- max(x,na.rm = T)
  range_val <- max_val - min_val + epsilon
  return((x - min_val) / range_val)
}

df$SRdesc <- log10(df$SRdesc+1)
df$SRnoloc <- log10(df$SRnoloc+1)
df$SRnoseq <- log10(df$SRnoseq+1)

df$SRdesc <- min_max_normalize_safe(df$SRdesc)
df$SRnoloc <- min_max_normalize_safe(df$SRnoloc)
df$SRnoseq <- min_max_normalize_safe(df$SRnoseq)



inland <- readRDS("input/raw/country_20251212.rds")
dt <- inland %>% left_join(df, by = "iso3") %>%
  st_wrap_dateline(options = c("WRAPDATELINE=YES"))
moll_proj <- st_crs("+proj=moll")
lat_points <- data.frame(lon = 0, lat = c(-60, 90)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
projected_points <- st_transform(lat_points, crs = moll_proj)
y_limits <- st_coordinates(projected_points)[, "Y"]


leg <- bi_pal(pal = "PurpleOr", dim = 4)
grid_dat <- leg$layers[[1]][["data"]]

dt$SRdesc_bin <- cut(dt$SRdesc, breaks = classIntervals(dt$SRdesc, n = 4, style = "kmeans")$brks, include.lowest = TRUE)
dt$SRnoloc_bin <- cut(dt$SRnoloc, breaks = classIntervals(dt$SRnoloc, n = 4, style = "kmeans")$brks, include.lowest = TRUE)

d1 <- dt %>% dplyr::select(iso3,SRdesc,SRnoloc,SRdesc_bin,SRnoloc_bin) %>%
  tidyr::drop_na() %>% bi_class(x= SRdesc_bin,y = SRnoloc_bin, style = "quantile",dim = 4)


d1_pts <- d1 %>% st_drop_geometry() %>%
  left_join(
    leg$layers[[1]][["data"]] %>% select(bi_class, x, y),
    by = "bi_class"
  )

cor_res <- cor.test(
  d1_pts$SRdesc,
  d1_pts$SRnoloc,
  method = "spearman",
  use    = "complete.obs"
)

rho   <- unname(cor_res$estimate)
pval  <- cor_res$p.value
label_expr <- if (pval < 0.001) {
  bquote(
    "r = " ~ .(sprintf("%.2f", rho)) * ", " *
      italic(p) < 0.001
  )
} else {
  bquote(
    "r = " ~ .(sprintf("%.2f", rho)) * ", " *
      italic(p) == .(signif(pval, 2))
  )
}

count_df_1 <- d1_pts %>%
  count(bi_class, name = "n_cells") %>%
  right_join(
    grid_dat %>% select(bi_class, x, y),
    by = "bi_class"
  ) %>%
  replace_na(list(n_cells = 0))


bi_legend_scatter_1 <- ggplot() +
  geom_tile(data = grid_dat,aes(x = x, y = y, fill = bi_fill),colour = "white") +
  scale_fill_identity() +
  geom_text(data = count_df_1,aes(x = x, y = y, label = n_cells),colour   = "white",fontface = "bold", size = 1.8)+
  scale_x_continuous(breaks = 1:4,labels = c("Low", "", "", "High"),expand = c(0.015, 0.015)) +
  scale_y_continuous(breaks = 1:4,labels = c("Low", "", "", "High"),expand = c(0.015, 0.015)) +
  labs(
    # x = expression(S[Linnaean] %->% ""), 
    # y = expression(S[Wallacean] %->% "")
    x = expression(LS %->% ""), 
    y = expression(WS %->% "")
  ) +
  coord_equal(clip = "off") +
  annotate("text",x = Inf, y = 4.8, hjust = 1.2, vjust = 0,size  = 2.2,label = deparse(label_expr),parse = TRUE)+
  theme_void(base_size = 8) +
  theme(
    axis.title  = element_text(size = 8, colour = "black", face = "bold"),
    axis.text   = element_text(size = 6, colour = "black"),
    axis.title.y.left = element_text(angle = 90),
    axis.text.y.left = element_text(angle = 90),
    plot.margin = margin(t = 20, r = 3, b = 3, l = 3),
    panel.ontop = TRUE
  )

bi_legend_scatter_1



map1 <- ggplot(data = d1) +
  ggrastr::rasterise(geom_sf(aes(fill = bi_class),colour = "white", linewidth = 0.03, show.legend = FALSE),dpi = 300)+
  bi_scale_fill(pal = "PurpleOr",dim = 4)+
  bi_theme()+
  coord_sf(crs = "+proj=moll +lon_0=0", ylim = y_limits, expand = FALSE)

pp1 <- ggdraw() +
  draw_plot(map1, 0, 0, 1, 1) +
  draw_plot(bi_legend_scatter_1, -0.1, 0.1, 0.55, 0.55)


################################################################################

dt$SRnoseq_bin <- cut(dt$SRnoseq, breaks = classIntervals(dt$SRnoseq, n = 4, style = "kmeans")$brks, include.lowest = TRUE)

d2 <- dt %>% select(iso3,SRdesc,SRnoseq,SRdesc_bin,SRnoseq_bin) %>%
  na.omit() %>% bi_class(x= SRdesc_bin,y = SRnoseq_bin, style = "quantile",dim = 4)

d2_pts <- d2 %>% st_drop_geometry() %>%
  left_join(
    leg$layers[[1]][["data"]] %>% select(bi_class, x, y),
    by = "bi_class"
  )

cor_res <- cor.test(
  d2_pts$SRdesc,
  d2_pts$SRnoseq,
  method = "spearman",
  use    = "complete.obs"
)

rho   <- unname(cor_res$estimate)
pval  <- cor_res$p.value

label_expr <- if (pval < 0.001) {
  bquote(
    "r = " ~ .(sprintf("%.2f", rho)) * ", " *
      italic(p) < 0.001
  )
} else {
  bquote(
    "r = " ~ .(sprintf("%.2f", rho)) * ", " *
      italic(p) == .(signif(pval, 2))
  )
}

count_df_2 <- d2_pts %>%
  count(bi_class, name = "n_cells") %>%
  right_join(
    grid_dat %>% select(bi_class, x, y),
    by = "bi_class"
  ) %>%
  replace_na(list(n_cells = 0))

bi_legend_scatter_2 <- ggplot() +
  geom_tile(data = grid_dat,aes(x = x, y = y, fill = bi_fill),colour = "white") +
  scale_fill_identity() +
  geom_text(data = count_df_2,aes(x = x, y = y, label = n_cells),colour   = "white",fontface = "bold", size = 1.8)+
  scale_x_continuous(breaks = 1:4,labels = c("Low", "", "", "High"),expand = c(0.015, 0.015)) +
  scale_y_continuous(breaks = 1:4,labels = c("Low", "", "", "High"),expand = c(0.015, 0.015)) +
  labs(
    # x = expression(S[Linnaean] %->% ""), 
    # y = expression(S[Darwinian] %->% "")
    x = expression(LS %->% ""), 
    y = expression(DS %->% "")
  ) +
  coord_equal(clip = "off") +
  annotate("text",x = Inf, y = 4.8, hjust = 1.2, vjust = 0,size  = 2.2,label = deparse(label_expr),parse = TRUE)+
  theme_void(base_size = 8) +
  theme(
    axis.title  = element_text(size = 8, colour = "black", face = "bold"),
    axis.text   = element_text(size = 6, colour = "black"),
    axis.title.y.left = element_text(angle = 90),
    axis.text.y.left = element_text(angle = 90),
    plot.margin = margin(t = 20, r = 3, b = 3, l = 3),
    panel.ontop = TRUE
  )

bi_legend_scatter_2


map2 <- ggplot(data = d2) +
  ggrastr::rasterise(geom_sf(aes(fill = bi_class),colour = "white", linewidth = 0.03, show.legend = FALSE),dpi = 300)+
  bi_scale_fill(pal = "PurpleOr",dim = 4)+
  bi_theme()+
  coord_sf(crs = "+proj=moll +lon_0=0", ylim = y_limits, expand = FALSE)


pp2 <- ggdraw() +
  draw_plot(map2, 0, 0, 1, 1) +
  draw_plot(bi_legend_scatter_2, -0.1, 0.1, 0.55, 0.55)

################################################################################

d3 <- dt %>% select(iso3,SRnoloc,SRnoseq,SRnoloc_bin,SRnoseq_bin) %>%
  na.omit() %>% bi_class(x= SRnoloc_bin,y = SRnoseq_bin, style = "quantile",dim = 4)

d3_pts <- d3 %>% st_drop_geometry() %>%
  left_join(
    leg$layers[[1]][["data"]] %>% select(bi_class, x, y),
    by = "bi_class"
  )

cor_res <- cor.test(
  d3_pts$SRnoloc,
  d3_pts$SRnoseq,
  method = "spearman",
  use    = "complete.obs"
)

rho   <- unname(cor_res$estimate)
pval  <- cor_res$p.value
label_expr <- if (pval < 0.001) {
  bquote(
    "r = " ~ .(sprintf("%.2f", rho)) * ", " *
      italic(p) < 0.001
  )
} else {
  bquote(
    "r = " ~ .(sprintf("%.2f", rho)) * ", " *
      italic(p) == .(signif(pval, 2))
  )
}

count_df_3 <- d3_pts %>%
  count(bi_class, name = "n_cells") %>%
  right_join(
    grid_dat %>% select(bi_class, x, y),
    by = "bi_class"
  ) %>%
  replace_na(list(n_cells = 0))

bi_legend_scatter_3 <- ggplot() +
  geom_tile(data = grid_dat,aes(x = x, y = y, fill = bi_fill),colour = "white") +
  scale_fill_identity() +
  geom_text(data = count_df_3,aes(x = x, y = y, label = n_cells),colour   = "white",fontface = "bold", size = 1.8)+
  scale_x_continuous(breaks = 1:4,labels = c("Low", "", "", "High"),expand = c(0.015, 0.015)) +
  scale_y_continuous(breaks = 1:4,labels = c("Low", "", "", "High"),expand = c(0.015, 0.015)) +
  labs(
    # x = expression(S[Wallacean] %->% ""), 
    # y = expression(S[Darwinian] %->% "")
    x = expression(WS %->% ""), 
    y = expression(DS %->% "")
  ) +
  coord_equal(clip = "off") +
  annotate("text",x = Inf, y = 4.8, hjust = 1.2, vjust = 0,size  = 2.2,label = deparse(label_expr),parse = TRUE)+
  theme_void(base_size = 8) +
  theme(
    axis.title  = element_text(size = 8, colour = "black", face = "bold"),
    axis.text   = element_text(size = 6, colour = "black"),
    axis.title.y.left = element_text(angle = 90),
    axis.text.y.left = element_text(angle = 90),
    plot.margin = margin(t = 20, r = 3, b = 3, l = 3),
    panel.ontop = TRUE
  )

bi_legend_scatter_3


map3 <- ggplot(data = d3) +
  ggrastr::rasterise(geom_sf(aes(fill = bi_class),colour = "white", linewidth = 0.03, show.legend = FALSE),dpi = 300)+
  bi_scale_fill(pal = "PurpleOr",dim = 4)+
  bi_theme()+
  coord_sf(crs = "+proj=moll +lon_0=0", ylim = y_limits, expand = FALSE)

pp3 <- ggdraw() +
  draw_plot(map3, 0, 0, 1, 1) +
  draw_plot(bi_legend_scatter_3, -0.1, 0.1, 0.55, 0.55)
################################################################################
library(cowplot)
pp1 <- pp1 + theme(plot.margin = margin(0, 0, 0, 0))
pp2 <- pp2 + theme(plot.margin = margin(0, 0, 0, 0))
pp3 <- pp3 + theme(plot.margin = margin(0, 0, 0, 0))

p_combined <- plot_grid(
  pp1, pp2, pp3,
  ncol = 1,
  align = "v",
  axis  = "lr",
  labels = c("A", "B", "C"),
  label_size = 10,
  label_fontface = "bold",
  label_x = 0.1,   
  label_y = 0.95   
)

p_combined

ggsave("figures/supplement/Figure_S1.png",width = 15, height = 20, units = "cm",dpi = 300)

