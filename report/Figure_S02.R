# ------------------------------------------------------------------------------
# Supplementary Figure S2
# Relationship between basin-level probability of lacking sequences
# (Darwinian shortfall component) and phylogenetic deficit
# (PD-based evolutionary incompleteness).
# ------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)   # coarse global land polygons
library(cowplot)         # plot composition
library(ggrastr)         # rasterise dense vector layers for faster plotting
options(warms = -1)
options(sf_use_s2 = FALSE)

# ------------------------------------------------------------------------------
# 1. Load and prepare basin-level Darwinian shortfall and PD deficit data
# ------------------------------------------------------------------------------

darwinian_tbl <- read.csv("output/tables/basin_darwinian_shortfall.csv")
pd_deficit_tbl <- read.csv("output/tables/pd_deficits_basin_id.csv")

# Compute mean PD deficit per basin
pd_deficit_basin <- pd_deficit_tbl %>%
  group_by(basin_id) %>%
  summarise(
    mean_pd_deficit = mean(PD_deficit, na.rm = TRUE),
    .groups = "drop"
  )

# Simple 0–1 rescaling with a small offset to avoid division by zero
min_max_normalize_safe <- function(x, epsilon = 1e-4) {
  min_val <- min(x,na.rm = T)
  max_val <- max(x,na.rm = T)
  range_val <- max_val - min_val + epsilon
  return((x - min_val) / range_val)
}

# Join PD deficit and Darwinian shortfall; normalise to common 0–1 scale
basin_stats <- pd_deficit_basin %>%
  left_join(darwinian_tbl[, c(1, 6)], by = "basin_id") %>%
  mutate(
    mean_pd_deficit  = min_max_normalize_safe(mean_pd_deficit),
    SRnoseq  = min_max_normalize_safe(log10(SRnoseq+1))
  ) %>%
  tidyr::drop_na()

# ------------------------------------------------------------------------------
# 2. Load inland basin polygons and harmonise biogeographic realm labels
# ------------------------------------------------------------------------------

inland_basins <- readRDS("input/raw/basin/basin_sf_v1.rds")

inland_basins <- inland_basins %>%
  mutate(
    biogeographic_realm = dplyr::case_when(
      biogeographic_realm == "Australasia" ~ "Australasian",
      biogeographic_realm == "Oceania"     ~ "Oceanian",
      TRUE ~ biogeographic_realm
    )
  )

realm_levels <- c(
  "Nearctic",
  "Neotropic",
  "Palearctic",
  "Afrotropic",
  "Indomalayan",
  "Australasian",
  "Oceanian"
)

inland_basins$biogeographic_realm <- factor(
  inland_basins$biogeographic_realm,
  levels = rev(realm_levels)
)

# Shared colour palette for biogeographic realms
realm_cols <- c(
  "Afrotropic"   = "#D55E00",
  "Australasian" = "#0072B2",
  "Indomalayan"  = "#009E73",
  "Nearctic"     = "#F0E442",
  "Neotropic"    = "#E69F00",
  "Oceanian"     = "#56B4E9",
  "Palearctic"   = "#CC79A7"
)

# ------------------------------------------------------------------------------
# 3. Basemap for inset: Mollweide projection with biogeographic realms
# ------------------------------------------------------------------------------

world_sf <- ne_countries(scale = "small", returnclass = "sf") %>%
  st_wrap_dateline(options = c("WRAPDATELINE=YES"))

moll_crs <- st_crs("+proj=moll")

# Compute y-limits in Mollweide to constrain map extent
lat_band <- data.frame(lon = 0, lat = c(-60, 90)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs = moll_crs)

y_limits <- st_coordinates(lat_band)[, "Y"]

map_inset <- ggplot() +
  rasterise(geom_sf(data = world_sf, fill = "gray80", colour = NA, linewidth = 0.2),dpi = 300) +
  rasterise(geom_sf(data  = inland_basins,aes(fill = biogeographic_realm),alpha = 0.7,colour = NA,linewidth = 0.5,show.legend = FALSE),dpi = 300) +
  scale_fill_manual(values = realm_cols) +
  coord_sf(crs = "+proj=moll +lon_0=0",ylim  = y_limits,expand = FALSE) +
  theme_void()

# ------------------------------------------------------------------------------
# 4. Join basin statistics to geometries and compute correlations
# ------------------------------------------------------------------------------

basin_sf <- inland_basins %>%
  left_join(basin_stats, by = "basin_id")

# Global Spearman correlation between PD deficit and probability of missing sequences
cor_global_basin <- cor.test(
  basin_stats$mean_pd_deficit,
  basin_stats$SRnoseq,
  method = "spearman",
  exact  = FALSE,
  use    = "complete.obs"
)

# Realm-wise Spearman correlations
cor_by_realm <- basin_sf %>%
  st_drop_geometry() %>%
  group_by(biogeographic_realm) %>%
  summarise(
    rho = suppressWarnings(
      cor(mean_pd_deficit, SRnoseq,
          method = "spearman", use = "complete.obs")
    ),
    p   = suppressWarnings(
      cor.test(mean_pd_deficit, SRnoseq,
               method = "spearman", use = "complete.obs")$p.value
    ),
    n   = sum(complete.cases(mean_pd_deficit, SRnoseq)),
    .groups = "drop"
  )

# Label positions for realm-wise correlation summary
# (placed in the upper-left quadrant of the scatter plot)
label_df <- cor_by_realm %>%
  arrange(factor(biogeographic_realm, levels = realm_levels)) %>%
  mutate(
    x = 0,                                        # fixed left margin in 0–1 space
    y = seq(0.95, 0.75, length.out = n()),          # vertically staggered labels
    label = paste0(
      biogeographic_realm,
      ": r = ", sprintf("%.2f", rho),
      "; p = ", format.pval(p, digits = 2),
      "; n = ", n
    )
  )

# ------------------------------------------------------------------------------
# 5. Scatter plot of PD deficit vs. probability of missing sequences
# ------------------------------------------------------------------------------

scatter_df <- basin_sf %>% st_drop_geometry() %>%
  filter(
    is.finite(mean_pd_deficit),
    is.finite(SRnoseq)
  )
scatter_df$biogeographic_realm <- factor(scatter_df$biogeographic_realm, levels = realm_levels)
label_df$biogeographic_realm <- factor(label_df$biogeographic_realm, levels = realm_levels)
p_scatter <- ggplot(scatter_df,aes(x = mean_pd_deficit, y = SRnoseq)) +
  geom_jitter(aes(fill = biogeographic_realm),width  = 0.01,height = 0.01,alpha  = 0.5,
              size   = 1.8,shape  = 21,colour = "black",show.legend = FALSE) +
  geom_smooth(aes(colour = biogeographic_realm),method      = "lm",
              se = FALSE,show.legend = FALSE) +
  geom_text(data = label_df,aes(x = x,y = y,label = label,colour = biogeographic_realm),
            inherit.aes = FALSE,hjust = 0,vjust = 1,size = 2,show.legend = FALSE) +
  scale_fill_manual(values = realm_cols) +
  scale_colour_manual(values = realm_cols) +
  labs(
    x     = "Phylogenetic diversity deficit",
    y     = "Undocumented species (molecular data)"
    #title = "Complementary dimensions of the Darwinian shortfall"
  ) +
  annotate(
    "text",
    x      = Inf,
    y      = Inf,
    hjust  = 1.1,
    vjust  = 1.5,
    size   = 3,
    label  = paste0(
      "Global: r = ",
      sprintf("%.2f", cor_global_basin$estimate),
      "; p < 0.001",
      "; n = ", nrow(stats::na.omit(scatter_df))
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 10)

# ------------------------------------------------------------------------------
# 6. Compose scatter plot with realm map inset and save to file
# ------------------------------------------------------------------------------

fig_darwinian_basin <- ggdraw() +
  draw_plot(p_scatter, 0, 0, 1, 1) +
  draw_plot(map_inset, 0.55, 0.65, 0.4, 0.4)


# ------------------------------------------------------------------------------
# Interpretation (for caption / text)
# ------------------------------------------------------------------------------
# Basins show a weak but significant negative relationship between
# phylogenetic diversity deficit and the probability of missing sequences
# (Spearman ρ ≈ –0.17). The broad scatter and shallow trend indicate that
# phylogenetic and molecular knowledge gaps are only weakly aligned, highlighting
# their complementary spatial patterns.
#
# Phylogenetic diversity loss and molecular data gaps do not systematically
# co-occur, underscoring the multidimensional nature of the Darwinian shortfall.
# ------------------------------------------------------------------------------
min_max_normalize_safe <- function(x, epsilon = 1e-4) {
  min_val <- min(x,na.rm = T)
  max_val <- max(x,na.rm = T)
  range_val <- max_val - min_val + epsilon
  return((x - min_val) / range_val)
}

country_darwinian_tbl <- read.csv("output/tables/country_darwinian_shortfall.csv")
country_pd_deficit_tbl <- read.csv("output/tables/pd_deficits_iso3.csv")

# Compute mean PD deficit per basin
pd_deficit_country <- country_pd_deficit_tbl %>%
  group_by(iso3) %>%
  summarise(
    mean_pd_deficit = mean(PD_deficit, na.rm = TRUE),
    .groups = "drop"
  )


# Join PD deficit and Darwinian shortfall; normalise to common 0–1 scale
country_stats <- pd_deficit_country %>%
  left_join(country_darwinian_tbl[, c(1, 6)], by = "iso3") %>%
  mutate(
    mean_pd_deficit  = min_max_normalize_safe(mean_pd_deficit),
    SRnoseq  = min_max_normalize_safe(log10(SRnoseq+1))
  ) %>%
  tidyr::drop_na()

inland_country <- readRDS("input/raw/country.rds")
country <- read.csv("input/raw/country_list.csv")
inland_country <- inland_country %>% left_join(country[,c(1,5)], by = "iso3")
continent_levels <- c(
  "North America",
  "South America",
  "Europe",
  "Africa",
  "Asia",
  "Oceania"
)

inland_country$continent <- factor(inland_country$continent,levels = continent_levels)

# Shared colour palette for biogeographic realms
continent_cols <- c(
  "Africa"   = "#D55E00",
  "Oceania" = "#0072B2",
  "Asia"  = "#009E73",
  "North America"     = "#F0E442",
  "South America"    = "#E69F00",
  "Europe"   = "#CC79A7"
)

moll_crs <- st_crs("+proj=moll")

# Compute y-limits in Mollweide to constrain map extent
lat_band <- data.frame(lon = 0, lat = c(-60, 90)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(crs = moll_crs)

y_limits <- st_coordinates(lat_band)[, "Y"]

map_inset_country <- ggplot() +
  rasterise(geom_sf(data= inland_country,aes(fill = continent),alpha = 0.7,
                    colour = NA,linewidth = 0.5,show.legend = FALSE),dpi = 300) +
  scale_fill_manual(values = continent_cols) +
  coord_sf(crs = "+proj=moll +lon_0=0",ylim = y_limits,expand = FALSE) +
  theme_void()

country_sf <- inland_country %>%
  left_join(country_stats, by = "iso3")

# Global Spearman correlation between PD deficit and probability of missing sequences
cor_global_country <- cor.test(
  country_stats$mean_pd_deficit,
  country_stats$SRnoseq,
  method = "spearman",
  exact  = FALSE,
  use    = "complete.obs"
)

# Realm-wise Spearman correlations
cor_by_continent <- country_sf %>%
  st_drop_geometry() %>%
  select(continent,mean_pd_deficit,SRnoseq) %>%
  na.omit() %>%
  group_by(continent) %>%
  summarise(
    rho = suppressWarnings(
      cor(mean_pd_deficit, SRnoseq,
          method = "spearman", use = "complete.obs")
    ),
    p   = suppressWarnings(
      cor.test(mean_pd_deficit, SRnoseq,
               method = "spearman", use = "complete.obs")$p.value
    ),
    n   = sum(complete.cases(mean_pd_deficit, SRnoseq)),
    .groups = "drop"
  )

# Label positions for realm-wise correlation summary
# (placed in the upper-left quadrant of the scatter plot)
label_df_continent <- cor_by_continent %>%
  arrange(factor(continent, levels = continent_levels)) %>%
  mutate(
    x = 0,                                        # fixed left margin in 0–1 space
    y = seq(0.95, 0.75, length.out = n()),          # vertically staggered labels
    label = paste0(
      continent,
      ": r = ", sprintf("%.2f", rho),
      "; p = ", format.pval(p, digits = 2),
      "; n = ", n
    )
  )


scatter_df_continent <- country_sf %>%
  st_drop_geometry() %>%
  filter(
    is.finite(mean_pd_deficit),
    is.finite(SRnoseq)
  )
scatter_df_continent$continent <- factor(scatter_df_continent$continent, levels = continent_levels)
label_df_continent$continent <- factor(label_df_continent$continent, levels = continent_levels)
p_scatter_cpuntry <- ggplot(scatter_df_continent,aes(x = mean_pd_deficit, y = SRnoseq)) +
  geom_jitter(aes(fill = continent),width  = 0.01,height = 0.01,
              alpha  = 0.5,size   = 1.8,shape  = 21,colour = "black",show.legend = FALSE ) +
  geom_smooth(aes(colour = continent),method  = "lm",se = FALSE,show.legend = FALSE ) +
  geom_text(data = label_df_continent,
            aes(x= x,y  = y,label = label,colour = continent),
            inherit.aes = FALSE,hjust= 0,vjust = 1,size = 2,show.legend = FALSE) +
  scale_fill_manual(values = continent_cols) +
  scale_colour_manual(values = continent_cols) +
  labs(
    x     = "Phylogenetic diversity deficit",
    y     = "Undocumented species (molecular data)"
    #title = "Complementary dimensions of the Darwinian shortfall"
  ) +
  annotate("text",x = Inf,y = Inf,hjust  = 1.1,vjust  = 1.5,size   = 3,
           label  = paste0(
             "Global: r = ",
             sprintf("%.2f", cor_global_country$estimate),
             "; p < 0.001",
             "; n = ", nrow(stats::na.omit(scatter_df_continent))
           )
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 10)

fig_darwinian_country <- ggdraw() +
  draw_plot(p_scatter_cpuntry, 0, 0, 1, 1) +
  draw_plot(map_inset_country, 0.55, 0.65, 0.4, 0.4)

library(cowplot)
plot_grid(
  fig_darwinian_basin, 
  fig_darwinian_country,
  ncol = 2,                 
  align = "h",              
  axis  = "tb",             
  labels = c("A", "B"),
  label_size = 10,
  label_fontface = "bold",
  label_x = 0.02,          
  label_y = 0.98
)


ggsave(
  filename = "figures/supplement/Figure_S2.png",
  width    = 21,
  height   = 16,
  units    = "cm",
  dpi      = 300
)



