# ------------------------------------------------------------------------------
# Supplementary Figure S4
# Counts of significant predictor effects across biogeographic realms

library(brms)
library(dplyr)
library(purrr)
library(tidybayes)
library(ggplot2)
options(warn = -1)


compute_realm_slopes <- function(v, draws, cont_levels, baseline_realm = "Afrotropic") {
  # Name of the main-effect coefficient for predictor v
  base_name <- paste0("b_", v)
  
  if (!base_name %in% names(draws)) {
    stop("No coefficient found for ", base_name)
  }
  
  # Baseline slope β_v (effect in the reference realm)
  base <- draws[[base_name]]
  
  # Keep .draw (and optionally .chain, .iteration) to preserve posterior structure
  meta_cols <- draws %>% dplyr::select(.chain, .iteration, .draw)
  
  # For each realm, construct the posterior slope
  out <- map_dfr(cont_levels, function(ct) {
    if (ct == baseline_realm) {
      # Baseline realm: slope = β_v
      slope <- base
    } else {
      # Non-baseline realms: slope = β_v + β_v:realm=ct
      inter_name <- paste0("b_", v, ":biogeographic_realm", ct)
      if (inter_name %in% names(draws)) {
        slope <- base + draws[[inter_name]]
      } else {
        # Safety fallback: no interaction term found (should not occur given your formula)
        slope <- base
      }
    }
    
    tibble(
      realm = ct,
      slope     = slope
    )
  })
  
  # Bind back posterior meta info and add variable name
  bind_cols(
    meta_cols[rep(1:nrow(meta_cols), times = length(cont_levels)), ],
    out
  ) %>%
    mutate(variable = v)
}


compute_continent_slopes <- function(v, draws, cont_levels, baseline_continent = "Africa") {
  # Name of the main-effect coefficient for predictor v
  base_name <- paste0("b_", v)
  
  if (!base_name %in% names(draws)) {
    stop("No coefficient found for ", base_name)
  }
  
  # Baseline slope β_v (effect in the reference continent)
  base <- draws[[base_name]]
  
  # Keep .draw (and optionally .chain, .iteration) to preserve posterior structure
  meta_cols <- draws %>% dplyr::select(.chain, .iteration, .draw)
  
  # For each continent, construct the posterior slope
  out <- map_dfr(cont_levels, function(ct) {
    if (ct == baseline_continent) {
      # Baseline continent: slope = β_v
      slope <- base
    } else {
      # Non-baseline continents: slope = β_v + β_v:continent=ct
      inter_name <- paste0("b_", v, ":continent", ct)
      if (inter_name %in% names(draws)) {
        slope <- base + draws[[inter_name]]
      } else {
        # Safety fallback: no interaction term found (should not occur given your formula)
        slope <- base
      }
    }
    
    tibble(
      continent = ct,
      slope     = slope
    )
  })
  
  # Bind back posterior meta info and add variable name
  bind_cols(
    meta_cols[rep(1:nrow(meta_cols), times = length(cont_levels)), ],
    out
  ) %>%
    mutate(variable = v)
}

################################################################################
basin_linnaean <- readRDS("output/model/basin_linnaean_all.rds")$lognormal
draws_basin_linnaean <- as_draws_df(basin_linnaean)
dat_basin_linnaean <- basin_linnaean$data
levels(dat_basin_linnaean$biogeographic_realm)
cont_levels <- levels(dat_basin_linnaean$biogeographic_realm)
vars_basin_linnaean <- c(
  "body_size", "taxonmic_effort", "taxonomic_activity", "watershed_area",
  "range_size", "elevation", "latitude", "discharge", "watertemp",
  "preserved_specimen", "sampling_effort","sequencing_effort", "range_rarity", "population_density"
)


slopes_draws_basin_linnaean <- map_dfr(
  vars_basin_linnaean,
  ~ compute_realm_slopes(
    v             = .x,
    draws         = draws_basin_linnaean,
    cont_levels   = cont_levels,
    baseline_realm = "Afrotropic"  # 
  )
)

table_basin_linnaean <- slopes_draws_basin_linnaean %>%
  mutate(slope = as.numeric(slope)) %>%
  group_by(variable, realm) %>%
  summarise(
    median   = median(slope, na.rm = TRUE),
    lower_95 = quantile(slope, 0.025, na.rm = TRUE),
    upper_95 = quantile(slope, 0.975, na.rm = TRUE),
    prob_pos = mean(slope > 0, na.rm = TRUE),   # P(β > 0)
    .groups  = "drop"
  ) %>%
  arrange(variable, realm)


sig_basin_linnaean <- table_basin_linnaean %>%
  mutate(
    effect_sign = case_when(
      lower_95 > 0 ~ "positive",
      upper_95 < 0 ~ "negative",
      TRUE         ~ "nonsignificant"
    )
  )

effect_counts_basin_linnaean <- sig_basin_linnaean %>%
  group_by(variable, effect_sign) %>%
  summarise(
    n = n(), 
    .groups = "drop"
  ) %>%
  filter(effect_sign != "nonsignificant") %>%
  mutate(
    n_signed = if_else(effect_sign == "positive", n, -n),
    shortfall = "Linnaean"
  )
rm(basin_linnaean,dat_basin_linnaean,draws_basin_linnaean,
   sig_basin_linnaean,slopes_draws_basin_linnaean,table_basin_linnaean)

################################################################################
basin_wallacean <- readRDS("output/model/basin_wallacean_all.rds")$lognormal
draws_basin_wallacean <- as_draws_df(basin_wallacean)
dat_basin_wallacean <- basin_wallacean$data
levels(dat_basin_wallacean$biogeographic_realm)
cont_levels <- levels(dat_basin_wallacean$biogeographic_realm)
vars_basin_wallacean <- c(
  "body_size", "taxonmic_effort", "taxonomic_activity", "watershed_area",
  "range_size", "elevation", "latitude", "discharge", "watertemp",
  "preserved_specimen", "sequencing_effort", "range_rarity", "population_density"
)


slopes_draws_basin_wallacean <- map_dfr(
  vars_basin_wallacean,
  ~ compute_realm_slopes(
    v             = .x,
    draws         = draws_basin_wallacean,
    cont_levels   = cont_levels,
    baseline_realm = "Afrotropic"  # 
  )
)

table_basin_wallacean <- slopes_draws_basin_wallacean %>%
  mutate(slope = as.numeric(slope)) %>%
  group_by(variable, realm) %>%
  summarise(
    median   = median(slope, na.rm = TRUE),
    lower_95 = quantile(slope, 0.025, na.rm = TRUE),
    upper_95 = quantile(slope, 0.975, na.rm = TRUE),
    prob_pos = mean(slope > 0, na.rm = TRUE),   # P(β > 0)
    .groups  = "drop"
  ) %>%
  arrange(variable, realm)


sig_basin_wallacean <- table_basin_wallacean %>%
  mutate(
    effect_sign = case_when(
      lower_95 > 0 ~ "positive",
      upper_95 < 0 ~ "negative",
      TRUE         ~ "nonsignificant"
    )
  )

effect_counts_basin_wallacean <- sig_basin_wallacean %>%
  group_by(variable, effect_sign) %>%
  summarise(
    n = n(), 
    .groups = "drop"
  ) %>%
  filter(effect_sign != "nonsignificant") %>%
  mutate(
    n_signed = if_else(effect_sign == "positive", n, -n),
    shortfall = "Wallacean"
  )
rm(basin_wallacean,dat_basin_wallacean,draws_basin_wallacean,
   sig_basin_wallacean,slopes_draws_basin_wallacean,table_basin_wallacean)

################################################################################
basin_darwinian <- readRDS("output/model/basin_darwinian_all.rds")$gamma
draws_basin_darwinian <- as_draws_df(basin_darwinian)
dat_basin_darwinian <- basin_darwinian$data
levels(dat_basin_darwinian$biogeographic_realm)
cont_levels <- levels(dat_basin_darwinian$biogeographic_realm)
vars_basin_darwinian <- c(
  "body_size", "taxonmic_effort", "taxonomic_activity", "watershed_area",
  "range_size", "elevation", "latitude", "discharge", "watertemp",
  "preserved_specimen", "sampling_effort", "range_rarity", "population_density"
)


slopes_draws_basin_darwinian <- map_dfr(
  vars_basin_darwinian,
  ~ compute_realm_slopes(
    v             = .x,
    draws         = draws_basin_darwinian,
    cont_levels   = cont_levels,
    baseline_realm = "Afrotropic"  # 
  )
)

table_basin_darwinian <- slopes_draws_basin_darwinian %>%
  mutate(slope = as.numeric(slope)) %>%
  group_by(variable, realm) %>%
  summarise(
    median   = median(slope, na.rm = TRUE),
    lower_95 = quantile(slope, 0.025, na.rm = TRUE),
    upper_95 = quantile(slope, 0.975, na.rm = TRUE),
    prob_pos = mean(slope > 0, na.rm = TRUE),   # P(β > 0)
    .groups  = "drop"
  ) %>%
  arrange(variable, realm)


sig_basin_darwinian <- table_basin_darwinian %>%
  mutate(
    effect_sign = case_when(
      lower_95 > 0 ~ "positive",
      upper_95 < 0 ~ "negative",
      TRUE         ~ "nonsignificant"
    )
  )

effect_counts_basin_darwinian <- sig_basin_darwinian %>%
  group_by(variable, effect_sign) %>%
  summarise(
    n = n(), 
    .groups = "drop"
  ) %>%
  filter(effect_sign != "nonsignificant") %>%
  mutate(
    n_signed = if_else(effect_sign == "positive", n, -n),
    shortfall = "Darwinian"
  )
rm(basin_darwinian,dat_basin_darwinian,draws_basin_darwinian,
   sig_basin_darwinian,slopes_draws_basin_darwinian,table_basin_darwinian)

#################################################################################
effect_counts_basin <- rbind(effect_counts_basin_linnaean,
                             effect_counts_basin_wallacean,
                             effect_counts_basin_darwinian)

rename_map_basin <- c(
  body_size            = "Body size",
  preserved_specimen   = "N. preserved specimens",
  range_size           = "Range size",
  range_rarity         = "Range rarity",
  latitude             = "Latitude",
  watershed_area       = "Watershed area",
  discharge            = "Streamflow",
  watertemp            = "Water temperature",
  elevation            = "Elevation",
  population_density   = "Human density",
  taxonmic_effort      = "Taxonomic effort",
  taxonomic_activity   = "Taxonomic activity",
  sampling_effort      = "Sampling effort",
  sequencing_effort    = "Sequencing effort"
)

effect_counts_plot_1 <- effect_counts_basin %>%
  mutate(
    variable = recode(variable, !!!rename_map_basin),
    variable = factor(variable, levels = rev(unname(rename_map_basin))),
    shortfall = factor(shortfall,
                       levels = c("Linnaean", "Wallacean", "Darwinian")),
    count = abs(n_signed),
    y_pos = n_signed + 0.25 * sign(n_signed)
  )


nc_colors <- c(
  positive = "#4C78A8",  
  negative = "#F58518"   
)


custom_labels <- c("Linnaean" = "Description", "Wallacean" = "Geolocation", "Darwinian" = "Sequencing")
p1 <- ggplot(effect_counts_plot_1,
             aes(x = variable, y = n_signed, fill = effect_sign)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linetype = 2,
             color = "grey50", linewidth = 0.4) +
  geom_text(
    aes(y = y_pos, label = count),
    size = 3
  ) +
  scale_y_continuous(labels = function(x) abs(x)) +
  coord_flip() +
  scale_fill_manual(
    values = nc_colors,
    breaks = c("positive", "negative"),
    labels = c("Positive", "Negative"),
    name   = "Direction"
  ) +
  labs(
    x = NULL,
    y = expression("Number of biogeographic realm with significant effects (" * italic(p) * " < 0.05)")
  ) +
  annotate("text", x = 0.8, y = 6, label = "+",colour = "#4C78A8", face = "bold",size = 5)+
  annotate("text", x = 0.8, y = -6, label = "-",colour = "#F58518", face = "bold",size = 5)+
  facet_wrap(. ~ shortfall, labeller = labeller(shortfall = custom_labels)) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.title.y  = element_text(size = 8,colour = "black"),
    axis.text.y   = element_text(size = 8,colour = "black"),
    legend.position = "top",
    legend.title    = element_text(size = 8,face = "bold",colour = "black"),
    legend.text     = element_text(size = 8)
  )


################################################################################
country_linnaean <- readRDS("output/model/country_linnaean_all.rds")$weibull
draws_country_linnaean <- as_draws_df(country_linnaean)
dat_country_linnaean <- country_linnaean$data
levels(dat_country_linnaean$continent)
cont_levels <- levels(dat_country_linnaean$continent)
vars_country_linnaean <- c(
  "body_size", "taxonmic_effort", "taxonomic_activity", "country_area",
  "range_size", "elevation", "latitude", "discharge", "watertemp",
  "preserved_specimen", "sampling_effort","sequencing_effort", "range_rarity", "population_density"
)


slopes_draws_country_linnaean <- map_dfr(
  vars_country_linnaean,
  ~ compute_continent_slopes(
    v             = .x,
    draws         = draws_country_linnaean,
    cont_levels   = cont_levels,
    baseline_continent = "Africa"  # 如果 baseline 不是 Africa，可以改这里
  )
)

table_country_linnaean <- slopes_draws_country_linnaean %>%
  mutate(slope = as.numeric(slope)) %>%
  group_by(variable, continent) %>%
  summarise(
    median   = median(slope, na.rm = TRUE),
    lower_95 = quantile(slope, 0.025, na.rm = TRUE),
    upper_95 = quantile(slope, 0.975, na.rm = TRUE),
    prob_pos = mean(slope > 0, na.rm = TRUE),   # P(β > 0)
    .groups  = "drop"
  ) %>%
  arrange(variable, continent)


sig_country_linnaean <- table_country_linnaean %>%
  mutate(
    effect_sign = case_when(
      lower_95 > 0 ~ "positive",
      upper_95 < 0 ~ "negative",
      TRUE         ~ "nonsignificant"
    )
  )

effect_counts_country_linnaean <- sig_country_linnaean %>%
  group_by(variable, effect_sign) %>%
  summarise(
    n = n(), 
    .groups = "drop"
  ) %>%
  filter(effect_sign != "nonsignificant") %>%
  mutate(
    n_signed = if_else(effect_sign == "positive", n, -n),
    shortfall = "Linnaean"
  )
rm(country_linnaean,dat_country_linnaean,draws_country_linnaean,
   sig_country_linnaean,slopes_draws_country_linnaean,table_country_linnaean)
################################################################################
country_wallacean <- readRDS("output/model/country_wallacean_all.rds")$lognormal
draws_country_wallacean <- as_draws_df(country_wallacean)
dat_country_wallacean <- country_wallacean$data
levels(dat_country_wallacean$continent)
cont_levels <- levels(dat_country_wallacean$continent)
vars_country_wallacean <- c(
  "body_size", "taxonmic_effort", "taxonomic_activity", "country_area",
  "range_size", "elevation", "latitude", "discharge", "watertemp",
  "preserved_specimen", "sequencing_effort", "range_rarity", "population_density"
)


slopes_draws_country_wallacean <- map_dfr(
  vars_country_wallacean,
  ~ compute_continent_slopes(
    v             = .x,
    draws         = draws_country_wallacean,
    cont_levels   = cont_levels,
    baseline_continent = "Africa"  
  )
)

table_country_wallacean <- slopes_draws_country_wallacean %>%
  mutate(slope = as.numeric(slope)) %>%
  group_by(variable, continent) %>%
  summarise(
    median   = median(slope, na.rm = TRUE),
    lower_95 = quantile(slope, 0.025, na.rm = TRUE),
    upper_95 = quantile(slope, 0.975, na.rm = TRUE),
    prob_pos = mean(slope > 0, na.rm = TRUE),   # P(β > 0)
    .groups  = "drop"
  ) %>%
  arrange(variable, continent)

sig_country_wallacean <- table_country_wallacean %>%
  mutate(
    effect_sign = case_when(
      lower_95 > 0 ~ "positive",
      upper_95 < 0 ~ "negative",
      TRUE         ~ "nonsignificant"
    )
  )

effect_counts_country_wallacean <- sig_country_wallacean %>%
  group_by(variable, effect_sign) %>%
  summarise(
    n = n(), 
    .groups = "drop"
  ) %>%
  filter(effect_sign != "nonsignificant") %>%
  mutate(
    n_signed = if_else(effect_sign == "positive", n, -n),
    shortfall = "Wallacean"
  )
rm(country_wallacean,dat_country_wallacean,draws_country_wallacean,
   sig_country_wallacean,slopes_draws_country_wallacean,table_country_wallacean)

################################################################################
country_darwinian <- readRDS("output/model/country_darwinian_all.rds")$weibull
draws_country_darwinian <- as_draws_df(country_darwinian)
dat_country_darwinian <- country_darwinian$data
levels(dat_country_darwinian$continent)
cont_levels <- levels(dat_country_darwinian$continent)
vars_country_darwinian <- c(
  "body_size", "taxonmic_effort", "taxonomic_activity", "country_area",
  "range_size", "elevation", "latitude", "discharge", "watertemp",
  "preserved_specimen", "sampling_effort", "range_rarity", "population_density"
)


slopes_draws_country_darwinian <- map_dfr(
  vars_country_darwinian,
  ~ compute_continent_slopes(
    v             = .x,
    draws         = draws_country_darwinian,
    cont_levels   = cont_levels,
    baseline_continent = "Africa"  
  )
)

table_country_darwinian <- slopes_draws_country_darwinian %>%
  mutate(slope = as.numeric(slope)) %>%
  group_by(variable, continent) %>%
  summarise(
    median   = median(slope, na.rm = TRUE),
    lower_95 = quantile(slope, 0.025, na.rm = TRUE),
    upper_95 = quantile(slope, 0.975, na.rm = TRUE),
    prob_pos = mean(slope > 0, na.rm = TRUE),   # P(β > 0)
    .groups  = "drop"
  ) %>%
  arrange(variable, continent)


sig_country_darwinian <- table_country_darwinian %>%
  mutate(
    effect_sign = case_when(
      lower_95 > 0 ~ "positive",
      upper_95 < 0 ~ "negative",
      TRUE         ~ "nonsignificant"
    )
  )

effect_counts_country_darwinian <- sig_country_darwinian %>%
  group_by(variable, effect_sign) %>%
  summarise(
    n = n(), 
    .groups = "drop"
  ) %>%
  filter(effect_sign != "nonsignificant") %>%
  mutate(
    n_signed = if_else(effect_sign == "positive", n, -n),
    shortfall = "Darwinian"
  )
rm(country_darwinian,dat_country_darwinian,draws_country_darwinian,
   sig_country_darwinian,slopes_draws_country_darwinian,table_country_darwinian)

################################################################################
effect_counts_country <- rbind(effect_counts_country_linnaean,
                               effect_counts_country_wallacean,
                               effect_counts_country_darwinian)

rename_map_country <- c(
  body_size            = "Body size",
  preserved_specimen   = "N. preserved specimens",
  range_size           = "Range size",
  range_rarity         = "Range rarity",
  latitude             = "Latitude",
  country_area         = "Country area",
  discharge            = "Streamflow",
  watertemp            = "Water temperature",
  elevation            = "Elevation",
  population_density   = "Human density",
  taxonmic_effort      = "Taxonomic effort",
  taxonomic_activity   = "Taxonomic activity",
  sampling_effort      = "Sampling effort",
  sequencing_effort    = "Sequencing effort"
)

effect_counts_plot_2 <- effect_counts_country %>%
  mutate(
    variable = recode(variable, !!!rename_map_country),
    variable = factor(variable, levels = rev(unname(rename_map_country))),
    shortfall = factor(shortfall,
                       levels = c("Linnaean", "Wallacean", "Darwinian")),
    count = abs(n_signed),
    y_pos = n_signed + 0.25 * sign(n_signed)
  )


nc_colors <- c(
  positive = "#4C78A8",  
  negative = "#F58518"   
)



p2 <- ggplot(effect_counts_plot_2,
             aes(x = variable, y = n_signed, fill = effect_sign)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 0, linetype = 2,
             color = "grey50", linewidth = 0.4) +
  geom_text(
    aes(y = y_pos, label = count),
    size = 3
  ) +
  scale_y_continuous(labels = function(x) abs(x)) +
  coord_flip() +
  scale_fill_manual(
    values = nc_colors,
    breaks = c("positive", "negative"),
    labels = c("Positive", "Negative"),
    name   = "Direction"
  ) +
  labs(
    x = NULL,
    y = expression("Number of continent with significant effects (" * italic(p) * " < 0.05)")
  ) +
  annotate("text", x = 0.8, y = 6, label = "+",colour = "#4C78A8", face = "bold",size = 5)+
  annotate("text", x = 0.8, y = -6, label = "-",colour = "#F58518", face = "bold",size = 5)+
  facet_wrap(. ~ shortfall, labeller = labeller(shortfall = custom_labels)) +
  theme_bw() +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.title.y  = element_text(size = 8,colour = "black"),
    axis.text.y   = element_text(size = 8,colour = "black"),
    legend.position = "none",
    legend.title    = element_text(size = 8,face = "bold",colour = "black"),
    legend.text     = element_text(size = 8)
  )

library(patchwork)  
p1 / p2 + plot_annotation(tag_levels = "A", tag_prefix = "(", tag_suffix = ")")


ggsave("figures/supplement/Figure_S4.png",
       units = "cm",dpi = 300,
       width = 18, height = 18)
