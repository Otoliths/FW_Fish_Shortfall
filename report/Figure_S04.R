# Supplementary Figure S4
# Counts of significant predictor effects across biogeographic realms and continents
rm(list = ls())
library(brms)
library(dplyr)
library(purrr)
library(tidybayes)
library(ggplot2)
options(warn = -1)

# ---- helper: validate family and transform beta to HR ----
.beta_to_hr <- function(beta,
                        draws,
                        family = c("weibull", "gompertz"),
                        shape_col = "shape") {
  family <- match.arg(family)
  
  if (family == "gompertz") {
    log_hr <- beta
  } else if (family == "weibull") {
    if (!shape_col %in% names(draws)) {
      stop("No Weibull shape column found: ", shape_col)
    }
    shape  <- draws[[shape_col]]
    log_hr <- -shape * beta
  }
  
  list(
    log_hr = log_hr,
    hr     = exp(log_hr)
  )
}

# ---- helper: construct realm-specific beta ----
.compute_realm_beta <- function(v,
                                draws,
                                cont_levels,
                                baseline_realm = "Afrotropic",
                                realm_var = "biogeographic_realm") {
  base_name <- paste0("b_", v)
  
  if (!base_name %in% names(draws)) {
    stop("No coefficient found for ", base_name)
  }
  
  base <- draws[[base_name]]
  
  purrr::map_dfr(cont_levels, function(ct) {
    if (ct == baseline_realm) {
      beta <- base
    } else {
      inter_name <- paste0("b_", v, ":", realm_var, ct)
      
      # fallback: try reversed interaction naming if needed
      inter_name_alt <- paste0("b_", realm_var, ct, ":", v)
      
      if (inter_name %in% names(draws)) {
        beta <- base + draws[[inter_name]]
      } else if (inter_name_alt %in% names(draws)) {
        beta <- base + draws[[inter_name_alt]]
      } else {
        beta <- base
      }
    }
    
    tibble::tibble(
      realm = ct,
      beta  = beta
    )
  })
}

# ---- realm-specific HR ----
compute_realm_hr <- function(v,
                             draws,
                             cont_levels,
                             baseline_realm = "Afrotropic",
                             family = c("weibull", "gompertz"),
                             shape_col = "shape",
                             realm_var = "biogeographic_realm") {
  family <- match.arg(family)
  
  req_meta <- c(".chain", ".iteration", ".draw")
  has_meta <- req_meta[req_meta %in% names(draws)]
  
  if (length(has_meta) == 0) {
    stop("draws must contain at least one posterior meta column such as .draw")
  }
  
  meta_cols <- draws %>% dplyr::select(dplyr::all_of(has_meta))
  
  realm_beta <- .compute_realm_beta(
    v = v,
    draws = draws,
    cont_levels = cont_levels,
    baseline_realm = baseline_realm,
    realm_var = realm_var
  )
  
  trans <- .beta_to_hr(
    beta = realm_beta$beta,
    draws = draws[rep(seq_len(nrow(draws)), times = length(cont_levels)), , drop = FALSE],
    family = family,
    shape_col = shape_col
  )
  
  out <- realm_beta %>%
    dplyr::mutate(
      log_hr = trans$log_hr,
      hr     = trans$hr
    )
  
  dplyr::bind_cols(
    meta_cols[rep(seq_len(nrow(meta_cols)), times = length(cont_levels)), , drop = FALSE],
    out
  ) %>%
    dplyr::mutate(
      variable = v,
      family   = family
    ) %>%
    dplyr::relocate(variable, family, realm, .after = dplyr::last_col())
}

################################################################################
fit <- readRDS("output/model/basin_linnaean_weibull_full.rds")
draws <- as_draws_df(fit)
dat <- fit$data
levels(dat$biogeographic_realm)
w <- prop.table(table(dat$biogeographic_realm))
cont_levels <- names(w)
n_realm <- as.data.frame(table(dat$biogeographic_realm))
names(n_realm) <- c("realm","n")

vars_linnaean <- c(
  "body_size", "taxonomic_effort", "taxonomic_activity", "watershed_area",
  "range_size", "elevation", "latitude", "discharge", "watertemp",
  "preserved_specimen", "sampling_effort","sequencing_effort", "range_rarity", "population_density"
)

effect_counts_basin_linnaean <- map_dfr(
  vars_linnaean,
  ~ compute_realm_hr(
    v             = .x,
    draws         = draws,
    cont_levels   = cont_levels,
    baseline_realm = "Afrotropic",
    family = "weibull"
  )
) %>%
  mutate(
    hr    = as.numeric(hr),
    realm = factor(realm, levels = cont_levels),
    variable  = factor(variable, levels = vars_linnaean)
  ) %>%
  group_by(variable, realm) %>%
  summarise(
    median   = median(hr, na.rm = TRUE),
    lower_95 = quantile(hr, 0.025, na.rm = TRUE),
    upper_95 = quantile(hr, 0.975, na.rm = TRUE),
    prob_pos = mean(hr > 1, na.rm = TRUE),   
    .groups  = "drop"
  ) %>%
  arrange(variable, realm) %>%
  mutate(
    effect_sign = case_when(
      lower_95 > 1 ~ "positive",
      upper_95 < 1 ~ "negative",
      TRUE         ~ "nonsignificant"
    )
  ) %>%
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

rm(fit,w,vars_linnaean,cont_levels,draws,dat,n_realm)

################################################################################
fit <- readRDS("output/model/basin_wallacean_gompertz_full.rds")
draws <- as_draws_df(fit)
dat <- fit$data
levels(dat$biogeographic_realm)
w <- prop.table(table(dat$biogeographic_realm))
cont_levels <- names(w)
n_realm <- as.data.frame(table(dat$biogeographic_realm))
names(n_realm) <- c("realm","n")

vars_wallacean <- c(
  "body_size", "taxonomic_effort", "taxonomic_activity", "watershed_area",
  "range_size", "elevation", "latitude", "discharge", "watertemp",
  "preserved_specimen", "sequencing_effort", "range_rarity", "population_density"
)

effect_counts_basin_wallacean <- map_dfr(
  vars_wallacean,
  ~ compute_realm_hr(
    v             = .x,
    draws         = draws,
    cont_levels   = cont_levels,
    baseline_realm = "Afrotropic",
    family = "gompertz"
  )
) %>%
  mutate(
    hr    = as.numeric(hr),
    realm = factor(realm, levels = cont_levels),
    variable  = factor(variable, levels = vars_wallacean)
  ) %>%
  group_by(variable, realm) %>%
  summarise(
    median   = median(hr, na.rm = TRUE),
    lower_95 = quantile(hr, 0.025, na.rm = TRUE),
    upper_95 = quantile(hr, 0.975, na.rm = TRUE),
    prob_pos = mean(hr > 1, na.rm = TRUE),   
    .groups  = "drop"
  ) %>%
  arrange(variable, realm) %>%
  mutate(
    effect_sign = case_when(
      lower_95 > 1 ~ "positive",
      upper_95 < 1 ~ "negative",
      TRUE         ~ "nonsignificant"
    )
  ) %>%
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

rm(fit,w,vars_wallacean,cont_levels,draws,dat,n_realm)

################################################################################
fit <- readRDS("output/model/basin_darwinian_gompertz_full.rds")
draws <- as_draws_df(fit)
dat <- fit$data
levels(dat$biogeographic_realm)
w <- prop.table(table(dat$biogeographic_realm))
cont_levels <- names(w)
n_realm <- as.data.frame(table(dat$biogeographic_realm))
names(n_realm) <- c("realm","n")

vars_darwinian <- c(
  "body_size", "taxonomic_effort", "taxonomic_activity", "watershed_area",
  "range_size", "elevation", "latitude", "discharge", "watertemp",
  "preserved_specimen", "sampling_effort","range_rarity", "population_density"
)

effect_counts_basin_darwinian <- map_dfr(
  vars_darwinian,
  ~ compute_realm_hr(
    v             = .x,
    draws         = draws,
    cont_levels   = cont_levels,
    baseline_realm = "Afrotropic",
    family = "gompertz"
  )
) %>%
  mutate(
    hr    = as.numeric(hr),
    realm = factor(realm, levels = cont_levels),
    variable  = factor(variable, levels = vars_darwinian)
  ) %>%
  group_by(variable, realm) %>%
  summarise(
    median   = median(hr, na.rm = TRUE),
    lower_95 = quantile(hr, 0.025, na.rm = TRUE),
    upper_95 = quantile(hr, 0.975, na.rm = TRUE),
    prob_pos = mean(hr > 1, na.rm = TRUE),   
    .groups  = "drop"
  ) %>%
  arrange(variable, realm) %>%
  mutate(
    effect_sign = case_when(
      lower_95 > 1 ~ "positive",
      upper_95 < 1 ~ "negative",
      TRUE         ~ "nonsignificant"
    )
  ) %>%
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

rm(fit,w,vars_darwinian,cont_levels,draws,dat,n_realm)
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
  taxonomic_effort     = "Taxonomic effort",
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
    labels = c("HR>1", "HR<1"),
    name   = "Direction"
  ) +
  labs(
    x = NULL,
    y = expression("Number of biogeographic realm with significant effects (" * italic(p) * " < 0.05)")
  ) +
  # annotate("text", x = 0.8, y = 6, label = "HR>1",colour = "#4C78A8", face = "bold",size = 3)+
  # annotate("text", x = 0.8, y = -6, label = "HR<1",colour = "#F58518", face = "bold",size = 3)+
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
# ---- helper: validate family and transform beta to HR ----
.beta_to_hr <- function(beta,
                        draws,
                        family = c("weibull", "gompertz"),
                        shape_col = "shape") {
  family <- match.arg(family)
  
  if (family == "gompertz") {
    log_hr <- beta
  } else if (family == "weibull") {
    if (!shape_col %in% names(draws)) {
      stop("No Weibull shape column found: ", shape_col)
    }
    shape  <- draws[[shape_col]]
    log_hr <- -shape * beta
  }
  
  list(
    log_hr = log_hr,
    hr     = exp(log_hr)
  )
}

# ---- helper: construct continent-specific beta ----
.compute_continent_beta <- function(v,
                                    draws,
                                    cont_levels,
                                    baseline_continent = "Africa",
                                    continent_var = "continent") {
  base_name <- paste0("b_", v)
  
  if (!base_name %in% names(draws)) {
    stop("No coefficient found for ", base_name)
  }
  
  base <- draws[[base_name]]
  
  purrr::map_dfr(cont_levels, function(ct) {
    if (ct == baseline_continent) {
      beta <- base
    } else {
      inter_name <- paste0("b_", v, ":", continent_var, ct)
      
      # fallback: try reversed interaction naming if needed
      inter_name_alt <- paste0("b_", continent_var, ct, ":", v)
      
      if (inter_name %in% names(draws)) {
        beta <- base + draws[[inter_name]]
      } else if (inter_name_alt %in% names(draws)) {
        beta <- base + draws[[inter_name_alt]]
      } else {
        beta <- base
      }
    }
    
    tibble::tibble(
      continent = ct,
      beta  = beta
    )
  })
}

# ---- continent-specific HR ----
compute_continent_hr <- function(v,
                                 draws,
                                 cont_levels,
                                 baseline_continent = "Africa",
                                 family = c("weibull", "gompertz"),
                                 shape_col = "shape",
                                 continent_var = "continent") {
  family <- match.arg(family)
  
  req_meta <- c(".chain", ".iteration", ".draw")
  has_meta <- req_meta[req_meta %in% names(draws)]
  
  if (length(has_meta) == 0) {
    stop("draws must contain at least one posterior meta column such as .draw")
  }
  
  meta_cols <- draws %>% dplyr::select(dplyr::all_of(has_meta))
  
  continent_beta <- .compute_continent_beta(
    v = v,
    draws = draws,
    cont_levels = cont_levels,
    baseline_continent = baseline_continent,
    continent_var = continent_var
  )
  
  trans <- .beta_to_hr(
    beta = continent_beta$beta,
    draws = draws[rep(seq_len(nrow(draws)), times = length(cont_levels)), , drop = FALSE],
    family = family,
    shape_col = shape_col
  )
  
  out <- continent_beta %>%
    dplyr::mutate(
      log_hr = trans$log_hr,
      hr     = trans$hr
    )
  
  dplyr::bind_cols(
    meta_cols[rep(seq_len(nrow(meta_cols)), times = length(cont_levels)), , drop = FALSE],
    out
  ) %>%
    dplyr::mutate(
      variable = v,
      family   = family
    ) %>%
    dplyr::relocate(variable, family, continent, .after = dplyr::last_col())
}

################################################################################
fit <- readRDS("output/model/country_linnaean_gompertz_full.rds")
draws <- as_draws_df(fit)
dat <- fit$data
levels(dat$continent)
w <- prop.table(table(dat$continent))
cont_levels <- names(w)
n_continent <- as.data.frame(table(dat$continent))
names(n_continent) <- c("continent","n")

vars_linnaean <- c(
  "body_size", "taxonomic_effort", "taxonomic_activity", "country_area",
  "range_size", "elevation", "latitude", "discharge", "watertemp",
  "preserved_specimen", "sampling_effort","sequencing_effort", "range_rarity", "population_density"
)

effect_counts_country_linnaean <- map_dfr(
  vars_linnaean,
  ~ compute_continent_hr(
    v             = .x,
    draws         = draws,
    cont_levels   = cont_levels,
    baseline_continent = "Africa",
    family = "gompertz"
  )
) %>%
  mutate(
    hr    = as.numeric(hr),
    continent = factor(continent, levels = cont_levels),
    variable  = factor(variable, levels = vars_linnaean)
  ) %>%
  group_by(variable, continent) %>%
  summarise(
    median   = median(hr, na.rm = TRUE),
    lower_95 = quantile(hr, 0.025, na.rm = TRUE),
    upper_95 = quantile(hr, 0.975, na.rm = TRUE),
    prob_pos = mean(hr > 1, na.rm = TRUE),   
    .groups  = "drop"
  ) %>%
  arrange(variable, continent) %>%
  mutate(
    effect_sign = case_when(
      lower_95 > 1 ~ "positive",
      upper_95 < 1 ~ "negative",
      TRUE         ~ "nonsignificant"
    )
  ) %>%
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

rm(fit,w,vars_linnaean,cont_levels,draws,dat,n_continent)
################################################################################
fit <- readRDS("output/model/country_wallacean_gompertz_full.rds")
draws <- as_draws_df(fit)
dat <- fit$data
levels(dat$continent)
w <- prop.table(table(dat$continent))
cont_levels <- names(w)
n_continent <- as.data.frame(table(dat$continent))
names(n_continent) <- c("continent","n")

vars_wallacean <- c(
  "body_size", "taxonomic_effort", "taxonomic_activity", "country_area",
  "range_size", "elevation", "latitude", "discharge", "watertemp",
  "preserved_specimen", "sequencing_effort", "range_rarity", "population_density"
)

effect_counts_country_wallacean <- map_dfr(
  vars_wallacean,
  ~ compute_continent_hr(
    v             = .x,
    draws         = draws,
    cont_levels   = cont_levels,
    baseline_continent = "Africa",
    family = "gompertz"
  )
) %>%
  mutate(
    hr    = as.numeric(hr),
    continent = factor(continent, levels = cont_levels),
    variable  = factor(variable, levels = vars_wallacean)
  ) %>%
  group_by(variable, continent) %>%
  summarise(
    median   = median(hr, na.rm = TRUE),
    lower_95 = quantile(hr, 0.025, na.rm = TRUE),
    upper_95 = quantile(hr, 0.975, na.rm = TRUE),
    prob_pos = mean(hr > 1, na.rm = TRUE),   
    .groups  = "drop"
  ) %>%
  arrange(variable, continent) %>%
  mutate(
    effect_sign = case_when(
      lower_95 > 1 ~ "positive",
      upper_95 < 1 ~ "negative",
      TRUE         ~ "nonsignificant"
    )
  ) %>%
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
rm(fit,w,vars_wallacean,cont_levels,draws,dat,n_continent)

################################################################################
fit <- readRDS("output/model/country_darwinian_gompertz_full.rds")
draws <- as_draws_df(fit)
dat <- fit$data
levels(dat$continent)
w <- prop.table(table(dat$continent))
cont_levels <- names(w)
n_continent <- as.data.frame(table(dat$continent))
names(n_continent) <- c("continent","n")

vars_darwinian <- c(
  "body_size", "taxonomic_effort", "taxonomic_activity", "country_area",
  "range_size", "elevation", "latitude", "discharge", "watertemp",
  "preserved_specimen", "sampling_effort","range_rarity", "population_density"
)

effect_counts_country_darwinian <- map_dfr(
  vars_darwinian,
  ~ compute_continent_hr(
    v             = .x,
    draws         = draws,
    cont_levels   = cont_levels,
    baseline_continent = "Africa",
    family = "gompertz"
  )
) %>%
  mutate(
    hr    = as.numeric(hr),
    continent = factor(continent, levels = cont_levels),
    variable  = factor(variable, levels = vars_darwinian)
  ) %>%
  group_by(variable, continent) %>%
  summarise(
    median   = median(hr, na.rm = TRUE),
    lower_95 = quantile(hr, 0.025, na.rm = TRUE),
    upper_95 = quantile(hr, 0.975, na.rm = TRUE),
    prob_pos = mean(hr > 1, na.rm = TRUE),   
    .groups  = "drop"
  ) %>%
  arrange(variable, continent) %>%
  mutate(
    effect_sign = case_when(
      lower_95 > 1 ~ "positive",
      upper_95 < 1 ~ "negative",
      TRUE         ~ "nonsignificant"
    )
  ) %>%
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
rm(fit,w,vars_darwinian,cont_levels,draws,dat,n_continent)

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
  taxonomic_effort    = "Taxonomic effort",
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
    labels = c("HR>1", "HR<1"),
    name   = "Direction"
  ) +
  labs(
    x = NULL,
    y = expression("Number of continent with significant effects (" * italic(p) * " < 0.05)")
  ) +
  # annotate("text", x = 0.8, y = 6, label = "HR>1",colour = "#4C78A8", face = "bold",size = 3)+
  # annotate("text", x = 0.8, y = -6, label = "HR<1",colour = "#F58518", face = "bold",size = 3)+
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
