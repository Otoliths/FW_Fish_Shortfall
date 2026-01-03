# ------------------------------------------------------------------------------
# Supplementary Figure S6
# Variable-wise and group-wise relative importance of predictors in survival models of species description across continents
# Variable importance for species discovery delay (Linnaean)
# Lognormal survival model + DALEX permutation importance
# ------------------------------------------------------------------------------
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)
library(brms)       # fitted lognormal survival model
library(posterior)  # ndraws()
library(DALEX)      # model-agnostic explanations
library(furrr)      # parallel purrr
library(progressr)
library(ggplot2)
library(ggh4x)
library(sf)
library(rnaturalearth)
library(cowplot)
packageVersion("progressr")
# ------------------------------------------------------------------------------
# 1. Load fitted Linnaean model and data
#    time  = species discovery delay (e.g. year_description - 1758)
#    event = 1 for described species, 0 for still-undescribed (censored)
# ------------------------------------------------------------------------------

fit_L <- readRDS("output/model/country_linnaean_all.rds")$weibull
dat_L <- fit_L$data   # or readRDS("output/stan_survdata_basin_linnaean.rds") if you prefer

continent <- levels(dat_L$continent)


draw_ids <- sample(
  x    = 1:posterior::ndraws(fit_L),
  size = 100
)

# ------------------------------------------------------------------------------
# 2. Predictor set and variable groups (Linnaean version)
# ------------------------------------------------------------------------------

predictors_L <- c(
  "body_size",
  "taxonmic_effort",
  "taxonomic_activity",
  "country_area",
  "range_size",
  "elevation",
  "latitude",
  "discharge",
  "watertemp",
  "preserved_specimen",
  "sequencing_effort",
  "sampling_effort",
  "range_rarity",
  "population_density",
  "continent"
)

var_groups_L <- list(
  Biology     = c("body_size", "preserved_specimen"),
  Geography   = c("range_size", "range_rarity", "latitude", "country_area"),
  Environment = c("discharge", "watertemp"),
  Access      = c("elevation", "population_density"),
  Activity    = c("taxonmic_effort", "taxonomic_activity",
                  "sampling_effort", "sequencing_effort")
)


# ------------------------------------------------------------------------------
# 3.Posterior-draw-specific prediction for expected discovery delay
# Weibull in brms (default parameterization):
#   mu    = scale, link usually log
#   shape = shape, link usually log
# ------------------------------------------------------------------------------

predict_L_delay_single_draw <- function(fit,
                                        newdata,
                                        draw_id,
                                        re_formula = NA) {
  
  # linear predictor for mu (scale)
  eta_mu <- posterior_linpred(
    fit,
    newdata    = newdata,
    dpar       = "mu",
    draw_ids   = draw_id,
    re_formula = re_formula,
    transform  = FALSE
  )
  eta_mu <- as.numeric(eta_mu)
  
  # linear predictor for shape
  eta_shape <- posterior_linpred(
    fit,
    newdata    = newdata,
    dpar       = "shape",
    draw_ids   = draw_id,
    re_formula = re_formula,
    transform  = FALSE
  )
  eta_shape <- as.numeric(eta_shape)
  
  # inverse links (assume log-links, which is brms default for weibull dpars)
  mu    <- exp(eta_mu)
  shape <- exp(eta_shape)
  
  # expected value for Weibull(scale=mu, shape=shape)
  mean_delay <- mu * gamma(1 + 1 / shape)
  
  as.numeric(mean_delay)
}

make_predict_fun_L_delay <- function(fit, draw_id) {
  force(fit); force(draw_id)
  function(model, newdata) {
    predict_L_delay_single_draw(
      fit        = fit,
      newdata    = newdata,
      draw_id    = draw_id,
      re_formula = NA
    )
  }
}


# ------------------------------------------------------------------------------
# 4. Pre-split data by continent (use only species with observed discovery times)
# ------------------------------------------------------------------------------

dat_L_by_continent <- dat_L %>%
  split(.$continent)

job_grid_L <- tidyr::expand_grid(
  continent   = continent,
  draw_id = draw_ids
)

# ------------------------------------------------------------------------------
# 5. Parallel plan
# ------------------------------------------------------------------------------

plan(multisession, workers = 60L)

# ------------------------------------------------------------------------------
# 6. Variable-wise importance for discovery delay
# ------------------------------------------------------------------------------
# > **Note:** The following code block is commented out by default due to its computational cost (~2 hours).  
# > To reproduce the results, please **uncomment** the code.

# with_progress({
#   p <- progressor(steps = nrow(job_grid_L))
#   results_L_varwise <- job_grid_L %>%
#     furrr::future_pmap_dfr(function(continent, draw_id) {
#       p(message = paste("continent:", continent, "| Draw:", draw_id))
#       dat_r <- dat_L_by_continent[[continent]]
#       
#       if (is.null(dat_r) || nrow(dat_r) == 0) {
#         return(tibble(
#           continent        = character(0),
#           draw_id      = integer(0),
#           variable     = character(0),
#           dropout_loss = numeric(0)
#         ))
#       }
#       
#       x_r <- dat_r[, predictors_L, drop = FALSE]
#       y_r <- dat_r$time   # observed species discovery delay
#       
#       pred_fun <- make_predict_fun_L_delay(
#         fit     = fit_L,
#         draw_id = draw_id
#       )
#       
#       expl <- DALEX::explain(
#         model            = fit_L,
#         data             = x_r,
#         y                = y_r,
#         label            = paste0("Linnaean_delay_", continent, "_draw_", draw_id),
#         predict_function = pred_fun,
#         verbose          = FALSE
#       )
#       
#       mp <- DALEX::model_parts(
#         explainer = expl,
#         type      = "variable_importance",
#         B         = 10
#       )
#       
#       mp_df <- if ("result" %in% names(mp)) mp$result else mp
#       
#       mp_df %>%
#         dplyr::filter(variable != "_full_model_") %>%
#         transmute(
#           continent        = continent,
#           draw_id      = draw_id,
#           variable     = variable,
#           dropout_loss = dropout_loss
#         )
#     })
# })
# saveRDS(results_L_varwise,"output/tables/results_L_varwise_country.rds")
# ------------------------------------------------------------------------------
# 7. Group-wise importance (Biology / Geography / Environment / Access / Activity)
# ------------------------------------------------------------------------------
# > **Note:** The following code block is commented out by default due to its computational cost (~2 hours).  
# > To reproduce the results, please **uncomment** the code.

# with_progress({
#   
#   p <- progressor(steps = nrow(job_grid_L))
#   results_L_groupwise <- job_grid_L %>%
#     furrr::future_pmap_dfr(function(continent, draw_id) {
#       p(message = paste("continent:", continent, "| Draw:", draw_id))
#       dat_r <- dat_L_by_continent[[continent]]
#       
#       if (is.null(dat_r) || nrow(dat_r) == 0) {
#         return(tibble(
#           continent        = character(0),
#           draw_id      = integer(0),
#           group        = character(0),
#           dropout_loss = numeric(0)
#         ))
#       }
#       
#       x_r <- dat_r[, predictors_L, drop = FALSE]
#       y_r <- dat_r$time
#       
#       pred_fun <- make_predict_fun_L_delay(
#         fit     = fit_L,
#         draw_id = draw_id
#       )
#       
#       expl <- DALEX::explain(
#         model            = fit_L,
#         data             = x_r,
#         y                = y_r,
#         label            = paste0("Linnaean_delay_", continent, "_draw_", draw_id),
#         predict_function = pred_fun,
#         verbose          = FALSE
#       )
#       
#       mp_group <- DALEX::model_parts(
#         explainer       = expl,
#         type            = "variable_importance",
#         B               = 10,
#         variable_groups = var_groups_L
#       )
#       
#       mp_group_df <- if ("result" %in% names(mp_group)) mp_group$result else mp_group
#       
#       mp_group_df %>%
#         filter(variable %in% names(var_groups_L)) %>%
#         transmute(
#           continent        = continent,
#           draw_id      = draw_id,
#           group        = variable,
#           dropout_loss = dropout_loss
#         )
#     })
# })
# 
# saveRDS(results_L_groupwise,"output/tables/results_L_groupwise_country.rds")
# summary_L_varwise    # variable-level importance for discovery delay
# summary_L_groupwise  # grouped importance

################################################################################
results_L_varwise <- readRDS("output/tables/results_L_varwise_country.rds")

summary_L_varwise <- results_L_varwise %>%
  filter(variable != "_baseline_") %>%
  filter(variable != "continent") %>%
  group_by(continent, variable) %>%
  mutate(
    loss_median = median(dropout_loss),
    loss_l95  = quantile(dropout_loss, 0.025),
    loss_u95  = quantile(dropout_loss, 0.975),
    .groups      = "drop"
  ) %>%
  arrange(continent, desc(loss_median)) %>%
  select(continent,variable,loss_median,loss_l95,loss_u95) %>%
  distinct()

summary_L_varwise_rel <- summary_L_varwise %>%
  group_by(continent) %>%
  mutate(
    rel_imp_median = loss_median / sum(loss_median, na.rm = TRUE),
    rel_imp_l95  = loss_l95  / sum(loss_l95,  na.rm = TRUE),
    rel_imp_u95  = loss_u95  / sum(loss_u95,  na.rm = TRUE)
  ) %>%
  ungroup()

rename_map <- c(
  body_size            = "Body size",
  preserved_specimen   = "N. preserved\nspecimens",
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



summary_L_varwise_rel <- summary_L_varwise_rel %>%
  mutate(
    variable = recode(variable, !!!rename_map),
    variable = factor(variable, levels = rev(unname(rename_map))),
    group = case_when(
      variable %in% c("Sequencing effort","Sampling effort",
                      "Taxonomic activity","Taxonomic effort") ~ "Activity",
      variable %in% c("Human density","Elevation") ~ "Access",
      variable %in% c("Water temperature","Streamflow") ~ "Environment",
      variable %in% c("Country area","Latitude",
                      "Range rarity","Range size") ~ "Geography",
      variable %in% c("N. preserved\nspecimens","Body size") ~ "Biology",
      TRUE ~ NA_character_
    ),
    group = factor(group,
                   levels = rev(c("Activity","Access",
                                  "Environment","Geography","Biology")))
  ) %>% as.data.frame()

summary_L_varwise_rel$continent <- factor(summary_L_varwise_rel$continent,levels = c( "North America",
                                                                                      "South America",
                                                                                      "Europe",
                                                                                      "Africa",
                                                                                      "Asia",
                                                                                      "Oceania"))

p1 <- ggplot(summary_L_varwise_rel, aes(x = rel_imp_median, y = variable, fill = continent)) +
  geom_col(alpha = 0.7, width = 0.6, position = position_dodge(0.6)) +
  ggh4x::facet_grid2(  
    group ~ .,
    scales = "free_y",
    space = "free",
    switch = "y",
    strip = strip_themed(
      background_y = elem_list_rect( colour = NA,
                                     fill = c(
                                       Biology     = "#6C8F5D80",
                                       Geography   = "#6D4D8080",
                                       Environment = "#A3473E80",
                                       Access      = "#4E6E8E80",
                                       Activity    = "#8A6A3F80"
                                     ) 
      )
    )
  ) +
  labs(y = "", x = "Variable-wise relative importance (%)") +
  scale_fill_manual(
    "Region",
    values = c(
      "Africa"   = "#D55E00",
      "Oceania" = "#0072B2",
      "Asia"  = "#009E73",
      "North America"     = "#F0E442",
      "South America"    = "#E69F00",
      "Europe"   = "#CC79A7"
    )
  ) +
  theme_minimal() +
  scale_x_continuous(expand = c(0, 0)) +
  theme(
    strip.background = element_blank(),  
    strip.placement = "outside",
    legend.key.size = unit(0.5,"cm"),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 8),
    strip.text = element_text(face = "bold", size = 8),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(colour = "black", size = 8)
  )

world <- ne_countries(scale = "small", returnclass = "sf") %>%
  st_wrap_dateline(options = c("WRAPDATELINE=YES"))
moll_proj <- st_crs("+proj=moll")
lat_points <- data.frame(lon = 0, lat = c(-60, 90)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
projected_points <- st_transform(lat_points, crs = moll_proj)
y_limits <- st_coordinates(projected_points)[, "Y"]

inland_country <- readRDS("input/raw/country.rds")
country <- read.csv("input/raw/country_list.csv")
inland_country <- inland_country %>% left_join(country[,c(1,5)], by = "iso3")

map <- ggplot() +
  ggrastr::rasterise(geom_sf(data = world, fill = "gray80", color = NA, size = 0.2),dpi = 300) +  
  ggrastr::rasterise(geom_sf(data = inland_country,aes(fill = continent), alpha = 0.7,
                             color = NA, linewidth = 0.5, show.legend = F),dpi=300) +  
  scale_fill_manual(values = c(
    "Africa"   = "#D55E00",
    "Oceania" = "#0072B2",
    "Asia"  = "#009E73",
    "North America"     = "#F0E442",
    "South America"    = "#E69F00",
    "Europe"   = "#CC79A7"
  ))+
  coord_sf(crs = "+proj=moll +lon_0=0", ylim = y_limits, expand = FALSE)+
  theme_void() +
  theme(plot.background = element_blank())

p11 <- ggdraw(p1)+
  draw_plot(map,x = 0.15, y = 0.01, scale = 0.3)


results_L_groupwise <- readRDS("output/tables/results_L_groupwise_country.rds")

summary_L_groupwise <- results_L_groupwise %>%
  group_by(continent, group) %>%
  mutate(
    loss_median = median(dropout_loss),
    loss_l95  = quantile(dropout_loss, 0.025),
    loss_u95  = quantile(dropout_loss, 0.975),
    .groups      = "drop"
  ) %>%
  arrange(continent, desc(loss_median)) %>%
  select(continent,group,loss_median) %>%
  distinct()

summary_L_groupwise_rel <- summary_L_groupwise %>%
  group_by(continent) %>%
  mutate(
    rel_imp_median = loss_median / sum(loss_median, na.rm = TRUE)
  ) %>%
  ungroup()

global <- summary_L_groupwise_rel %>%
  group_by(group) %>%
  summarise(
    rel_imp_median = mean(rel_imp_median),
    loss_median = mean(loss_median),
  ) %>%
  mutate(continent = "Global")

dt <-  dplyr::bind_rows(summary_L_groupwise_rel,global)

group_colors <- c(
  Biology     = "#6C8F5D80",
  Geography   = "#6D4D8080",
  Environment = "#A3473E80",
  Access      = "#4E6E8E80",
  Activity    = "#8A6A3F80"
)


dt$continent <- factor(dt$continent,levels = c( "North America",
                                                "South America",
                                                "Europe",
                                                "Africa",
                                                "Asia",
                                                "Oceania","Global"))

p2 <- ggplot(dt, aes(x = "", y = rel_imp_median, fill = group, group = group)) +
  geom_bar(width = 1, stat = "identity", color = "white", linewidth = 0.5,
           alpha = 0.8, show.legend = F) +
  coord_polar(theta = "y",start=0) + 
  scale_fill_manual("Group",values = group_colors)+
  facet_grid(.~ continent)+
  theme_minimal()+
  labs(y= "Group-wise relative Importance (%)", x= "")+
  theme(axis.text = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(0,0,0,0),
        legend.position = "left",
        strip.text = element_text(size = 9),
        axis.title = element_text(face = "bold", size = 10))
plot_grid(
  p11,
  p2,
  ncol = 1,
  rel_heights = c(0.8, 0.2)
)


ggsave("figures/supplement/Figure_S6.png",dpi = 300, units = "cm", width = 20, height = 20)

