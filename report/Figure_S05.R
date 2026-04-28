# ------------------------------------------------------------------------------
# Supplementary Figure S5
# Variable-wise and group-wise relative importance of predictors in survival models of species description across biogeographic realms
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

fit_L <- readRDS("output/model/basin_linnaean_weibull_full.rds")
dat_L <- fit_L$data   

realms <- levels(dat_L$biogeographic_realm)


draw_ids <- sample(
  x    = 1:posterior::ndraws(fit_L),
  size = 100
)

# ------------------------------------------------------------------------------
# 2. Predictor set and variable groups (Linnaean version)
# ------------------------------------------------------------------------------

predictors_L <- c(
  "body_size",
  "taxonomic_effort",
  "taxonomic_activity",
  "watershed_area",
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
  "biogeographic_realm"
)

var_groups_L <- list(
  Biology     = c("body_size", "preserved_specimen"),
  Geography   = c("range_size", "range_rarity", "latitude", "watershed_area"),
  Environment = c("discharge", "watertemp"),
  Access      = c("elevation", "population_density"),
  Activity    = c("taxonomic_effort", "taxonomic_activity",
                  "sampling_effort", "sequencing_effort")
)


# ------------------------------------------------------------------------------
# 3. Posterior-draw-specific prediction for expected discovery delay
#    Weibull parameterisation in brms:
#      mu    = mean (E[T])
#      shape = shape parameter
#
#    IMPORTANT:
#      In brms Weibull, mu is already the expected value:
#        E[T] = mu
# ------------------------------------------------------------------------------

predict_L_delay_single_draw <- function(fit,
                                        newdata,
                                        draw_id,
                                        re_formula = NA) {
  mu <- posterior_linpred(
    fit,
    newdata    = newdata,
    dpar       = "mu",
    draw_ids   = draw_id,
    re_formula = re_formula,
    transform  = TRUE
  )
  
  # Extracted for completeness; not needed for the mean in brms Weibull
  shape <- posterior_linpred(
    fit,
    newdata    = newdata,
    dpar       = "shape",
    draw_ids   = draw_id,
    re_formula = re_formula,
    transform  = TRUE
  )
  
  mu    <- as.numeric(mu)
  shape <- as.numeric(shape)
  
  mean_delay <- mu
  
  as.numeric(mean_delay)
}

# Use `model` passed by DALEX rather than capturing `fit` in the closure
make_predict_fun_L_delay <- function(draw_id) {
  force(draw_id)
  function(model, newdata) {
    predict_L_delay_single_draw(
      fit        = model,
      newdata    = newdata,
      draw_id    = draw_id,
      re_formula = NA
    )
  }
}

# ------------------------------------------------------------------------------
# 4. Pre-split data by realm
# ------------------------------------------------------------------------------

dat_L_by_realm <- dat_L %>%
  split(.$biogeographic_realm)

job_grid_L <- tidyr::expand_grid(
  realm   = realms,
  draw_id = draw_ids
)

# ------------------------------------------------------------------------------
# 5. Parallel plan
# ------------------------------------------------------------------------------

plan(multisession, workers = 8L)

# ------------------------------------------------------------------------------
# 5a. MSE loss function
# ------------------------------------------------------------------------------

mse_loss <- function(observed, predicted) {
  ok <- is.finite(observed) & is.finite(predicted)
  if (!any(ok)) return(NA_real_)
  mean((observed[ok] - predicted[ok])^2)
}
attr(mse_loss, "loss_name") <- "MSE"

# ------------------------------------------------------------------------------
# 6. Variable-wise importance for discovery delay
#    Reported as increase in mean squared error (ΔMSE)
# ------------------------------------------------------------------------------

# with_progress({
#   p <- progressor(steps = nrow(job_grid_L))
#   
#   results_L_varwise <- job_grid_L %>%
#     furrr::future_pmap_dfr(function(realm, draw_id) {
#       p(message = paste("Realm:", realm, "| Draw:", draw_id))
#       
#       dat_r <- dat_L_by_realm[[realm]]
#       
#       if (is.null(dat_r) || nrow(dat_r) == 0) {
#         return(tibble(
#           realm     = character(0),
#           draw_id   = integer(0),
#           variable  = character(0),
#           delta_mse = numeric(0)
#         ))
#       }
#       
#       x_r <- dat_r[, predictors_L, drop = FALSE]
#       y_r <- dat_r$time   # observed species discovery delay
#       
#       pred_fun <- make_predict_fun_L_delay(
#         draw_id = draw_id
#       )
#       
#       expl <- DALEX::explain(
#         model            = fit_L,
#         data             = x_r,
#         y                = y_r,
#         label            = paste0("Linnaean_delay_", realm, "_draw_", draw_id),
#         predict_function = pred_fun,
#         verbose          = FALSE
#       )
#       
#       mp <- DALEX::model_parts(
#         explainer     = expl,
#         loss_function = mse_loss,
#         type          = "difference",
#         B             = 10
#       )
#       
#       mp_df <- if ("result" %in% names(mp)) mp$result else mp
#       
#       mp_df %>%
#         dplyr::filter(variable != "_full_model_") %>%
#         dplyr::transmute(
#           realm     = realm,
#           draw_id   = draw_id,
#           variable  = variable,
#           delta_mse = dropout_loss
#         )
#     })
# })
# 
# saveRDS(results_L_varwise, "output/tables/results_L_varwise.rds")

# ------------------------------------------------------------------------------
# 7. Group-wise importance (Biology / Geography / Environment / Access / Activity)
#    Reported as increase in mean squared error (ΔMSE)
# ------------------------------------------------------------------------------

# with_progress({
#   
#   p <- progressor(steps = nrow(job_grid_L))
#   
#   results_L_groupwise <- job_grid_L %>%
#     furrr::future_pmap_dfr(function(realm, draw_id) {
#       p(message = paste("Realm:", realm, "| Draw:", draw_id))
#       
#       dat_r <- dat_L_by_realm[[realm]]
#       
#       if (is.null(dat_r) || nrow(dat_r) == 0) {
#         return(tibble(
#           realm     = character(0),
#           draw_id   = integer(0),
#           group     = character(0),
#           delta_mse = numeric(0)
#         ))
#       }
#       
#       x_r <- dat_r[, predictors_L, drop = FALSE]
#       y_r <- dat_r$time
#       
#       pred_fun <- make_predict_fun_L_delay(
#         draw_id = draw_id
#       )
#       
#       expl <- DALEX::explain(
#         model            = fit_L,
#         data             = x_r,
#         y                = y_r,
#         label            = paste0("Linnaean_delay_", realm, "_draw_", draw_id),
#         predict_function = pred_fun,
#         verbose          = FALSE
#       )
#       
#       mp_group <- DALEX::model_parts(
#         explainer       = expl,
#         loss_function   = mse_loss,
#         type            = "difference",
#         B               = 10,
#         variable_groups = var_groups_L
#       )
#       
#       mp_group_df <- if ("result" %in% names(mp_group)) mp_group$result else mp_group
#       
#       mp_group_df %>%
#         dplyr::filter(variable %in% names(var_groups_L)) %>%
#         dplyr::transmute(
#           realm     = realm,
#           draw_id   = draw_id,
#           group     = variable,
#           delta_mse = dropout_loss
#         )
#     })
# })
# 
# saveRDS(results_L_groupwise, "output/tables/results_L_groupwise.rds")
# summary_L_varwise    # variable-level importance for discovery delay
# summary_L_groupwise  # grouped importance
################################################################################
results_L_varwise <- readRDS("output/tables/results_L_varwise.rds")

final_L_importance <- results_L_varwise %>%
  filter(!variable %in% c("_baseline_", "_full_model_", "biogeographic_realm")) %>%
  group_by(realm, draw_id) %>%
  mutate(
    delta_pos = pmax(delta_mse, 0),
    denom_pos = sum(delta_pos, na.rm = TRUE),
    rel_imp_pos = if_else(denom_pos > 0, delta_pos / denom_pos, NA_real_)
  ) %>%
  ungroup() %>%
  group_by(realm, variable) %>%
  summarise(
    delta_mse_median   = median(delta_mse, na.rm = TRUE),
    delta_mse_l95      = quantile(delta_mse, 0.025, na.rm = TRUE),
    delta_mse_u95      = quantile(delta_mse, 0.975, na.rm = TRUE),
    rel_imp_pos_median = median(rel_imp_pos, na.rm = TRUE),
    rel_imp_pos_l95    = quantile(rel_imp_pos, 0.025, na.rm = TRUE),
    rel_imp_pos_u95    = quantile(rel_imp_pos, 0.975, na.rm = TRUE),
    p_negative         = mean(delta_mse < 0, na.rm = TRUE),
    p_positive         = mean(delta_mse > 0, na.rm = TRUE),
    sign_class = case_when(
      delta_mse_l95 > 0 ~ "robust_positive",
      delta_mse_u95 < 0 ~ "robust_negative",
      TRUE              ~ "sign_uncertain"
    ),
    .groups = "drop"
  ) %>%
  arrange(realm, desc(delta_mse_median))

rename_map <- c(
  body_size            = "Body size",
  preserved_specimen   = "N. preserved\nspecimens",
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



summary_L_varwise <- final_L_importance %>%
  mutate(
    realm     = recode(
      realm,
      "Australasia" = "Australasian",
      "Oceania"     = "Oceanian"
    ),
    variable = recode(variable, !!!rename_map),
    variable = factor(variable, levels = rev(unname(rename_map))),
    group = case_when(
      variable %in% c("Sequencing effort","Sampling effort",
                      "Taxonomic activity","Taxonomic effort") ~ "Activity",
      variable %in% c("Human density","Elevation") ~ "Access",
      variable %in% c("Water temperature","Streamflow") ~ "Environment",
      variable %in% c("Watershed area","Latitude",
                      "Range rarity","Range size") ~ "Geography",
      variable %in% c("N. preserved\nspecimens","Body size") ~ "Biology",
      TRUE ~ NA_character_
    ),
    group = factor(group,
                   levels = rev(c("Activity","Access",
                                  "Environment","Geography","Biology")))
  ) %>% as.data.frame()

summary_L_varwise$realm <- factor(summary_L_varwise$realm,levels = c( "Nearctic","Neotropic","Palearctic",
                                                                      "Afrotropic","Indomalayan","Australasian","Oceanian"))

p1 <- ggplot(data = summary_L_varwise, aes(x = delta_mse_median, y = variable, colour = realm)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = delta_mse_l95, xmax = delta_mse_u95), width = 0.18, linewidth = 0.5, position = position_dodge(width = 0.6)) +
  geom_point(size = 2.2,position = position_dodge(width = 0.6)) +
  ggh4x::facet_grid2(
    group ~ .,
    scales = "free_y",
    space = "free",
    switch = "y",
    strip = strip_themed(
      background_y = elem_list_rect(
        colour = NA,
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
  labs(
    y = "",
    x = expression(Delta * MSE)
  ) +
  scale_colour_manual(
    "Realm",
    values = c(
      "Afrotropic" = "#D55E00",
      "Australasian" = "#0072B2",
      "Indomalayan" = "#009E73",
      "Nearctic" = "#F0E442",
      "Neotropic" = "#E69F00",
      "Oceanian" = "#56B4E9",
      "Palearctic" = "#CC79A7"
    )
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05)),limits = c(-300,4000)) +
  theme_minimal() +
  theme(
    strip.background = element_blank(),
    strip.placement = "outside",
    legend.key.size = unit(0.5, "cm"),
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 8),
    strip.text = element_text(face = "bold", size = 8),
    axis.title = element_text(face = "bold", size = 10),
    axis.text = element_text(colour = "black", size = 8),
    panel.grid.minor = element_blank()
  )


world <- ne_countries(scale = "small", returnclass = "sf") %>%
  st_wrap_dateline(options = c("WRAPDATELINE=YES"))
moll_proj <- st_crs("+proj=moll")
lat_points <- data.frame(lon = 0, lat = c(-60, 90)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
projected_points <- st_transform(lat_points, crs = moll_proj)
y_limits <- st_coordinates(projected_points)[, "Y"]
inland_basins <- readRDS("input/raw/basin/basin_sf_v1.rds")

inland_basins <- inland_basins %>%
  mutate(
    biogeographic_realm = dplyr::case_when(
      biogeographic_realm == "Australasia" ~ "Australasian",
      biogeographic_realm == "Oceania"     ~ "Oceanian",
      TRUE ~ biogeographic_realm
    )
  )

map <- ggplot() +
  ggrastr::rasterise(geom_sf(data = world, fill = "gray80", color = NA, size = 0.2),dpi = 300) +  
  ggrastr::rasterise(geom_sf(data = inland_basins,aes(fill = biogeographic_realm), alpha = 0.7,
                             color = NA, linewidth = 0.5, show.legend = F),dpi=300) +  
  scale_fill_manual(values = c(
    "Afrotropic" = "#D55E00",
    "Australasian" = "#0072B2",
    "Indomalayan" = "#009E73",
    "Nearctic" = "#F0E442",
    "Neotropic" = "#E69F00",
    "Oceanian" = "#56B4E9",
    "Palearctic" = "#CC79A7"
  ))+
  coord_sf(crs = "+proj=moll +lon_0=0", ylim = y_limits, expand = FALSE)+
  theme_void() +
  theme(plot.background = element_blank())

p11 <- ggdraw(p1)+
  draw_plot(map,x = 0.2, y = 0.4, scale = 0.3)


results_L_groupwise <- readRDS("output/tables/results_L_groupwise.rds")

results_L_groupwise <- results_L_groupwise %>% mutate(realm = recode(
  realm,
  "Australasia" = "Australasian",
  "Oceania"     = "Oceanian"
))

results_L_groupwise_clean <- results_L_groupwise %>%
  group_by(realm, draw_id, group) %>%
  summarise(
    delta_mse = sum(delta_mse, na.rm = TRUE),
    .groups = "drop"
  )

groupwise_draw_rel <- results_L_groupwise_clean %>%
  group_by(realm, draw_id) %>%
  mutate(
    delta_pos = pmax(delta_mse, 0),
    denom_pos = sum(delta_pos, na.rm = TRUE),
    rel_imp = if_else(denom_pos > 0, delta_pos / denom_pos, NA_real_)
  ) %>%
  ungroup()

summary_L_groupwise_rel <- groupwise_draw_rel %>%
  group_by(realm, group) %>%
  summarise(
    rel_imp_median = median(rel_imp, na.rm = TRUE),
    rel_imp_l95    = quantile(rel_imp, 0.025, na.rm = TRUE),
    rel_imp_u95    = quantile(rel_imp, 0.975, na.rm = TRUE),
    delta_mse_median = median(delta_mse, na.rm = TRUE),
    p_negative = mean(delta_mse < 0, na.rm = TRUE),
    .groups = "drop"
  )

global <- summary_L_groupwise_rel %>%
  group_by(group) %>%
  summarise(
    rel_imp_median = mean(rel_imp_median, na.rm = TRUE),
    rel_imp_l95    = mean(rel_imp_l95, na.rm = TRUE),
    rel_imp_u95    = mean(rel_imp_u95, na.rm = TRUE),
    delta_mse_median = mean(delta_mse_median, na.rm = TRUE),
    p_negative = mean(p_negative, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(realm = "Global")

dt <- bind_rows(summary_L_groupwise_rel, global) %>%
  mutate(
    rel_imp_pct = rel_imp_median * 100
  )

dt <- dt %>%
  group_by(realm) %>%
  mutate(
    rel_imp_pct = 100 * rel_imp_pct / sum(rel_imp_pct, na.rm = TRUE)
  ) %>%
  ungroup()

group_colors <- c(
  Biology     = "#6C8F5D80",
  Geography   = "#6D4D8080",
  Environment = "#A3473E80",
  Access      = "#4E6E8E80",
  Activity    = "#8A6A3F80"
)


dt$realm <- factor(dt$realm,levels = c( "Nearctic","Neotropic","Palearctic",
                                        "Afrotropic","Indomalayan","Australasian","Oceanian","Global"))


dt_plot <- dt %>%
  mutate(
    group = factor(
      group,
      levels = c("Biology", "Geography", "Environment", "Access", "Activity")
    ),
    label = ifelse(rel_imp_pct >= 5, paste0(round(rel_imp_pct), "%"), "")
  )

p2 <- ggplot(dt_plot, aes(x = 1, y = rel_imp_pct, fill = group)) +
  geom_col(
    width = 1,
    color = "white",
    linewidth = 0.5,
    alpha = 0.85,
    show.legend = FALSE
  ) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),
    size = 2.2,
    color = "black"
  ) +
  coord_polar(theta = "y", start = 0) +
  scale_fill_manual(values = group_colors) +
  facet_grid(. ~ realm) +
  labs(
    x = NULL,
    y = "Group-wise relative importance (%)"
  ) +
  theme_void() +
  theme(
    legend.position = "left",
    plot.margin = margin(0, 0, 0, 0),
    strip.text = element_text(size = 9),
    axis.title = element_text(face = "bold", size = 10)
  )

plot_grid(
  p11,
  p2,
  ncol = 1,
  rel_heights = c(0.8, 0.2)
)


ggsave("figures/supplement/Figure_S5.png",dpi = 300, units = "cm", width = 20, height = 20)

