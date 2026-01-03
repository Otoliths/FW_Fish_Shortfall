# ------------------------------------------------------------------------------
# Supplementary Figure S7
# Variable importance for occurrence record delay (Wallacean)
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
# 1. Predictor set and variable groups (Wallacean version)
# ------------------------------------------------------------------------------

predictors_W <- c(
  "body_size",
  "taxonmic_effort",
  "taxonomic_activity",
  "watershed_area",
  "range_size",
  "elevation",
  "latitude",
  "discharge",
  "watertemp",
  "preserved_specimen",
  "sequencing_effort",
  "range_rarity",
  "population_density",
  "biogeographic_realm"
)

var_groups_W <- list(
  Biology     = c("body_size", "preserved_specimen"),
  Geography   = c("range_size", "range_rarity", "latitude", "watershed_area"),
  Environment = c("discharge", "watertemp"),
  Access      = c("elevation", "population_density"),
  Activity    = c("taxonmic_effort", "taxonomic_activity", "sequencing_effort")
)

# ------------------------------------------------------------------------------
# 2. Load fitted Wallacean model and survival data
#    time  = occurrence record delay
#    event = 1 for species with at least one georeferenced record
#          = 0 for right-censored species (no record by cutoff year)
# ------------------------------------------------------------------------------

fit_W <- readRDS("output/model/basin_wallacean_all.rds")$lognormal
dat_W <- readRDS("output/stan_survdata_basin_wallacean.rds")

realms <- levels(dat_W$biogeographic_realm)


draw_ids <- sample(
  x    = 1:posterior::ndraws(fit_W),
  size = 100
)

# ------------------------------------------------------------------------------
# 3. Posterior-draw-specific prediction for expected delay (lognormal mean)
#    Lognormal parameterisation in brms:
#      mu    = meanlog
#      sigma = sdlog
#    E[T] = exp(mu + 0.5 * sigma^2)
# ------------------------------------------------------------------------------

predict_W_delay_single_draw <- function(fit,
                                        newdata,
                                        draw_id,
                                        re_formula = NA) {
  # linear predictors for mu (meanlog) and sigma (sdlog)
  mu <- posterior_linpred(
    fit,
    newdata    = newdata,
    dpar       = "mu",
    draw_ids   = draw_id,
    re_formula = re_formula,
    transform  = FALSE
  )
  
  sigma <- posterior_linpred(
    fit,
    newdata    = newdata,
    dpar       = "sigma",
    draw_ids   = draw_id,
    re_formula = re_formula,
    transform  = FALSE
  )
  
  mu    <- as.numeric(mu)
  sigma <- as.numeric(sigma)
  
  # expected occurrence delay under the lognormal distribution
  mean_delay <- exp(mu + 0.5 * sigma^2)
  
  as.numeric(mean_delay)
}

make_predict_fun_W_delay <- function(fit, draw_id) {
  force(fit); force(draw_id)
  function(model, newdata) {
    predict_W_delay_single_draw(
      fit        = fit,
      newdata    = newdata,
      draw_id    = draw_id,
      re_formula = NA
    )
  }
}

# ------------------------------------------------------------------------------
# 4. Pre-split data by realm (only species with observed occurrences, event == 1)
# ------------------------------------------------------------------------------

dat_W_by_realm <- dat_W %>%
  filter(event == 1) %>%           # use only species with realised occurrence delay
  split(.$biogeographic_realm)

# Grid of (realm, draw_id) jobs for parallel execution
job_grid_W <- tidyr::expand_grid(
  realm   = realms,
  draw_id = draw_ids
)

# ------------------------------------------------------------------------------
# 5. Parallel plan (adjust workers according to your machine)
# ------------------------------------------------------------------------------

plan(multisession, workers = 72L)

# ------------------------------------------------------------------------------
# 6. Variable-wise importance for occurrence record delay
# ------------------------------------------------------------------------------
# > **Note:** The following code block is commented out by default due to its computational cost (~4 hours).  
# > To reproduce the results, please **uncomment** the code.

# with_progress({
#   
#   p <- progressor(steps = nrow(job_grid_W))
#   results_W_varwise <- job_grid_W %>%
#     furrr::future_pmap_dfr(function(realm, draw_id) {
#       p(message = paste("Realm:", realm, "| Draw:", draw_id))
#       dat_r <- dat_W_by_realm[[realm]]
#       
#       # if this realm has no observed events, return empty tibble
#       if (is.null(dat_r) || nrow(dat_r) == 0) {
#         return(tibble(
#           realm        = character(0),
#           draw_id      = integer(0),
#           variable     = character(0),
#           dropout_loss = numeric(0)
#         ))
#       }
#       
#       x_r <- dat_r[, predictors_W, drop = FALSE]
#       y_r <- dat_r$time   # observed occurrence delay
#       
#       pred_fun <- make_predict_fun_W_delay(
#         fit     = fit_W,
#         draw_id = draw_id
#       )
#       
#       expl <- DALEX::explain(
#         model            = fit_W,
#         data             = x_r,
#         y                = y_r,
#         label            = paste0("Wallacean_delay_", realm, "_draw_", draw_id),
#         predict_function = pred_fun,
#         verbose          = TRUE
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
#           realm        = realm,
#           draw_id      = draw_id,
#           variable     = variable,
#           dropout_loss = dropout_loss
#         )
#     })
# })
# 
# saveRDS(results_W_varwise,"output/tables/results_W_varwise.rds")
# ------------------------------------------------------------------------------
# 7. Group-wise importance (Biology / Geography / Environment / Access / Activity)
# ------------------------------------------------------------------------------
# > **Note:** The following code block is commented out by default due to its computational cost (~4 hours).  
# > To reproduce the results, please **uncomment** the code.

# with_progress({
#   
#   p <- progressor(steps = nrow(job_grid_W))
#   results_W_groupwise <- job_grid_W %>%
#     furrr::future_pmap_dfr(function(realm, draw_id) {
#       p(message = paste("Realm:", realm, "| Draw:", draw_id))
#       dat_r <- dat_W_by_realm[[realm]]
#       
#       if (is.null(dat_r) || nrow(dat_r) == 0) {
#         return(tibble(
#           realm        = character(0),
#           draw_id      = integer(0),
#           group        = character(0),
#           dropout_loss = numeric(0)
#         ))
#       }
#       
#       x_r <- dat_r[, predictors_W, drop = FALSE]
#       y_r <- dat_r$time
#       
#       pred_fun <- make_predict_fun_W_delay(
#         fit     = fit_W,
#         draw_id = draw_id
#       )
#       
#       expl <- DALEX::explain(
#         model            = fit_W,
#         data             = x_r,
#         y                = y_r,
#         label            = paste0("Wallacean_delay_", realm, "_draw_", draw_id),
#         predict_function = pred_fun,
#         verbose          = TRUE
#       )
#       
#       mp_group <- DALEX::model_parts(
#         explainer       = expl,
#         type            = "variable_importance",
#         B               = 10,
#         variable_groups = var_groups_W
#       )
#       
#       mp_group_df <- if ("result" %in% names(mp_group)) mp_group$result else mp_group
#       
#       mp_group_df %>%
#         filter(variable %in% names(var_groups_W)) %>%
#         transmute(
#           realm        = realm,
#           draw_id      = draw_id,
#           group        = variable,
#           dropout_loss = dropout_loss
#         )
#     })
# })
# 
# saveRDS(results_W_groupwise,"output/tables/results_W_groupwise.rds")
# summary_W_varwise    # variable-level importance for occurrence record delay
# summary_W_groupwise  # grouped importance (Biology / Geography / Environment / Access / Activity)
################################################################################
results_W_varwise <- readRDS("output/tables/results_W_varwise.rds")

summary_W_varwise <- results_W_varwise %>%
  filter(variable != "_baseline_") %>%
  filter(variable != "biogeographic_realm") %>%
  group_by(realm, variable) %>%
  mutate(
    loss_median = median(dropout_loss),
    loss_l95  = quantile(dropout_loss, 0.025),
    loss_u95  = quantile(dropout_loss, 0.975),
    .groups      = "drop"
  ) %>%
  arrange(realm, desc(loss_median)) %>%
  select(realm,variable,loss_median,loss_l95,loss_u95) %>%
  distinct()

summary_W_varwise_rel <- summary_W_varwise %>%
  group_by(realm) %>%
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
  watershed_area       = "Watershed area",
  discharge            = "Streamflow",
  watertemp            = "Water temperature",
  elevation            = "Elevation",
  population_density   = "Human density",
  taxonmic_effort      = "Taxonomic effort",
  taxonomic_activity   = "Taxonomic activity",
  #sampling_effort      = "Sampling effort",
  sequencing_effort    = "Sequencing effort"
)



summary_W_varwise_rel <- summary_W_varwise_rel %>%
  mutate(
    realm     = recode(
      realm,
      "Australasia" = "Australasian",
      "Oceania"     = "Oceanian"
    ),
    variable = recode(variable, !!!rename_map),
    variable = factor(variable, levels = rev(unname(rename_map))),
    group = case_when(
      variable %in% c("Sequencing effort",#"Sampling effort",
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

summary_W_varwise_rel$realm <- factor(summary_W_varwise_rel$realm,levels = c( "Nearctic","Neotropic","Palearctic",
                                                                              "Afrotropic","Indomalayan","Australasian","Oceanian"))

p1 <- ggplot(summary_W_varwise_rel, aes(x = rel_imp_median, y = variable, fill = realm)) +
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
      "Afrotropic" = "#D55E00",
      "Australasian" = "#0072B2",
      "Indomalayan" = "#009E73",
      "Nearctic" = "#F0E442",
      "Neotropic" = "#E69F00",
      "Oceanian" = "#56B4E9",
      "Palearctic" = "#CC79A7"
    )
  ) +
  theme_minimal() +
  scale_x_continuous(expand = c(0, 0)) +
  theme(
    strip.background = element_blank(),  # 必须禁用默认背景
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
  draw_plot(map,x = 0.15, y = 0.01, scale = 0.3)


results_W_groupwise <- readRDS("output/tables/results_W_groupwise.rds")

summary_W_groupwise <- results_W_groupwise %>%
  group_by(realm, group) %>%
  mutate(
    loss_median = median(dropout_loss),
    loss_l95  = quantile(dropout_loss, 0.025),
    loss_u95  = quantile(dropout_loss, 0.975),
    .groups      = "drop"
  ) %>%
  arrange(realm, desc(loss_median)) %>%
  select(realm,group,loss_median) %>%
  distinct()

summary_W_groupwise_rel <- summary_W_groupwise %>%
  group_by(realm) %>%
  mutate(
    rel_imp_median = loss_median / sum(loss_median, na.rm = TRUE)
  ) %>%
  ungroup()

global <- summary_W_groupwise_rel %>%
  group_by(group) %>%
  summarise(
    rel_imp_median = mean(rel_imp_median),
    loss_median = mean(loss_median),
  ) %>%
  mutate(realm = "Global")

dt <-  dplyr::bind_rows(summary_W_groupwise_rel,global)

group_colors <- c(
  Biology     = "#6C8F5D80",
  Geography   = "#6D4D8080",
  Environment = "#A3473E80",
  Access      = "#4E6E8E80",
  Activity    = "#8A6A3F80"
)

dt <- dt %>% mutate(realm     = recode(
  realm,
  "Australasia" = "Australasian",
  "Oceania"     = "Oceanian"
))


dt$realm <- factor(dt$realm,levels = c( "Nearctic","Neotropic","Palearctic",
                                        "Afrotropic","Indomalayan","Australasian","Oceanian","Global"))

p2 <- ggplot(dt, aes(x = "", y = rel_imp_median, fill = group, group = group)) +
  geom_bar(width = 1, stat = "identity", color = "white", linewidth = 0.5,
           alpha = 0.8, show.legend = F) +
  coord_polar(theta = "y",start=0) + 
  scale_fill_manual("Group",values = group_colors)+
  facet_grid(.~ realm)+
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


ggsave("figures/supplement/Figure_S7.png",dpi = 300, units = "cm", width = 20, height = 20)
