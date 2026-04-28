# ------------------------------------------------------------------------------
# Supplementary Figure S10 (continents)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Darwinian shortfall: variable importance for sequence acquisition delay
# Gamma survival model + DALEX permutation importance
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
# 1. Predictor set and variable groups (Darwinian version)
# ------------------------------------------------------------------------------

predictors_D <- c(
  "body_size",
  "taxonomic_effort",
  "taxonomic_activity",
  "country_area",
  "range_size",
  "elevation",
  "latitude",
  "discharge",
  "watertemp",
  "preserved_specimen",
  "sampling_effort",
  "range_rarity",
  "population_density",
  "continent"
)

var_groups_D <- list(
  Biology     = c("body_size", "preserved_specimen"),
  Geography   = c("range_size", "range_rarity", "latitude", "country_area"),
  Environment = c("discharge", "watertemp"),
  Access      = c("elevation", "population_density"),
  Activity    = c("taxonomic_effort", "taxonomic_activity", "sampling_effort")
)

# ------------------------------------------------------------------------------
# 2. Load fitted Weibull model and survival data
#    time  = sequence acquisition delay
#    event = 1 for sequenced species, 0 for right-censored
# ------------------------------------------------------------------------------

fit_D  <- readRDS("output/model/country_darwinian_gompertz_full.rds")
dat_D  <- readRDS("output/stan_survdata_country_darwinian_single.rds")

continent <- levels(dat_D$continent)

draw_ids <- sample(
  x    = 1:posterior::ndraws(fit_D),
  size = 100
)

# ------------------------------------------------------------------------------
# 3. Gompertz survival and mean
#    brms Gompertz parameterisation:
#      mu    = rate parameter, usually with log link
#      gamma = shape parameter, usually with identity link
#
#    Survival:
#      S(t) = exp(-(mu / gamma) * (exp(gamma * t) - 1)),  gamma != 0
#      S(t) = exp(-mu * t),                               gamma -> 0
#
#    Mean:
#      E[T] = integral_0^Inf S(t) dt
# ------------------------------------------------------------------------------

gompertz_surv <- function(t, mu, gamma, eps = 1e-8) {
  out <- numeric(length(t))
  near0 <- abs(gamma) < eps
  
  if (any(near0)) {
    out[near0] <- exp(-mu[near0] * t[near0])
  }
  
  if (any(!near0)) {
    g  <- gamma[!near0]
    m  <- mu[!near0]
    tt <- t[!near0]
    out[!near0] <- exp(-(m / g) * (exp(g * tt) - 1))
  }
  
  out
}

gompertz_mean <- function(mu, gamma, eps = 1e-8,
                          rel.tol = 1e-8,
                          subdivisions = 1000L,
                          tail_tol = 1e-10,
                          max_upper = 1e6) {
  stopifnot(length(mu) == length(gamma))
  
  out <- rep(NA_real_, length(mu))
  
  for (i in seq_along(mu)) {
    m <- mu[i]
    g <- gamma[i]
    
    if (!is.finite(m) || !is.finite(g) || m <= 0) {
      out[i] <- NA_real_
      next
    }
    
    # Under this parameterisation, if gamma < 0, the mean is not finite
    if (g < -eps) {
      out[i] <- Inf
      next
    }
    
    # Exponential limit
    if (abs(g) < eps) {
      out[i] <- 1 / m
      next
    }
    
    surv_i <- function(t) {
      exp(-(m / g) * (exp(g * t) - 1))
    }
    
    upper <- 1
    while (upper < max_upper && surv_i(upper) > tail_tol) {
      upper <- upper * 2
    }
    
    integ <- try(
      stats::integrate(
        surv_i,
        lower = 0,
        upper = upper,
        rel.tol = rel.tol,
        subdivisions = subdivisions
      ),
      silent = TRUE
    )
    
    if (inherits(integ, "try-error")) {
      out[i] <- NA_real_
    } else {
      out[i] <- integ$value
    }
  }
  
  out
}

# ------------------------------------------------------------------------------
# 4. Posterior-draw-specific prediction for expected sequencing delay
#    Documentation -> delay framework:
#      response  = observed sequencing delay
#      prediction = expected delay E[T] under the Gompertz model
# ------------------------------------------------------------------------------

predict_D_delay_single_draw <- function(fit,
                                        newdata,
                                        draw_id,
                                        re_formula = NA) {
  # Linear predictor for mu
  eta_mu <- posterior_linpred(
    fit,
    newdata    = newdata,
    dpar       = "mu",
    draw_ids   = draw_id,
    re_formula = re_formula,
    transform  = FALSE
  )
  eta_mu <- as.numeric(eta_mu)
  
  # Linear predictor for gamma
  eta_gamma <- posterior_linpred(
    fit,
    newdata    = newdata,
    dpar       = "gamma",
    draw_ids   = draw_id,
    re_formula = re_formula,
    transform  = FALSE
  )
  eta_gamma <- as.numeric(eta_gamma)
  
  # Inverse links
  mu    <- exp(eta_mu)   # log link
  gamma <- eta_gamma     # identity link
  
  mean_delay <- gompertz_mean(mu = mu, gamma = gamma)
  
  as.numeric(mean_delay)
}

# Use the model passed by DALEX, rather than capturing fit_D in the closure
make_predict_fun_D_delay <- function(draw_id) {
  force(draw_id)
  function(model, newdata) {
    predict_D_delay_single_draw(
      fit        = model,
      newdata    = newdata,
      draw_id    = draw_id,
      re_formula = NA
    )
  }
}

# ------------------------------------------------------------------------------
# 5. Pre-split data by continent
#    Documentation logic: use only species with realised sequencing delay
# ------------------------------------------------------------------------------

continents <- sort(unique(dat_D$continent))

dat_D_by_continent <- dat_D %>%
  dplyr::filter(event == 1) %>%
  split(.$continent)

job_grid_D <- tidyr::expand_grid(
  continent = continent,
  draw_id   = draw_ids
)

# ------------------------------------------------------------------------------
# 6. Parallel plan
# ------------------------------------------------------------------------------

plan(multisession, workers = 10L)

# ------------------------------------------------------------------------------
# 7. MSE loss function
#    With continuous delay outcome, this is ordinary MSE
# ------------------------------------------------------------------------------

mse_loss <- function(observed, predicted) {
  ok <- is.finite(observed) & is.finite(predicted)
  if (!any(ok)) return(NA_real_)
  mean((observed[ok] - predicted[ok])^2)
}
attr(mse_loss, "loss_name") <- "MSE"

# ------------------------------------------------------------------------------
# 8. Variable-wise importance for sequence acquisition delay
#    Reported as increase in mean squared error (ΔMSE)
# ------------------------------------------------------------------------------

# with_progress({
#   p <- progressor(steps = nrow(job_grid_D))
#   
#   results_D_varwise <- job_grid_D %>%
#     furrr::future_pmap_dfr(function(continent, draw_id) {
#       p(message = paste("continent:", continent, "| Draw:", draw_id))
#       dat_r <- dat_D_by_continent[[continent]]
#       
#       if (is.null(dat_r) || nrow(dat_r) == 0) {
#         return(tibble(
#           continent = character(0),
#           draw_id   = integer(0),
#           variable  = character(0),
#           delta_mse = numeric(0)
#         ))
#       }
#       
#       x_r <- dat_r[, predictors_D, drop = FALSE]
#       y_r <- dat_r$time   # observed sequence acquisition delay
#       
#       pred_fun <- make_predict_fun_D_delay(
#         draw_id = draw_id
#       )
#       
#       expl <- DALEX::explain(
#         model            = fit_D,
#         data             = x_r,
#         y                = y_r,
#         label            = paste0("Darwinian_delay_", continent, "_draw_", draw_id),
#         predict_function = pred_fun,
#         verbose          = TRUE
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
#           continent = continent,
#           draw_id   = draw_id,
#           variable  = variable,
#           delta_mse = dropout_loss
#         )
#     })
# })
# 
# summary_D_varwise <- results_D_varwise %>%
#   dplyr::group_by(continent, draw_id) %>%
#   dplyr::mutate(rel_importance = delta_mse / sum(delta_mse, na.rm = TRUE)) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(continent, variable) %>%
#   dplyr::summarise(
#     rel_imp_mean = mean(rel_importance, na.rm = TRUE),
#     rel_imp_l95  = quantile(rel_importance, 0.025, na.rm = TRUE),
#     rel_imp_u95  = quantile(rel_importance, 0.975, na.rm = TRUE),
#     .groups      = "drop"
#   ) %>%
#   dplyr::arrange(continent, desc(rel_imp_mean))
# 
# saveRDS(results_D_varwise, "output/tables/results_D_varwise_country.rds")

# ------------------------------------------------------------------------------
# 9. Group-wise importance (Biology / Geography / Environment / Access / Activity)
#    Reported as increase in mean squared error (ΔMSE)
# ------------------------------------------------------------------------------

# with_progress({
#   p <- progressor(steps = nrow(job_grid_D))
#   
#   results_D_groupwise <- job_grid_D %>%
#     furrr::future_pmap_dfr(function(continent, draw_id) {
#       p(message = paste("continent:", continent, "| Draw:", draw_id))
#       dat_r <- dat_D_by_continent[[continent]]
#       
#       if (is.null(dat_r) || nrow(dat_r) == 0) {
#         return(tibble(
#           continent = character(0),
#           draw_id   = integer(0),
#           group     = character(0),
#           delta_mse = numeric(0)
#         ))
#       }
#       
#       x_r <- dat_r[, predictors_D, drop = FALSE]
#       y_r <- dat_r$time
#       
#       pred_fun <- make_predict_fun_D_delay(
#         draw_id = draw_id
#       )
#       
#       expl <- DALEX::explain(
#         model            = fit_D,
#         data             = x_r,
#         y                = y_r,
#         label            = paste0("Darwinian_delay_", continent, "_draw_", draw_id),
#         predict_function = pred_fun,
#         verbose          = TRUE
#       )
#       
#       mp_group <- DALEX::model_parts(
#         explainer       = expl,
#         loss_function   = mse_loss,
#         type            = "difference",
#         B               = 10,
#         variable_groups = var_groups_D
#       )
#       
#       mp_group_df <- if ("result" %in% names(mp_group)) mp_group$result else mp_group
#       
#       mp_group_df %>%
#         dplyr::filter(variable %in% names(var_groups_D)) %>%
#         dplyr::transmute(
#           continent = continent,
#           draw_id   = draw_id,
#           group     = variable,
#           delta_mse = dropout_loss
#         )
#     })
# })
# 
# summary_D_groupwise <- results_D_groupwise %>%
#   dplyr::group_by(continent, draw_id) %>%
#   dplyr::mutate(rel_importance = delta_mse / sum(delta_mse, na.rm = TRUE)) %>%
#   dplyr::ungroup() %>%
#   dplyr::group_by(continent, group) %>%
#   dplyr::summarise(
#     rel_imp_mean = mean(rel_importance, na.rm = TRUE),
#     rel_imp_l95  = quantile(rel_importance, 0.025, na.rm = TRUE),
#     rel_imp_u95  = quantile(rel_importance, 0.975, na.rm = TRUE),
#     .groups      = "drop"
#   ) %>%
#   dplyr::arrange(continent, desc(rel_imp_mean))
# 
# saveRDS(results_D_groupwise, "output/tables/results_D_groupwise_country.rds")
# summary_D_varwise    # variable-level importance for sequence acquisition delay
# summary_D_groupwise  # grouped importance
################################################################################
results_D_varwise <- readRDS("output/tables/results_D_varwise_country.rds")

final_D_importance <- results_D_varwise %>%
  filter(!variable %in% c("_baseline_", "_full_model_", "continent")) %>%
  group_by(continent, draw_id) %>%
  mutate(
    delta_pos = pmax(delta_mse, 0),
    denom_pos = sum(delta_pos, na.rm = TRUE),
    rel_imp_pos = if_else(denom_pos > 0, delta_pos / denom_pos, NA_real_)
  ) %>%
  ungroup() %>%
  group_by(continent, variable) %>%
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
  arrange(continent, desc(delta_mse_median))

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
  taxonomic_effort     = "Taxonomic effort",
  taxonomic_activity   = "Taxonomic activity",
  sampling_effort      = "Sampling effort"
  #sequencing_effort    = "Sequencing effort"
)

summary_D_varwise <- final_D_importance %>%
  mutate(
    variable = recode(variable, !!!rename_map),
    variable = factor(variable, levels = rev(unname(rename_map))),
    group = case_when(
      variable %in% c("Sampling effort",
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

summary_D_varwise$continent <- factor(summary_D_varwise$continent,levels = c( "North America",
                                                                              "South America",
                                                                              "Europe",
                                                                              "Africa",
                                                                              "Asia",
                                                                              "Oceania"))

p1 <- ggplot(data = summary_D_varwise, aes(x = delta_mse_median, y = variable, colour = continent)) +
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
    "Continent",
    values = c(
      "Africa"        = "#D55E00",
      "Oceania"       = "#0072B2",
      "Asia"          = "#009E73",
      "North America" = "#F0E442",
      "South America" = "#E69F00",
      "Europe"        = "#CC79A7"
    )
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.05))) +
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
  draw_plot(map,x = 0.2, y = 0.4, scale = 0.3)


results_D_groupwise <- readRDS("output/tables/results_D_groupwise_country.rds")

results_D_groupwise_clean <- results_D_groupwise %>%
  group_by(continent, draw_id, group) %>%
  summarise(
    delta_mse = sum(delta_mse, na.rm = TRUE),
    .groups = "drop"
  )

groupwise_draw_rel <- results_D_groupwise_clean %>%
  group_by(continent, draw_id) %>%
  mutate(
    delta_pos = pmax(delta_mse, 0),
    denom_pos = sum(delta_pos, na.rm = TRUE),
    rel_imp = if_else(denom_pos > 0, delta_pos / denom_pos, NA_real_)
  ) %>%
  ungroup()

summary_D_groupwise_rel <- groupwise_draw_rel %>%
  group_by(continent, group) %>%
  summarise(
    rel_imp_median = median(rel_imp, na.rm = TRUE),
    rel_imp_l95    = quantile(rel_imp, 0.025, na.rm = TRUE),
    rel_imp_u95    = quantile(rel_imp, 0.975, na.rm = TRUE),
    delta_mse_median = median(delta_mse, na.rm = TRUE),
    p_negative = mean(delta_mse < 0, na.rm = TRUE),
    .groups = "drop"
  )

global <- summary_D_groupwise_rel %>%
  group_by(group) %>%
  summarise(
    rel_imp_median = mean(rel_imp_median, na.rm = TRUE),
    rel_imp_l95    = mean(rel_imp_l95, na.rm = TRUE),
    rel_imp_u95    = mean(rel_imp_u95, na.rm = TRUE),
    delta_mse_median = mean(delta_mse_median, na.rm = TRUE),
    p_negative = mean(p_negative, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(continent = "Global")

dt <- bind_rows(summary_D_groupwise_rel, global) %>%
  mutate(
    rel_imp_pct = rel_imp_median * 100
  )

dt <- dt %>%
  group_by(continent) %>%
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


dt$continent <- factor(dt$continent,levels = c( "North America",
                                                "South America",
                                                "Europe",
                                                "Africa",
                                                "Asia",
                                                "Oceania","Global"))
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
  facet_grid(. ~ continent) +
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


ggsave("figures/supplement/Figure_S10.png",dpi = 300, units = "cm", width = 20, height = 20)