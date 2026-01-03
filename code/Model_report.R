library(dplyr)        # Data manipulation verbs
library(magrittr)     # Pipe operator support
library(rstanarm)     # Bayesian regression modeling
library(brms)         # Bayesian multilevel models interface
library(broom.mixed)
library(sjlabelled)
library(sjPlot)

fit_L <- readRDS("output/model/basin_linnaean_all.rds")$lognormal 
fit_W <- readRDS("output/model/basin_wallacean_all.rds")$lognormal 
fit_D <- readRDS("output/model/basin_darwinian_all.rds")$gamma 

tab_model(
  fit_L, fit_W, fit_D,
  auto.label      = FALSE,
  show.ci         = TRUE,
  show.se         = FALSE,
  show.p          = FALSE,
  show.re.var     = TRUE,
  transform       = NULL,
  dv.labels       = c("Linnaean", "Wallacean", "Darwinian"),
  CSS = list(
    table = "border-top: 2px solid black; border-bottom: 2px solid black;",
    tfoot = "border-top: 1px solid black;",
    td = "padding: 4px;",
    th = "padding: 4px; border-bottom: 1px solid black;"
  ),
  file = "output/table_basin_shortfall_models.html"
)


################################################################################
fit_L <- readRDS("output/model/country_linnaean_all.rds")$weibull
fit_W <- readRDS("output/model/country_wallacean_all.rds")$lognormal 
fit_D <- readRDS("output/model/country_darwinian_all.rds")$weibull 

tab_model(
  fit_L, fit_W, fit_D,
  auto.label      = FALSE,
  show.ci         = TRUE,
  show.se         = FALSE,
  show.p          = FALSE,
  show.re.var     = TRUE,
  transform       = NULL,
  dv.labels       = c("Linnaean", "Wallacean", "Darwinian"),
  CSS = list(
    table = "border-top: 2px solid black; border-bottom: 2px solid black;",
    tfoot = "border-top: 1px solid black;",
    td = "padding: 4px;",
    th = "padding: 4px; border-bottom: 1px solid black;"
  ),
  file = "output/table_country_shortfall_models.html"
)


