# ============================================================
# Fig.4 Prioritization of freshwater fish biodiversity collections based on multi-objective trade-offs
# ============================================================

rm(list = ls())

library(dplyr)
library(ggplot2)
library(GGally)
library(rPref)
library(patchwork)
library(legendry)
library(ggh4x)
library(cowplot)
library(here)

# extra libs used later in your code
library(emoa)     # Pareto dominance (minimization) emoa::nds_rank()
library(ggrepel)  # geom_text_repel
library(sf)
library(rnaturalearth)
library(ggrastr)

sf_use_s2(FALSE)

# -----------------------------
# Helpers
# -----------------------------
read_rds_h  <- function(...) readRDS(here::here(...))
read_csv_h  <- function(...) read.csv(here::here(...), stringsAsFactors = FALSE)

min_max_normalize_safe <- function(x, epsilon = 1e-6) {
  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)
  (x - min_val) / (max_val - min_val + epsilon)
}

# Realm-specific Pareto frontier helper
add_realm_front <- function(df, obj_cols) {
  df %>%
    group_by(biogeographic_realm) %>%
    group_modify(~ {
      obj_mat <- as.matrix(.x[, obj_cols, drop = FALSE])
      fronts  <- emoa::nds_rank(t(obj_mat))
      stopifnot(length(fronts) == nrow(.x))
      .x %>%
        mutate(
          front        = fronts,
          front_rank   = 1 - percent_rank(front),              # 1 = best (front=1), 0 = worst
          front_scaled = rank(-front, ties.method = "min")     # larger = better
        )
    }) %>%
    ungroup()
}

# tag best/worst basins by q-quantiles of front_scaled
tag_strategy <- function(df, q = 0.10) {
  df %>%
    group_by(biogeographic_realm) %>%
    mutate(
      strategy = case_when(
        front_scaled >= quantile(front_scaled, 1 - q, na.rm = TRUE) ~ "best",
        front_scaled <= quantile(front_scaled, q,     na.rm = TRUE) ~ "worst",
        TRUE ~ "not_significant"
      )
    ) %>%
    ungroup()
}

# Build a consistent “Pareto scatter + frontier path + labels” plot
# NOTE: keep your aesthetics, just encapsulate to reduce repetition
plot_pareto_2d <- function(df, xvar, yvar,
                           best_pref, worst_pref,
                           xlab, ylab,
                           legend = TRUE,
                           legend_pos = c(0.2, 0.90)) {
  
  df_best  <- df %>% group_by(biogeographic_realm) %>% group_modify(~ psel(.x, best_pref))  %>% ungroup()
  df_worst <- df %>% group_by(biogeographic_realm) %>% group_modify(~ psel(.x, worst_pref)) %>% ungroup()
  
  # label each realm at its left-most point on the chosen frontier
  lab_best <- df_best %>%
    group_by(biogeographic_realm) %>%
    filter(.data[[xvar]] == min(.data[[xvar]], na.rm = TRUE)) %>%
    slice(1) %>%
    ungroup()
  
  lab_worst <- df_worst %>%
    group_by(biogeographic_realm) %>%
    filter(.data[[xvar]] == min(.data[[xvar]], na.rm = TRUE)) %>%
    slice(1) %>%
    ungroup()
  
  p <- ggplot() +
    geom_point(
      data = df,
      aes(x = .data[[xvar]], y = .data[[yvar]], colour = front_rank),
      size = 1, alpha = 0.5
    ) +
    # --- best frontier
    geom_point(
      data = df_best,
      aes(x = .data[[xvar]], y = .data[[yvar]], group = biogeographic_realm),
      shape = 21, fill = "#00B050", colour = "black", size = 1, stroke = 0.3
    ) +
    geom_path(
      data = df_best %>% group_by(biogeographic_realm) %>% arrange(.data[[xvar]], .by_group = TRUE) %>% ungroup(),
      aes(x = .data[[xvar]], y = .data[[yvar]], lty = biogeographic_realm),
      linewidth = 0.2, colour = "#00B050", show.legend = FALSE,
      arrow = arrow(type = "closed", length = unit(0.2, "cm"))
    ) +
    geom_text_repel(
      data = lab_best,
      aes(x = .data[[xvar]], y = .data[[yvar]], label = biogeographic_realm),
      nudge_x = 0.03,direction = "x", point.padding = 0.2,vjust = -2,
      segment.color = NA, size = 2, colour = "black"
    ) +
    # --- worst frontier
    geom_point(
      data = df_worst,
      aes(x = .data[[xvar]], y = .data[[yvar]], group = biogeographic_realm),
      shape = 21, fill = "#FFC000", colour = "black", size = 1, stroke = 0.3
    ) +
    geom_path(
      data = df_worst %>% group_by(biogeographic_realm) %>% arrange(.data[[xvar]], .by_group = TRUE) %>% ungroup(),
      aes(x = .data[[xvar]], y = .data[[yvar]], lty = biogeographic_realm),
      linewidth = 0.2, colour = "#FFC000", show.legend = FALSE,
      arrow = arrow(type = "closed", length = unit(0.2, "cm"))
    ) +
    geom_text_repel(
      data = lab_worst,
      aes(x = .data[[xvar]], y = .data[[yvar]], label = biogeographic_realm),
      nudge_x = 0.03,direction = "y", point.padding = 0.2,hjust = -2,
      segment.color = NA, size = 2, colour = "black"
    ) +
    scale_colour_gradientn(
      name = "Pareto frontier",
      colours = c("#FFC000", "#BFBFBF", "#00B050"),
      limits = c(0, 1), breaks = c(0, 1), labels = c("Worst", "Best"),
      guide = guide_colbar(show = c(TRUE, TRUE), oob = "squish")
    ) +
    labs(x = xlab, y = ylab) +
    theme_classic() +
    theme(
      axis.text  = element_text(colour = "black", size = 7),
      axis.title = element_text(colour = "black", size = 8),
      legend.title = element_text(size = 6, hjust = 0.5),
      legend.text  = element_text(size = 5),
      legend.key.width = unit(0.3, "cm"),
      legend.key.height = unit(0.15, "cm"),
      legend.background = element_blank(),
      legend.title.position = "top",
      legend.direction = "horizontal",
      plot.margin = margin(0, 0, 0, 0)
    )
  
  if (legend) {
    p + theme(legend.position = legend_pos)
  } else {
    p + theme(legend.position = "none")
  }
}

# Denominator for efficiency E = Shortfall / cost(scenario)
get_cost <- function(df, scenario) {
  if (scenario == "scenario_1") {
    df$Investment
  } else if (scenario == "scenario_2") {
    df$Socioeconomics
  } else if (scenario == "scenario_3") {
    df$Governance
  } else if (scenario == "scenario_4") {
    df$Investment * df$Socioeconomics * df$Governance
  } else {
    stop("Unknown scenario: ", scenario)
  }
}

#' Compute Best-Worst Ratio of Efficiency (BWR)
#'
#' This function calculates the ratio of mean efficiency between `best` and `worst` strategy groups,
#' along with key statistical metrics (sample size, standard error, coefficient of variation, 95% CI)
#' required for scientific reporting.
#' @param df A data frame containing strategy labels, scenario identifiers and target variable.
#' @param scenario Character string specifying the scenario to filter.
#' @param S_var Character string of the target variable name for efficiency calculation.
#' @param conf_level Numeric value of confidence level for CI calculation (default = 0.95).
#' @return A tibble with comprehensive statistics of BWR_S.
#' @export
compute_BWR <- function(df, scenario, S_var, conf_level = 0.95) {
  if (!is.data.frame(df)) stop("Argument 'df' must be a valid data frame")
  if (!scenario %in% df$scenario) stop("Scenario '", scenario, "' not found in df$scenario")
  if (!S_var %in% colnames(df)) stop("Variable '", S_var, "' not found in data frame columns")
  
  df_filtered <- df %>%
    filter(scenario == !!scenario) %>%
    tidyr::drop_na(dplyr::all_of(S_var))
  
  best_grp  <- df_filtered %>% filter(strategy == "best")
  worst_grp <- df_filtered %>% filter(strategy == "worst")
  
  E_best <- if (nrow(best_grp)  > 0) best_grp[[S_var]]  / get_cost(best_grp,  scenario) else numeric(0)
  E_worst <- if (nrow(worst_grp) > 0) worst_grp[[S_var]] / get_cost(worst_grp, scenario) else numeric(0)
  
  E_best  <- E_best[is.finite(E_best) & E_best > 0]
  E_worst <- E_worst[is.finite(E_worst) & E_worst > 0]
  
  if (length(E_best) == 0 | length(E_worst) == 0) {
    return(tibble::tibble(
      scenario = scenario, S_var = S_var,
      n_best = length(E_best), n_worst = length(E_worst),
      mean_best = NA_real_, mean_worst = NA_real_, mean_ratio = NA_real_,
      se_best = NA_real_, se_worst = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_,
      cv_best = NA_real_, cv_worst = NA_real_
    ))
  }
  
  mean_best  <- mean(E_best)
  mean_worst <- mean(E_worst)
  se_best    <- sd(E_best)  / sqrt(length(E_best))
  se_worst   <- sd(E_worst) / sqrt(length(E_worst))
  cv_best    <- (sd(E_best)  / mean_best)  * 100
  cv_worst   <- (sd(E_worst) / mean_worst) * 100
  mean_ratio <- mean_best / mean_worst
  
  se_log_ratio <- sqrt((se_best / mean_best)^2 + (se_worst / mean_worst)^2)
  df_total     <- length(E_best) + length(E_worst) - 2
  t_crit       <- qt((1 + conf_level) / 2, df = df_total)
  
  log_ratio <- log(mean_ratio)
  ci_lower  <- exp(log_ratio - t_crit * se_log_ratio)
  ci_upper  <- exp(log_ratio + t_crit * se_log_ratio)
  
  tibble::tibble(
    scenario = scenario,
    S_var = S_var,
    n_best = length(E_best),
    n_worst = length(E_worst),
    mean_best = mean_best,
    mean_worst = mean_worst,
    mean_ratio = mean_ratio,
    se_best = se_best,
    se_worst = se_worst,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    cv_best = cv_best,
    cv_worst = cv_worst
  )
}

# -----------------------------
# Data I/O (here)
# -----------------------------
cost      <- read_rds_h("input", "processed", "cost_basin_all.rds")
sdg       <- read_rds_h("input", "processed", "sdg_basin.rds")
shortfall <- read_rds_h("input", "processed", "basin_shortfall.rds")

biogeographic_list <- read_csv_h("input", "raw", "biogeographic_list.csv") %>%
  mutate(
    biogeographic_realm = recode(
      biogeographic_realm,
      "Afrotropic"   = "AFR",
      "Australasia"  = "AUS",
      "Indomalayan"  = "IND",
      "Nearctic"     = "NEA",
      "Neotropic"    = "NEO",
      "Palearctic"   = "PAL",
      "Oceania"      = "OCE"
    )
  )

# -----------------------------
# PCA: SDG → Socioeconomics / Governance
# -----------------------------
sdg.pca <- prcomp(sdg[, 2:18], center = TRUE, scale. = TRUE)
sdg.pca$rotation

# PC1(Socioeconomics)# ensure higher value = stronger socioeconomic capacity
# Governance = PC2  # already correct direction
sdg$PC1 <- -sdg.pca$x[, 1]
sdg$PC2 <-  sdg.pca$x[, 2]
names(sdg)[19:20] <- c("Socioeconomics", "Governance")

# -----------------------------
# Merge to analysis table
# -----------------------------
data <- shortfall %>%
  select(basin_id = basin_id, Shortfall = NSCI_norm) %>%
  left_join(cost %>% select(basin_id = basin_id, Investment = cost_scale), by = "basin_id") %>%
  left_join(sdg  %>% select(basin_id, Socioeconomics, Governance), by = "basin_id") %>%
  left_join(biogeographic_list[, c(1, 3)], by = "basin_id")

basin_shortfalls <- read_csv_h("output", "tables", "basin_shortfalls.csv") %>%
  select(basin_id, scenario) %>%
  na.omit() %>%
  pull(basin_id)

# data <- data %>% filter(basin_id %in% basin_shortfalls)
data <- data %>% filter(basin_id %in% basin_shortfalls)

# Normalize 0–1 (keep your intent)
data <- data %>%
  mutate(
    Investment     = min_max_normalize_safe(Investment),
    Shortfall      = min_max_normalize_safe(Shortfall),
    Socioeconomics = min_max_normalize_safe(Socioeconomics),
    Governance     = min_max_normalize_safe(Governance)
  )

rm(cost, sdg, sdg.pca, shortfall, biogeographic_list, basin_shortfalls)

################################################################################
# Realm-specific Pareto frontier

#High Shortfall × High Socioeconomics（Best）
#High Shortfall × Low Socioeconomics（Worst）
# Figure X. Trade-off between knowledge shortfall and socioeconomic feasibility, with Pareto-optimal strategic priorities highlighted.

# ------------------------------------------------------------
# Prepare data and specify objective directions
#    Shortfall  : higher is better  → multiply by -1
#    Investment : lower is better  → keep unchanged
# ------------------------------------------------------------
df1 <- data %>%
  select(basin_id, biogeographic_realm, Shortfall, Investment) %>%
  filter(!is.na(Investment), !is.na(Shortfall)) %>%
  mutate(
    # Convert objectives to minimization form for nds_rank():
    # High Shortfall (good) → -Shortfall (so smaller = better)
    Shortfall_obj  = -Shortfall,
    Investment_obj =  Investment # low Investment (good) → keep as original value
  )

# realm-conditional prioritization
# To avoid dominance of well-sampled and highly accessible regions in a single global ranking, 
# we conducted realm-conditioned prioritization and identified high-priority basins within each biogeographic realm, 
# ensuring representation across major biogeographic species pools while retaining comparable multi-objective trade-offs.
df_plot_realm_1 <- add_realm_front(df1, c("Shortfall_obj", "Investment_obj")) %>%
  tag_strategy(q = 0.10)

p1 <- plot_pareto_2d(
  df_plot_realm_1,
  xvar = "Investment", yvar = "Shortfall",
  best_pref  = high(Shortfall) * low(Investment),
  worst_pref = low(Shortfall)  * high(Investment),
  xlab = "Investment", ylab = "Shortfall (NSCI)",
  legend = TRUE, legend_pos = c(0.2, 0.90)
)

################################################################################
# Prepare data: Shortfall + Socioeconomics
df2 <- data %>%
  select(basin_id, biogeographic_realm, Shortfall, Socioeconomics) %>%
  filter(!is.na(Shortfall), !is.na(Socioeconomics)) %>%
  mutate(
    # Both objectives need to be maximized,
    # so convert them into minimization for nds_rank():
    Shortfall_obj = -Shortfall,       # higher shortfall = better priority
    Socio_obj     = -Socioeconomics   # higher socioeconomic capacity = better
  )

df_plot_realm_2 <- add_realm_front(df2, c("Shortfall_obj", "Socio_obj")) %>%
  tag_strategy(q = 0.10)

p2 <- plot_pareto_2d(
  df_plot_realm_2,
  xvar = "Socioeconomics", yvar = "Shortfall",
  best_pref  = high(Shortfall) * high(Socioeconomics),
  worst_pref = low(Shortfall)  * low(Socioeconomics),
  xlab = "Socioeconomics", ylab = "Shortfall (NSCI)",
  legend = FALSE
)

################################################################################
#    Goal: prioritize basins with
#      - high Shortfall   (large knowledge gaps)
#      - high Governance  (strong capacity to fill those gaps)
# ------------------------------------------------------------
df3 <- data %>%
  select(basin_id, biogeographic_realm, Governance, Shortfall) %>%
  filter(!is.na(Governance), !is.na(Shortfall)) %>%
  mutate(
    # Convert both objectives to "minimization" form for nds_rank():
    # We want to MAXIMIZE Shortfall → use -Shortfall
    Shortfall_obj  = -Shortfall,
    # We want to MAXIMIZE Governance → use -Governance
    Governance_obj = -Governance
  )

df_plot_realm_3 <- add_realm_front(df3, c("Shortfall_obj", "Governance_obj")) %>%
  tag_strategy(q = 0.10)

p3 <- plot_pareto_2d(
  df_plot_realm_3,
  xvar = "Governance", yvar = "Shortfall",
  best_pref  = high(Shortfall) * high(Governance),
  worst_pref = low(Shortfall)  * low(Governance),
  xlab = "Governance", ylab = "Shortfall (NSCI)",
  legend = FALSE
)

################################################################################
# ------------------------------------------------------------
#  Prepare data for 4-objective trade-off:
#    Goal: High Shortfall × Low Investment × High Socioeconomics × High Governance
#
#    Assumptions:
#    - Shortfall:       higher = larger knowledge gaps (priority ↑)
#    - Investment:      higher = higher total cost (priority ↓)
#    - Socioeconomics:  higher = stronger socioeconomic capacity
#    - Governance:      higher = stronger institutional/governance capacity
# ------------------------------------------------------------
df4 <- data %>%
  select(basin_id, biogeographic_realm, Investment, Socioeconomics, Governance, Shortfall) %>%
  na.omit() %>%
  mutate(
    # Convert all objectives to "minimize" form for nds_rank():
    # We want to MAXIMIZE Shortfall → minimize -Shortfall
    Shortfall_obj      = -Shortfall,
    # We want to MINIMIZE Investment → keep as is
    Investment_obj     =  Investment,
    # We want to MAXIMIZE Socioeconomics → minimize -Socioeconomics
    Socioeconomics_obj = -Socioeconomics,
    # We want to MAXIMIZE Governance → minimize -Governance
    Governance_obj     = -Governance
  )

df_plot_realm_4 <- add_realm_front(df4, c("Shortfall_obj", "Investment_obj", "Socioeconomics_obj", "Governance_obj")) %>%
  tag_strategy(q = 0.10)

p4 <- GGally::ggparcoord(
  df_plot_realm_4,
  columns = c(7, 8, 9, 10),          # keep your original column indices
  groupColumn = 12,                  # front_rank
  scale = "uniminmax",
  showPoints = FALSE,
  alphaLines = 0.15
) +
  scale_colour_gradientn(colours = c("#FFC000", "grey90", "#00B050")) +
  theme_minimal() +
  scale_x_discrete(
    expand = c(0.02, 0.05),
    labels = c("Shortfall \n (NSCI)", "Investment", "Socioeconomics", "Governance")
  ) +
  xlab("") +
  ylab("") +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey50", linewidth = 1),
    legend.position = "none",
    axis.text.x = element_text(colour = "black", size = 6.5),
    axis.text.y = element_text(colour = "black", size = 7),
    plot.margin = margin(0, 0, 0, 0)
  )

################################################################################
# Best and worst basins were defined as the upper and lower q-quantiles of the continuous 
# Pareto prioritization gradient (front_scaled), with sensitivity analyses across alternative thresholds.

front <- bind_rows(
  df_plot_realm_1 %>% select(biogeographic_realm, basin_id, front_scaled, strategy) %>% mutate(scenario = "scenario_1"),
  df_plot_realm_2 %>% select(biogeographic_realm, basin_id, front_scaled, strategy) %>% mutate(scenario = "scenario_2"),
  df_plot_realm_3 %>% select(biogeographic_realm, basin_id, front_scaled, strategy) %>% mutate(scenario = "scenario_3"),
  df_plot_realm_4 %>% select(biogeographic_realm, basin_id, front_scaled, strategy) %>% mutate(scenario = "scenario_4")
)

shortfall2 <- read_rds_h("input", "processed", "basin_shortfall.rds")
dt <- front %>% left_join(shortfall2[, c(1, 7:9)], by = "basin_id")
names(dt)[6:8] <- c("LS", "WS", "DS")

dt <- dt %>%
  mutate(
    LS = min_max_normalize_safe(LS, epsilon = 1e-3),
    WS = min_max_normalize_safe(WS, epsilon = 1e-3),
    DS = min_max_normalize_safe(DS, epsilon = 1e-3)
  )

dt <- dt %>% left_join(data[, 1:5], by = "basin_id")

biogeographic <- read_csv_h("input", "raw", "biogeographic_list.csv")
dt <- dt %>% left_join(biogeographic, by = "basin_id")

# -----------------------------
# BWR results (per scenario × LS/WS/DS)
# -----------------------------
scenarios  <- unique(dt$scenario)
shortfalls <- c("LS", "WS", "DS")

BWR_results <- tidyr::expand_grid(
  scenario = scenarios,
  S_var = shortfalls
) %>%
  mutate(stats = purrr::map2(scenario, S_var, ~ compute_BWR(dt, .x, .y))) %>%
  select(-scenario, -S_var) %>%
  tidyr::unnest(cols = stats) %>%
  rename(
    bwr_mean_ratio = mean_ratio,
    ci_95_lower    = ci_lower,
    ci_95_upper    = ci_upper,
    cv_best_pct    = cv_best,
    cv_worst_pct   = cv_worst
  )

BWR_results$S_var <- factor(BWR_results$S_var, levels = c("LS", "WS", "DS"))
BWR_results

# small barplot helper to avoid 4 copies
plot_bwr_bar <- function(df, scenario_id) {
  ggplot(
    df %>% filter(scenario == scenario_id),
    aes(x = S_var, y = bwr_mean_ratio, ymin = ci_95_lower, ymax = ci_95_upper)
  ) +
    geom_col(fill = "#3182bd", alpha = 0.5, width = 0.6) +
    geom_errorbar(width = 0.2, linewidth = 0.2, color = "grey") +
    labs(x = "Shortfall", y = "BWR (best-to-worst ratio)") +
    theme_classic() +
    scale_y_continuous(limits = c(0, 5), expand = c(0, 0)) +
    theme(
      axis.text = element_text(size = 5, color = "black"),
      axis.title = element_text(size = 6, color = "black"),
      axis.line = element_line(linewidth = 0.2),
      axis.ticks = element_line(linewidth = 0.2),
      axis.ticks.length = unit(0.05, "cm"),
      plot.background = element_blank(),
      panel.background = element_blank(),
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0)
    )
}

pp1_r <- plot_bwr_bar(BWR_results, "scenario_1")
pp2_r <- plot_bwr_bar(BWR_results, "scenario_2")
pp3_r <- plot_bwr_bar(BWR_results, "scenario_3")
pp4_r <- plot_bwr_bar(BWR_results, "scenario_4")

################################################################################
# Maps (use here for inputs, keep your styling)
library(sf)
library(rnaturalearth)
library(ggplot2)
sf_use_s2(FALSE)

inland <- read_rds_h("input", "raw", "basin", "basin_sf_v1.rds")
world_map <- rnaturalearth::ne_countries(scale = 50, type = "countries", returnclass = "sf")

moll_proj <- st_crs("+proj=moll")
lat_points <- data.frame(lon = 0, lat = c(-60, 90)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
projected_points <- st_transform(lat_points, crs = moll_proj)
y_limits <- st_coordinates(projected_points)[, "Y"]

# small helper for map panel (avoid repeating the same ggplot block 4 times)
make_front_map <- function(inland, world_map, df_front, fill_col = "front_rank") {
  case <- inland %>% left_join(df_front, by = "basin_id")
  ggplot() +
    ggrastr::rasterise(
      geom_sf(
        data = world_map %>% st_wrap_dateline(options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180")),
        fill = "#d3d3d3", colour = NA
      ),
      dpi = 300
    ) +
    ggrastr::rasterise(
      geom_sf(data = case, aes(fill = .data[[fill_col]]), colour = "white", linewidth = 0.03),
      dpi = 300
    ) +
    scale_fill_gradientn(
      name = "Pareto frontier",
      colours = c("#FFC000", "grey90", "#00B050"),
      limits = c(0, 1),
      breaks = c(0, 1),
      na.value = "grey70",
      labels = c("Worst", "Best"),
      guide = guide_colbar(show = c(TRUE, TRUE), oob = "squish")
    ) +
    theme_void() +
    coord_sf(crs = "+proj=moll +lon_0=0", ylim = y_limits, expand = FALSE) +
    theme(legend.position = "none", plot.margin = margin(0, 0, 0, 0))
}

# scenario maps (front_rank lives in each df_plot_realm_X)
map1 <- make_front_map(inland, world_map, df_plot_realm_1[, c("basin_id", "front_rank")])
map2 <- make_front_map(inland, world_map, df_plot_realm_2[, c("basin_id", "front_rank")])
map3 <- make_front_map(inland, world_map, df_plot_realm_3[, c("basin_id", "front_rank")])
map4 <- make_front_map(inland, world_map, df_plot_realm_4[, c("basin_id", "front_rank")])

# Combine: add your text labels exactly as before
pp1 <- ggdraw() +
  draw_plot(map1, 0, 0, 1, 1) +
  draw_plot(pp1_r, 0.06, 0.02, 0.2, 0.8) +
  draw_label(label = expression(Delta[NSCI] == -0.25), x = 0.58, y = 0.08, hjust = 1, vjust = 1, size = 7.2) +
  draw_label(label = expression("E-test, " * italic(p) * " < 0.01"), x = 0.8, y = 0.08, hjust = 1, vjust = 1, size = 7.2) +
  draw_label(label = expression(BWR[NSCI] * " = 5.78"), x = 0.3, y = 0.98, hjust = 1, vjust = 1, size = 7.5)

pp2 <- ggdraw() +
  draw_plot(map2, 0, 0, 1, 1) +
  draw_plot(pp2_r, 0.06, 0.02, 0.2, 0.8) +
  draw_label(label = expression(Delta[NSCI] == -0.35), x = 0.58, y = 0.08, hjust = 1, vjust = 1, size = 7.2) +
  draw_label(label = expression("E-test, " * italic(p) * " < 0.01"), x = 0.8, y = 0.08, hjust = 1, vjust = 1, size = 7.2) +
  draw_label(label = expression(BWR[NSCI] * " = 2.74"), x = 0.3, y = 0.98, hjust = 1, vjust = 1, size = 7.5)

pp3 <- ggdraw() +
  draw_plot(map3, 0, 0, 1, 1) +
  draw_plot(pp3_r, 0.06, 0.02, 0.2, 0.8) +
  draw_label(label = expression(Delta[NSCI] == -0.34), x = 0.58, y = 0.08, hjust = 1, vjust = 1, size = 7.2) +
  draw_label(label = expression("E-test, " * italic(p) * " < 0.01"), x = 0.8, y = 0.08, hjust = 1, vjust = 1, size = 7.2) +
  draw_label(label = expression(BWR[NSCI] * " = 4.80"), x = 0.3, y = 0.98, hjust = 1, vjust = 1, size = 7.5)

pp4 <- ggdraw() +
  draw_plot(map4, 0, 0, 1, 1) +
  draw_plot(pp4_r, 0.06, 0.02, 0.2, 0.8) +
  draw_label(label = expression(Delta[NSCI] == -0.20), x = 0.58, y = 0.08, hjust = 1, vjust = 1, size = 7.2) +
  draw_label(label = expression("E-test, " * italic(p) * " < 0.01"), x = 0.8, y = 0.08, hjust = 1, vjust = 1, size = 7.2) +
  draw_label(label = expression(BWR[NSCI] * " = 1.66"), x = 0.3, y = 0.98, hjust = 1, vjust = 1, size = 7.5)

################################################################################
# Arrange multi-panel figure
row_pair <- function(left_plot, right_plot, labels) {
  cowplot::plot_grid(
    left_plot  + theme(plot.margin = margin(5, 5, 0, 0, unit = "pt")),
    right_plot + theme(plot.margin = margin(5, 0, 0, 0, unit = "pt")),
    ncol = 2,
    rel_widths = c(0.4, 0.6),
    labels = labels,
    label_size = 9,
    label_fontface = "bold",
    align = "h",
    axis  = "tb"
  )
}

row1 <- row_pair(p1, pp1, labels = c("A", "B"))
row2 <- row_pair(p2, pp2, labels = c("C", "D"))
row3 <- row_pair(p3, pp3, labels = c("E", "F"))
row4 <- row_pair(p4, pp4, labels = c("G", "H"))

fig <- row1 / row2 / row3 / row4
fig

# Save with here()
ggsave(
  filename = here::here("figures", "main", "Figure_4.png"),
  plot = fig,
  dpi = 300,
  units = "cm",
  width = 18,
  height = 20
)
