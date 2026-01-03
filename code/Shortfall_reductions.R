#Reductions in shortfall were quantified as differences in mean shortfall between 
#baseline and best-performing strategies (B − A), with uncertainty estimated using non-parametric bootstrap confidence intervals.

library(dplyr)
library(purrr)
library(tidyr)
boot_diff_ci <- function(df,
                         value_col = "Shortfall",
                         group_col = "strategy",
                         A = "best",
                         B = "worst",
                         n_boot = 2000,
                         conf = 0.95,
                         seed = 123) {
  
  set.seed(seed)
  
  alpha <- (1 - conf) / 2
  
  
  obs <- df %>%
    filter(.data[[group_col]] %in% c(A, B)) %>%
    group_by(.data[[group_col]]) %>%
    summarise(mean_val = mean(.data[[value_col]], na.rm = TRUE),
              .groups = "drop") %>%
    pivot_wider(names_from = group_col, values_from = mean_val)
  
  delta_obs <- obs[[B]] - obs[[A]]
  
  # bootstrap
  boot_deltas <- replicate(n_boot, {
    df_boot <- df %>%
      filter(.data[[group_col]] %in% c(A, B)) %>%
      group_by(.data[[group_col]]) %>%
      slice_sample(prop = 1, replace = TRUE)
    
    m <- df_boot %>%
      group_by(.data[[group_col]]) %>%
      summarise(mean_val = mean(.data[[value_col]], na.rm = TRUE),
                .groups = "drop") %>%
      pivot_wider(names_from = group_col, values_from = mean_val)
    
    m[[B]] - m[[A]]
  })
  
  tibble(
    delta_mean = delta_obs,
    delta_lwr  = quantile(boot_deltas, probs = alpha, na.rm = TRUE),
    delta_upr  = quantile(boot_deltas, probs = 1 - alpha, na.rm = TRUE)
  )
}

rr <- read.csv("output/tables/Pareto_frontier_summary_basin_fix.csv")
results_boot <- rr %>%
  group_by(scenario) %>%
  group_modify(~ boot_diff_ci(
    df = .x,
    value_col = "Shortfall",
    group_col = "strategy",
    A = "best",
    B = "worst",
    n_boot = 2000
  )) %>%
  ungroup()

# A tibble: 4 × 4
# scenario   delta_mean delta_lwr delta_upr
# <chr>           <dbl>     <dbl>     <dbl>
# 1 scenario_1     -0.249    -0.297    -0.202
# 2 scenario_2     -0.349    -0.392    -0.304
# 3 scenario_3     -0.343    -0.390    -0.292
# 4 scenario_4     -0.195    -0.239    -0.150

rr %>%
  group_by(scenario) %>%
  group_modify(~ boot_diff_ci(.x, value_col = "LS"))

# A tibble: 4 × 4
# Groups:   scenario [4]
# scenario   delta_mean delta_lwr delta_upr
# <chr>           <dbl>     <dbl>     <dbl>
# 1 scenario_1     -0.145    -0.177   -0.114 
# 2 scenario_2     -0.203    -0.236   -0.171 
# 3 scenario_3     -0.211    -0.245   -0.174 
# 4 scenario_4     -0.108    -0.138   -0.0787


rr %>%
  group_by(scenario) %>%
  group_modify(~ boot_diff_ci(.x, value_col = "WS"))
# scenario   delta_mean delta_lwr delta_upr
# <chr>           <dbl>     <dbl>     <dbl>
# 1 scenario_1     -0.225    -0.267    -0.185
# 2 scenario_2     -0.271    -0.312    -0.230
# 3 scenario_3     -0.269    -0.310    -0.228
# 4 scenario_4     -0.201    -0.239    -0.163


rr %>%
  group_by(scenario) %>%
  group_modify(~ boot_diff_ci(.x, value_col = "DS"))

# A tibble: 4 × 4
# # Groups:   scenario [4]
#   scenario   delta_mean delta_lwr delta_upr
#   <chr>           <dbl>     <dbl>     <dbl>
# 1 scenario_1     -0.191    -0.238    -0.146
# 2 scenario_2     -0.265    -0.306    -0.223
# 3 scenario_3     -0.264    -0.308    -0.216
# 4 scenario_4     -0.169    -0.212    -0.125


dd <- read.csv("output/tables/Pareto_frontier_summary_country_fix.csv")
dd %>%
  group_by(scenario) %>%
  group_modify(~ boot_diff_ci(
    df = .x,
    value_col = "Shortfall",
    group_col = "strategy",
    A = "best",
    B = "worst",
    n_boot = 2000
  )) %>%
  ungroup()

# # A tibble: 4 × 4
# scenario   delta_mean delta_lwr delta_upr
# <chr>           <dbl>     <dbl>     <dbl>
# 1 scenario_1     0.0188    -0.189    0.203 
# 2 scenario_2     0.0940    -0.103    0.277 
# 3 scenario_3     0.112     -0.130    0.303 
# 4 scenario_4    -0.220     -0.392   -0.0429


dd %>%
  group_by(scenario) %>%
  group_modify(~ boot_diff_ci(.x, value_col = "LS"))

# scenario   delta_mean delta_lwr delta_upr
# <chr>           <dbl>     <dbl>     <dbl>
# 1 scenario_1      0.101   -0.0768    0.264 
# 2 scenario_2      0.162   -0.0186    0.337 
# 3 scenario_3      0.160   -0.0440    0.345 
# 4 scenario_4     -0.192   -0.345    -0.0353

dd %>%
  group_by(scenario) %>%
  group_modify(~ boot_diff_ci(.x, value_col = "WS"))
# scenario   delta_mean delta_lwr delta_upr
# <chr>           <dbl>     <dbl>     <dbl>
# 1 scenario_1    -0.140     -0.296  0.000428
# 2 scenario_2    -0.0398    -0.181  0.109   
# 3 scenario_3    -0.0557    -0.214  0.113   
# 4 scenario_4    -0.119     -0.244  0.00442 



dd %>%
  group_by(scenario) %>%
  group_modify(~ boot_diff_ci(.x, value_col = "DS"))

# scenario   delta_mean delta_lwr delta_upr
# <chr>           <dbl>     <dbl>     <dbl>
# 1 scenario_1    -0.125     -0.268    0.0191
# 2 scenario_2    -0.0143    -0.159    0.128 
# 3 scenario_3    -0.104     -0.253    0.0481
# 4 scenario_4    -0.0739    -0.200    0.0596








