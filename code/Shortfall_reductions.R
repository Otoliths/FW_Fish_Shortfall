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
  
  # 原始差值
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
results_boot
# A tibble: 5 × 4
# scenario   delta_mean delta_lwr delta_upr
# <chr>           <dbl>     <dbl>     <dbl>
# 1 scenario_2     -0.235    -0.278    -0.185
# 2 scenario_3     -0.338    -0.383    -0.296
# 3 scenario_4     -0.433    -0.472    -0.392
# 4 scenario_5     -0.278    -0.321    -0.233
# 5 scenario_6     -0.247    -0.293    -0.200
rr %>%
  group_by(scenario) %>%
  group_modify(~ boot_diff_ci(.x, value_col = "LS"))

# A tibble: 5 × 4
# Groups:   scenario [5]
# scenario   delta_mean delta_lwr delta_upr
# <chr>           <dbl>     <dbl>     <dbl>
# 1 scenario_2    -0.0871   -0.114    -0.0610
# 2 scenario_3    -0.156    -0.188    -0.129 
# 3 scenario_4    -0.168    -0.197    -0.138 
# 4 scenario_5    -0.0804   -0.103    -0.0572
# 5 scenario_6    -0.0629   -0.0885   -0.0379

rr %>%
  group_by(scenario) %>%
  group_modify(~ boot_diff_ci(.x, value_col = "WS"))
# A tibble: 5 × 4
# Groups:   scenario [5]
# scenario   delta_mean delta_lwr delta_upr
# <chr>           <dbl>     <dbl>     <dbl>
# 1 scenario_2     -0.119    -0.140   -0.0949
# 2 scenario_3     -0.141    -0.160   -0.122 
# 3 scenario_4     -0.193    -0.212   -0.175 
# 4 scenario_5     -0.144    -0.165   -0.122 
# 5 scenario_6     -0.147    -0.170   -0.123 


rr %>%
  group_by(scenario) %>%
  group_modify(~ boot_diff_ci(.x, value_col = "DS"))

# A tibble: 5 × 4
# Groups:   scenario [5]
# scenario   delta_mean delta_lwr delta_upr
# <chr>           <dbl>     <dbl>     <dbl>
# 1 scenario_2     -0.135    -0.164    -0.105
# 2 scenario_3     -0.209    -0.242    -0.178
# 3 scenario_4     -0.264    -0.292    -0.234
# 4 scenario_5     -0.171    -0.198    -0.143
# 5 scenario_6     -0.148    -0.177    -0.117


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

# A tibble: 5 × 4
# scenario   delta_mean delta_lwr delta_upr
# <chr>           <dbl>     <dbl>     <dbl>
# 1 scenario_2     0.184     0.0116    0.367 
# 2 scenario_3     0.101    -0.0318    0.230 
# 3 scenario_4     0.0142   -0.148     0.156 
# 4 scenario_5    -0.158    -0.378     0.0397
# 5 scenario_6    -0.0923   -0.289     0.0787


dd %>%
  group_by(scenario) %>%
  group_modify(~ boot_diff_ci(.x, value_col = "LS"))

# A tibble: 5 × 4
# scenario   delta_mean delta_lwr delta_upr
# <chr>           <dbl>     <dbl>     <dbl>
# 1 scenario_2     0.122    -0.0611     0.298
# 2 scenario_3     0.0581   -0.122      0.236
# 3 scenario_4     0.0412   -0.122      0.212
# 4 scenario_5    -0.0521   -0.193      0.100
# 5 scenario_6     0.0446   -0.118      0.213

dd %>%
  group_by(scenario) %>%
  group_modify(~ boot_diff_ci(.x, value_col = "WS"))
# A tibble: 5 × 4
# Groups:   scenario [5]
# scenario   delta_mean delta_lwr delta_upr
# <chr>           <dbl>     <dbl>     <dbl>
# 1 scenario_2     0.0974    -0.117    0.297 
# 2 scenario_3     0.0641    -0.131    0.252 
# 3 scenario_4     0.0209    -0.178    0.211 
# 4 scenario_5    -0.128     -0.324    0.0778
# 5 scenario_6    -0.0546    -0.270    0.159



dd %>%
  group_by(scenario) %>%
  group_modify(~ boot_diff_ci(.x, value_col = "DS"))

# A tibble: 5 × 4
# Groups:   scenario [5]
# scenario   delta_mean delta_lwr delta_upr
# <chr>           <dbl>     <dbl>     <dbl>
# 1 scenario_2     0.0720   -0.108     0.240 
# 2 scenario_3     0.0827   -0.0649    0.226 
# 3 scenario_4    -0.0182   -0.171     0.132 
# 4 scenario_5    -0.104    -0.301     0.0818
# 5 scenario_6    -0.0776   -0.258     0.0891
