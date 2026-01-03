# ------------------------------------------------------------------------------
# Supplementary Figure S17 Sensitivity of shortfall‐combination scenario classification to quantile thresholds
# ------------------------------------------------------------------------------
library(dplyr)
library(cowplot)
library(tidyr)
library(ggplot2)
library(purrr)
library(tidyr)
library(patchwork)
linnaean <- read.csv("output/tables/basin_linnaean_shortfall.csv")
wallacean <- read.csv("output/tables/basin_wallacean_shortfall.csv")
darwinian <- read.csv("output/tables/basin_darwinian_shortfall.csv")

df_basin <- linnaean[,c(1,10)] %>% 
  left_join(wallacean[,c(1,6)], by = "basin_id") %>%
  left_join(darwinian[,c(1,6)], by = "basin_id")


biogeographic <- read.csv("input/raw/biogeographic_list.csv")
df <- df_basin %>% left_join(biogeographic,by = "basin_id")
names(df)[2:4] <- c("LS","WS","DS")
head(df)
df$biogeographic_realm <- as.factor(df$biogeographic_realm)
df <- df %>% na.omit()

classify_within_realm_quantile <- function(data,
                                           value_col,
                                           realm_col = "biogeographic_realm",
                                           q = 0.25) {
  stopifnot(value_col %in% names(data), realm_col %in% names(data))
  stopifnot(q > 0 & q < 0.5)
  
  data %>%
    group_by(.data[[realm_col]]) %>%
    mutate(
      hi = quantile(.data[[value_col]], 1 - q, na.rm = TRUE, type = 7),
      lo = quantile(.data[[value_col]], q,     na.rm = TRUE, type = 7),
      state = case_when(
        .data[[value_col]] >= hi ~ "high",
        .data[[value_col]] <= lo ~ "low",
        TRUE ~ NA_character_
      )
    ) %>%
    ungroup() %>%
    pull(state)
}


qs <- seq(0.01, 0.35, by = 0.01)

scenario_by_q <- map_dfr(qs, function(qi) {
  
  df_tmp <- df %>%
    mutate(
      LS_state = classify_within_realm_quantile(., "LS", q = qi),
      WS_state = classify_within_realm_quantile(., "WS", q = qi),
      DS_state = classify_within_realm_quantile(., "DS", q = qi),
      scenario = case_when(
        LS_state == "high" & WS_state == "high" & DS_state == "high" ~ "LS↑ WS↑ DS↑",
        LS_state == "high" & WS_state == "high" & DS_state == "low"  ~ "LS↑ WS↑ DS↓",
        LS_state == "high" & WS_state == "low"  & DS_state == "high" ~ "LS↑ WS↓ DS↑",
        LS_state == "high" & WS_state == "low"  & DS_state == "low"  ~ "LS↑ WS↓ DS↓",
        LS_state == "low"  & WS_state == "high" & DS_state == "high" ~ "LS↓ WS↑ DS↑",
        LS_state == "low"  & WS_state == "high" & DS_state == "low"  ~ "LS↓ WS↑ DS↓",
        LS_state == "low"  & WS_state == "low"  & DS_state == "high" ~ "LS↓ WS↓ DS↑",
        TRUE ~ NA_character_
      )
    )
  
  df_tmp %>%
    count(scenario) %>%
    mutate(q = qi)
})

basin <- scenario_by_q %>%
  filter(!is.na(scenario)) %>%
  arrange(q, desc(n))

rm(biogeographic,linnaean,wallacean,darwinian,df,df_basin,scenario_by_q,classify_within_realm_quantile)

################################################################################
linnaean <- read.csv("output/tables/country_linnaean_shortfall.csv")
wallacean <- read.csv("output/tables/country_wallacean_shortfall.csv")
darwinian <- read.csv("output/tables/country_darwinian_shortfall.csv")

df_country <- linnaean[,c(1,10)] %>% 
  left_join(wallacean[,c(1,6)], by = "iso3") %>%
  left_join(darwinian[,c(1,6)], by = "iso3")

country_list <- read.csv("input/raw/country_list.csv")
df <- df_country %>% left_join(country_list,by = "iso3")
names(df)[2:4] <- c("LS","WS","DS")
head(df)
df$continent <- as.factor(df$continent)
df <- df %>% na.omit()
classify_within_continent_quantile <- function(data,
                                               value_col,
                                               realm_col = "continent",
                                               q = 0.25) {
  stopifnot(value_col %in% names(data), realm_col %in% names(data))
  stopifnot(q > 0 & q < 0.5)
  
  data %>%
    group_by(.data[[realm_col]]) %>%
    mutate(
      hi = quantile(.data[[value_col]], 1 - q, na.rm = TRUE, type = 7),
      lo = quantile(.data[[value_col]], q,     na.rm = TRUE, type = 7),
      state = case_when(
        .data[[value_col]] >= hi ~ "high",
        .data[[value_col]] <= lo ~ "low",
        TRUE ~ NA_character_
      )
    ) %>%
    ungroup() %>%
    pull(state)
}


scenario_by_q <- map_dfr(qs, function(qi) {
  
  df_tmp <- df %>%
    mutate(
      LS_state = classify_within_continent_quantile(., "LS", q = qi),
      WS_state = classify_within_continent_quantile(., "WS", q = qi),
      DS_state = classify_within_continent_quantile(., "DS", q = qi),
      scenario = case_when(
        LS_state == "high" & WS_state == "high" & DS_state == "high" ~ "LS↑ WS↑ DS↑",
        LS_state == "high" & WS_state == "high" & DS_state == "low"  ~ "LS↑ WS↑ DS↓",
        LS_state == "high" & WS_state == "low"  & DS_state == "high" ~ "LS↑ WS↓ DS↑",
        LS_state == "high" & WS_state == "low"  & DS_state == "low"  ~ "LS↑ WS↓ DS↓",
        LS_state == "low"  & WS_state == "high" & DS_state == "high" ~ "LS↓ WS↑ DS↑",
        LS_state == "low"  & WS_state == "high" & DS_state == "low"  ~ "LS↓ WS↑ DS↓",
        LS_state == "low"  & WS_state == "low"  & DS_state == "high" ~ "LS↓ WS↓ DS↑",
        TRUE ~ NA_character_
      )
    )
  
  df_tmp %>%
    count(scenario) %>%
    mutate(q = qi)
})

country <- scenario_by_q %>%
  filter(!is.na(scenario)) %>%
  arrange(q, desc(n))

rm(country_list,linnaean,wallacean,darwinian,df,df_country,scenario_by_q,classify_within_continent_quantile)


# Quantile thresholds for scenario classification were selected based on analytical scale. For drainage basins, 
# upper and lower quartiles (q = 0.25) were used to ensure broad spatial coverage. For country-level analyses, 
# a slightly more conservative threshold (q = 0.20) was adopted to account for smaller 
# sample sizes and to avoid over-expansion of dominant scenarios. Sensitivity analyses across alternative thresholds yielded consistent patterns.


scenario_cols <- c(
  "LS↑ WS↑ DS↑" = "#D73027",  
  "LS↑ WS↓ DS↓" = "#FC8D59",  
  "LS↑ WS↓ DS↑" = "#E6AB02",  
  "LS↑ WS↑ DS↓" = "#66A61E",  
  "LS↓ WS↑ DS↑" = "#1B9E77",  
  "LS↓ WS↑ DS↓" = "#7570B3",  
  "LS↓ WS↓ DS↑" = "#4DAF4A"   
)

p1 <- ggplot(
  basin %>% filter(!is.na(scenario)),
  aes(x = q, y = n, colour = scenario)
) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 0.25, linetype = "dashed", linewidth = 0.8) +
  scale_colour_manual(values = scenario_cols) +
  labs(
    x = "Quantile threshold (q)",
    y = "Number of basins",
    colour = "Scenario"
  ) +
  theme_minimal(base_size = 11)+
  annotate(
    "text",
    x = 0.25, y = Inf,
    label = "q = 0.25",
    vjust = 1.3, hjust = -0.1,
    size = 3.5
  ) 
p2 <- ggplot(
  country %>% filter(!is.na(scenario)),
  aes(x = q, y = n, colour = scenario)
) +
  geom_line(linewidth = 0.9, show.legend = F) +
  geom_vline(xintercept = 0.3, linetype = "dashed", linewidth = 0.8) +
  scale_colour_manual(values = scenario_cols) +
  labs(
    x = "Quantile threshold (q)",
    y = "Number of countries",
    colour = "Scenario"
  ) +
  theme_minimal(base_size = 11)+
  annotate(
    "text",
    x = 0.3, y = Inf,
    label = "q = 0.30",
    vjust = 1.3, hjust = -0.1,
    size = 3.5
  ) 

final_plot <- p1 + p2 +
  plot_layout(guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold")
    )
  )

final_plot
ggsave("figures/supplement/Figure_S16.png",width = 20, height = 10, units = "cm")

