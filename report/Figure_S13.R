# Supplementary Figure S13 
# Spatial validation of collection prioritization under species withholding


library(ggplot2)
library(dplyr)
basin_scenario <- readRDS("output/species_withholding/basin_scenario.rds")
basin_scenario$shortfall <- factor(basin_scenario$shortfall,levels = c("L","W","D"))
basin_scenario <- basin_scenario %>%
  mutate(
    whithheld = case_when(
      whithheld == "p005" ~ "5%",
      whithheld == "p010" ~ "10%",
      whithheld == "p020" ~ "20%",
      TRUE ~ whithheld
    ),
    whithheld = factor(whithheld, levels = c("5%", "10%", "20%"))
  )

sig_data_basin <- basin_scenario %>%
  group_by(scenario, whithheld, shortfall) %>%
  summarise(
    p_avg = mean(p_value, na.rm = TRUE),  
    y_pos = max(enrichment, na.rm = TRUE) * 1.05,
    .groups = "drop"
  ) %>%
  mutate(
    sig = case_when(
      p_avg < 0.001 ~ "***",
      p_avg < 0.01  ~ "**",
      p_avg < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

p1 <- ggplot(basin_scenario, aes(x = shortfall, y = enrichment, fill = shortfall)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  geom_boxplot(width = 0.5, linewidth = 0.2, outlier.shape = NA,colour = "grey40") +
  geom_jitter(shape = 21, size = 1, width = 0.15, alpha = 0.7, stroke = 0.2) +
  geom_text(
    data = sig_data_basin,
    aes(y = y_pos, label = sig,colour = shortfall),
    size = 4, fontface = "bold"
  ) +
  scale_colour_manual(values = c(
    "L" = "#D55E00",
    "W" = "#435792",
    "D" = "#009E73"
  )) +
  scale_fill_manual(values = c(
    "L" = "#D55E00",
    "W" = "#435792",
    "D" = "#009E73"
  )) +
  facet_grid(rows = vars(whithheld), cols = vars(scenario),scales = "free_y") +
  labs(x = "Shortfall", y = "Enrichment") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "gray95"),
    legend.position = "none"
  )


country_scenario <- readRDS("output/species_withholding/country_scenario.rds")
country_scenario$shortfall <- factor(country_scenario$shortfall,levels = c("L","W","D"))
country_scenario <- country_scenario %>%
  mutate(
    whithheld = case_when(
      whithheld == "p005" ~ "5%",
      whithheld == "p010" ~ "10%",
      whithheld == "p020" ~ "20%",
      TRUE ~ whithheld
    ),
    whithheld = factor(whithheld, levels = c("5%", "10%", "20%"))
  )

sig_data_country <- country_scenario %>%
  group_by(scenario, whithheld, shortfall) %>%
  summarise(
    p_avg = min(p_value, na.rm = TRUE),  
    y_pos = max(enrichment, na.rm = TRUE) * 1.05,
    .groups = "drop"
  ) %>%
  mutate(
    sig = case_when(
      p_avg < 0.001 ~ "***",
      p_avg < 0.01  ~ "**",
      p_avg < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

p2 <- ggplot(country_scenario, aes(x = shortfall, y = enrichment, fill = shortfall)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  geom_boxplot(width = 0.5, linewidth = 0.2, outlier.shape = NA,colour = "grey40") +
  geom_jitter(shape = 21, size = 1, width = 0.15, alpha = 0.7, stroke = 0.2) +
  geom_text(
    data = sig_data_country,
    aes(y = y_pos, label = sig,colour = shortfall),
    size = 4, fontface = "bold"
  ) +
  scale_colour_manual(values = c(
    "L" = "#D55E00",
    "W" = "#435792",
    "D" = "#009E73"
  )) +
  scale_fill_manual(values = c(
    "L" = "#D55E00",
    "W" = "#435792",
    "D" = "#009E73"
  )) +
  facet_grid(rows = vars(whithheld), cols = vars(scenario),scales = "free_y") +
  labs(x = "Shortfall", y = "Enrichment") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "gray95"),
    legend.position = "none"
  )


library(patchwork)  
p1 / p2 +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold"),
        axis.text = element_text(colour = "black",size = 9),
        axis.title = element_text(colour = "black",size = 10)
  )

ggsave("figures/supplement/Figure_S13.png",
       units = "cm",dpi = 300,
       width = 18, height = 20)
