# Supplementary Figure S13 
# Transitions in collection prioritization categories between the integrated baseline Scenario E and equity-adjusted Scenario F across spatial units

rm(list = ls())
library(ggalluvial)
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)

df_strategy_basin <- read.csv("output/tables/basin_strategy_transition.csv")
df_strategy_basin <- df_strategy_basin %>%
  mutate(
    strategy_E = recode(strategy_E,
                        "worst" = "Worst",
                        "not_significant" = "Not significant",
                        "best" = "Best"),
    strategy_F = recode(strategy_F,
                        "worst" = "Worst",
                        "not_significant" = "Not significant",
                        "best" = "Best"),
    strategy_E = factor(strategy_E, levels = c("Best","Not significant","Worst")),
    strategy_F = factor(strategy_F, levels = c("Best","Not significant","Worst"))
  )

df_sankey_basin <- df_strategy_basin %>%
  count(strategy_E, strategy_F, name = "n")


pal <- c(
  "Worst" = "#FFC000",
  "Not significant" = "grey90",
  "Best" = "#00B050"
)

p1 <- ggplot(
  df_sankey_basin,
  aes(axis1 = strategy_E, axis2 = strategy_F, y = n)
) +
  geom_alluvium(
    aes(fill = strategy_E),
    width = 0.22,
    alpha = 0.6,
    color = "grey70",
    linewidth = 0.3
  ) +
  geom_stratum(
    aes(fill = after_stat(stratum)),
    width = 0.22,
    color = "grey35",
    linewidth = 0.5
  ) +
  geom_text(
    stat = "stratum",
    aes(label = paste0(after_stat(stratum), "\n(n = ", after_stat(count), ")")),
    size = 3,
    fontface = "bold",
    lineheight = 0.95
  ) +
  scale_x_discrete(
    limits = c("Scenario E", "Scenario F"),
    expand = c(0.08, 0.08)
  ) +
  scale_fill_manual(
    values = pal,
    breaks = c("Worst", "Not significant", "Best")
  ) +
  labs(
    #title = "Strategy transitions from Scenario E to Scenario F",
    subtitle = paste0(
      "Stable category = ", percent(mean(df_strategy_basin$strategy_E == df_strategy_basin$strategy_F), accuracy = 0.1),
      "   |   Entered Best = ", sum(df_strategy_basin$strategy_E != "Best" & df_strategy_basin$strategy_F == "Best"),
      "   |   Exited Best = ", sum(df_strategy_basin$strategy_E == "Best" & df_strategy_basin$strategy_F != "Best")
    ),
    x = NULL,
    y = "Number of basins",
    fill = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "",
    panel.grid = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 11,hjust = 0.5,vjust = -3),
    axis.text.x = element_text(face = "bold"),
    legend.text = element_text(size = 11)
  )

p1


df_strategy_country <- read.csv("output/tables/country_strategy_transition.csv")

df_strategy_country <- df_strategy_country %>%
  mutate(
    strategy_E = recode(strategy_E,
                        "worst" = "Worst",
                        "not_significant" = "Not significant",
                        "best" = "Best"),
    strategy_F = recode(strategy_F,
                        "worst" = "Worst",
                        "not_significant" = "Not significant",
                        "best" = "Best"),
    strategy_E = factor(strategy_E, levels = c("Best","Not significant","Worst")),
    strategy_F = factor(strategy_F, levels = c("Best","Not significant","Worst"))
  )

df_sankey_country <- df_strategy_country %>%
  count(strategy_E, strategy_F, name = "n")


p2 <- suppressWarnings(suppressMessages(
  ggplot(
    df_sankey_country,
    aes(axis1 = strategy_E, axis2 = strategy_F, y = n)
  ) +
    geom_alluvium(
      aes(fill = strategy_E),
      width = 0.22,
      alpha = 0.6,
      color = "grey70",
      linewidth = 0.3
    ) +
    geom_stratum(
      aes(fill = after_stat(stratum)),
      width = 0.22,
      color = "grey35",
      linewidth = 0.5
    ) +
    geom_text(
      stat = "stratum",
      aes(label = paste0(after_stat(stratum), "\n(n = ", after_stat(count), ")")),
      size = 3,
      fontface = "bold",
      lineheight = 0.95
    ) +
    scale_x_discrete(
      limits = c("Scenario E", "Scenario F"),
      expand = c(0.08, 0.08)
    ) +
    scale_fill_manual(
      values = pal,
      breaks = c("Worst", "Not significant", "Best")
    ) +
    labs(
      #title = "Strategy transitions from Scenario E to Scenario F",
      subtitle = paste0(
        "Stable category = ", percent(mean(df_strategy_country$strategy_E == df_strategy_country$strategy_F), accuracy = 0.1),
        "   |   Entered Best = ", sum(df_strategy_country$strategy_E != "Best" & df_strategy_country$strategy_F == "Best"),
        "   |   Exited Best = ", sum(df_strategy_country$strategy_E == "Best" & df_strategy_country$strategy_F != "Best")
      ),
      x = NULL,
      y = "Number of countries",
      fill = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "",
      panel.grid = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 11,hjust = 0.5,vjust = -3),
      axis.text.x = element_text(face = "bold"),
      legend.text = element_text(size = 11)
    )
))
p2

pw <- (p1 / p2) + plot_annotation(tag_levels = "A")
pw & theme(
  plot.tag = element_text(face = "bold"),
  axis.text = element_text(colour = "black", size = 9),
  axis.title = element_text(colour = "black", size = 10),
  plot.margin = margin(0, 1.5, 0, 1.5),
  panel.spacing = unit(0.5, "mm")  
)
ggsave("figures/supplement/Figure_S13.png",
       units = "cm",dpi = 300,
       width = 20, height = 20)