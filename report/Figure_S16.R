# Supplementary Figure S16 
# Temporally independent validation of projected and subsequently observed species documentation across spatial units

rm(list = ls())
library(dplyr)
library(ggplot2)
data <- read.csv("output/tables/temporally_independent_validation.csv")
set2_colors <- c(
  "Description" = "#D55E00",
  "Geolocation" = "#009E73",
  "Sequencing"  = "#435792"
)

data <- data %>%
  mutate(
    sig_label = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

ggplot(data = data, aes(x = year, y = rho, colour = shortfall)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_line() +
  geom_text(
    aes(label = sig_label),
    vjust = 2,  
    size = 3,
    show.legend = FALSE
  ) +
  facet_grid(group ~ shortfall) +
  scale_colour_manual(values = set2_colors) +
  theme_bw() +
  labs(
    x = "Year",
    y = "Spearman correlation coefficient"
  ) +
  scale_x_continuous(breaks = seq(2015,2024,1)) +
  scale_y_continuous(breaks = seq(0,1,0.25), limits = c(0, 1)) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "gray95"),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, vjust = 0.6),
    axis.text = element_text(colour = "black", size = 9),
    axis.title = element_text(colour = "black", size = 10),
  )

ggsave("figures/supplement/Figure_S16.png",
       units = "cm",dpi = 300,
       width = 20, height = 10)
