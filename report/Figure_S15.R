# Supplementary Figure S15 
# Spatial correspondence between observed and predicted withheld species richness across validation datasets

rm(list = ls())
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggplot2)

basin <- read.csv("output/tables/basin_test_data_observed_predicted.csv")
country <- read.csv("output/tables/contry_test_data_observed_predicted.csv")


set2_colors <- c(
  "Description" = "#D55E00",
  "Geolocation" = "#009E73",
  "Sequencing"  = "#435792"
)
cor_basin <- basin %>%
  group_by(withheld, shortfall) %>%
  summarise(
    n = n(),
    ct = list(suppressWarnings(
      cor.test(observed_withheld_SR, predicted_withheld_SR, method = "spearman")
    )),
    r = round(ct[[1]]$estimate, 3),
    p_val = ct[[1]]$p.value,
    label = ifelse(
      p_val < 0.001,
      paste0("R==", r, "~~italic(P) < 0.001~~n==", n),
      paste0("R==", r, "~~italic(P)==", round(p_val, 3), "~~n==", n)
    ),
    .groups = "drop"
  )

p1 <- ggplot(data = basin,
             aes(x = observed_withheld_SR, y = predicted_withheld_SR, colour = shortfall)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.2) +
  geom_text(
    data = cor_basin,
    aes(label = label, colour = shortfall),
    x = Inf, y = 0,
    hjust = 1.1, vjust = 0,
    size = 2.5,
    parse = TRUE,
    inherit.aes = FALSE
  ) +
  scale_colour_manual(values = set2_colors) +
  theme_bw() +
  facet_grid(withheld ~ shortfall, scales = "free", space = "free_y") +
  labs(
    x = "Observed richness per basin",
    y = "Predicted richness per basin"
  ) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "gray95"),
    legend.position = "none"
  )

cor_country <- country %>%
  group_by(withheld, shortfall) %>%
  summarise(
    n = n(),
    ct = list(suppressWarnings(
      cor.test(observed_withheld_SR, predicted_withheld_SR, method = "spearman")
    )),
    r = round(ct[[1]]$estimate, 3),
    p_val = ct[[1]]$p.value,
    label = ifelse(
      p_val < 0.001,
      paste0("R==", r, "~~italic(P) < 0.001~~n==", n),
      paste0("R==", r, "~~italic(P)==", round(p_val, 3), "~~n==", n)
    ),
    .groups = "drop"
  )

p2 <- ggplot(data = country,
             aes(x = observed_withheld_SR, y = predicted_withheld_SR, colour = shortfall)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.2) +
  geom_text(
    data = cor_country,
    aes(label = label, colour = shortfall),
    x = Inf, y = 0,
    hjust = 1.1, vjust = 0,
    size = 2.5,
    parse = TRUE,
    inherit.aes = FALSE
  ) +
  scale_colour_manual(values = set2_colors) +
  theme_bw() +
  facet_grid(withheld ~ shortfall, scales = "free", space = "free_y") +
  labs(
    x = "Observed richness per country",
    y = "Predicted richness per country"
  ) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "gray95"),
    legend.position = "none"
  )


pw <- (p1 / p2) + plot_annotation(tag_levels = "A")
pw & theme(
  plot.tag = element_text(face = "bold"),
  axis.text = element_text(colour = "black", size = 9),
  axis.title = element_text(colour = "black", size = 10),
  plot.margin = margin(0, 1.5, 0, 1.5),
  panel.spacing = unit(0.5, "mm")  
)
ggsave("figures/supplement/Figure_S15.png",
       units = "cm",dpi = 300,
       width = 20, height = 20)
