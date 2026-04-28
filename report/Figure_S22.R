# Supplementary Figure S22 
# Sensitivity of the best-to-worst ratio (BWR) to the equity discount parameter θ

library(ggplot2)
data <- readRDS("output/equity-adjusted_prioritization.rds")


p1 <- ggplot(data$basin, aes(x = discount, y = BWR)) +
  geom_point(size = 3, color = "#2E86AB") +
  geom_linerange(
    aes(ymin = BWR_CI_95_lower, ymax = BWR_CI_95_upper),
    color = "#2E86AB", linewidth = 0.8
  ) +
  # geom_text(
  #   aes(label = sprintf("%.2f", BWR)),
  #   hjust = -0.5, size = 3.8, fontface = "bold"
  # ) +
 geom_text(
  aes(label = paste0("n[best]==", n_best)),
  hjust = -0.18,  
  vjust = -1,
  size = 3.2,
  color = "gray30",
  parse = TRUE,
  angle = 90
) +
 geom_text(
  aes(label = paste0("n[worst]==", n_worst)),
  hjust = -0.16, 
  vjust = 1.5,
  size = 3.2,
  color = "gray30",
  parse = TRUE,
  angle = 90
) +
  scale_x_continuous(breaks = data$basin$discount) +
  scale_y_continuous(limits = c(0,8))+
  labs(
    x = expression(paste("Discount parameter ", theta)),
    y = "BWR (95% CI)"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

p2 <- ggplot(data$country, aes(x = discount, y = BWR)) +
  geom_point(size = 3, color = "#2E86AB") +
  geom_linerange(
    aes(ymin = BWR_CI_95_lower, ymax = BWR_CI_95_upper),
    color = "#2E86AB", linewidth = 0.8
  ) +
 
  # geom_text(
  #   aes(label = sprintf("%.2f", BWR)),
  #   hjust = -0.5, size = 3.8, fontface = "bold"
  # ) +
  
geom_text(
  aes(label = paste0("n[best]==", n_best)),
  hjust = -0.5,  
  vjust = -1,
  size = 3.2,
  color = "gray30",
  parse = TRUE,
  angle = 90
) +
 
geom_text(
  aes(label = paste0("n[worst]==", n_worst)),
  hjust = -0.48, 
  vjust = 1.5,
  size = 3.2,
  color = "gray30",
  parse = TRUE,
  angle = 90
) +
  scale_x_continuous(breaks = data$country$discount) +
  scale_y_continuous(limits = c(0,4))+
  labs(
    x = expression(paste("Discount parameter ", theta)),
    y = "BWR (95% CI)"
  ) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

library(patchwork)  
p1 / p2 +
  plot_annotation(tag_levels = "A") +
  plot_layout(
    guides = "collect",    
    heights = c(1, 1)      
  ) &
  theme(
    plot.tag = element_text(colour = "black", face = "bold", size = 12),
    axis.title = element_text(colour = "black", size = 10),
    axis.text = element_text(colour = "black", size = 9)
  )

ggsave("figures/supplement/Figure_S22.png",
       units = "cm",dpi = 300,
       width = 18, height = 12)
