# Supplementary Figure S18 
# Sensitivity of effect estimates to imputation strategy


library(ggplot2)
library(ggtext)
data <- readRDS("output/body_size_sensitivity_analysis.rds")

data_1 <- data$basin

set2_colors <- c(
  "Description" = "#D55E00",
  "Geolocation" = "#009E73",
  "Sequencing"  = "#435792"
)

p1 <- ggplot(
  data = data_1,
  aes(
    x = delta_log_hr,
    y = variable,
    fill = shortfall,
    color = shortfall
  )
) +
  geom_vline(xintercept = 0, linetype = 2, size = 0.4, color = "gray30") +
  geom_pointrange(
    aes(xmin = .lower, xmax = .upper),
    position = position_dodge(0.2),
    linewidth = 0.5, shape = 21, stroke = 0.2, size = 0.46
  ) +
  scale_fill_manual(
    name = "First documentation",
    values = set2_colors,
    limits = c("Description", "Geolocation", "Sequencing")
  ) +
  scale_color_manual(
    name = "First documentation",
    values = set2_colors,
    limits = c("Description", "Geolocation", "Sequencing")
  ) +
  facet_grid(
    ~ realm
  ) +
  scale_x_continuous(limits = c(min(data_1$.lower),max(data_1$.upper)),breaks = c(-0.5,0,0.5))+
  scale_y_discrete(expand = expansion(c(0, 0))) +
  theme_bw()+
  theme(strip.background = element_blank(),  # Remove background for facet strips
        panel.background = element_rect(fill = "grey90"),
        # panel.grid.major  = element_blank(),
        legend.title = element_text(colour = "black", size = 9),  # Customize legend title
        legend.text = element_text(colour = "black", size = 8),  # Customize legend text
        legend.key.height = unit(0.5, "cm"),  # Adjust legend key height
        legend.key = element_rect(fill = NA, colour = NA),
        axis.title = element_text(colour = "black", size = 10),  # Customize axis titles
        axis.text.x = element_text(colour = "black", size = 8),  # Customize x-axis text
        axis.text.y = element_blank(),
        strip.text = element_markdown(size = 8)
  ) +  # Customize nested axis text appearance
  ylab("Body Size") +  # Remove y-axis label
  xlab("")

data_2 <- data$country


set2_colors <- c(
  "Description" = "#D55E00",
  "Geolocation" = "#009E73",
  "Sequencing"  = "#435792"
)


p2 <- ggplot(
  data = data_2,
  aes(
    x = delta_log_hr,
    y = variable,
    fill = shortfall,
    color = shortfall
  )
) +
  geom_vline(xintercept = 0, linetype = 2, size = 0.4, color = "gray30") +
  geom_pointrange(
    aes(xmin = .lower, xmax = .upper),
    position = position_dodge(0.2),
    linewidth = 0.5, shape = 21, stroke = 0.2, size = 0.46
  ) +
  scale_fill_manual(
    name = "First documentation",
    values = set2_colors,
    limits = c("Description", "Geolocation", "Sequencing")
  ) +
  scale_color_manual(
    name = "First documentation",
    values = set2_colors,
    limits = c("Description", "Geolocation", "Sequencing")
  ) +
  facet_grid(
    ~ continent
  ) +
  scale_x_continuous(limits = c(min(data_1$.lower),max(data_1$.upper)),breaks = c(-0.5,0,0.5))+
  scale_y_discrete(expand = expansion(c(0, 0))) +
  theme_bw()+
  theme(strip.background = element_blank(),  # Remove background for facet strips
        panel.background = element_rect(fill = "grey90"),
        # panel.grid.major  = element_blank(),
        legend.title = element_text(colour = "black", size = 9),  # Customize legend title
        legend.text = element_text(colour = "black", size = 8),  # Customize legend text
        legend.key.height = unit(0.5, "cm"),  # Adjust legend key height
        legend.key = element_rect(fill = NA, colour = NA),
        axis.title = element_text(colour = "black", size = 10),  # Customize axis titles
        axis.text.x = element_text(colour = "black", size = 8),  # Customize x-axis text
        axis.text.y = element_blank(),
        strip.text = element_markdown(size = 8)
  ) +  # Customize nested axis text appearance
  ylab("Body Size") +  # Remove y-axis label
  xlab(expression(Delta * log(HR) ~ "(95% CrI)"))

library(patchwork)  
p1 / p2 +
  plot_annotation(tag_levels = "A") +
  plot_layout(
    guides = "collect",    
    heights = c(1, 1)      
  ) &
  theme(
    legend.position = "top",        
    legend.justification = 0.5,     
    legend.box = "horizontal",      
    panel.spacing = unit(0.8, "pt"),       
    #plot.margin = unit(c(0,0,0,0), "pt"), 
    strip.margin = unit(0, "pt"),        
    legend.margin = margin(t=0, b=0, unit="pt"),
    axis.title = element_text(colour = "black", size = 10),
    axis.text.x = element_text(colour = "black", size = 7)
  )


ggsave("figures/supplement/Figure_S18.png",
       units = "cm",dpi = 300,
       width = 18, height = 12)
