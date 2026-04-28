# ------------------------------------------------------------------------------
# Supplementary Figure S16 
# Patterns of known species descriptions, spatial distribution, and phylogenetic coverage in global freshwater fishes
# ------------------------------------------------------------------------------

library(viridis)
library(rnaturalearth)
library(dplyr)
library(sf)
library(ggplot2)
library(tidyr)
library(stringr)
library(grid)  
library(patchwork)
library(minpack.lm) 
options(sf_use_s2 = FALSE)
data <- openxlsx::read.xlsx("input/raw/cas_freshwater_v1.xlsx")
################################################################################
# species discovery
sp <- data %>% dplyr::select(valid_name,year_description)
# Count species descriptions per year
species_per_year <- sp %>%
  count(year_description) %>%
  complete(year_description = 1758:2024, fill = list(n = 0)) %>%  # Fill missing years with 0
  arrange(year_description) %>%  # Ensure chronological order
  mutate(cumulative_species = cumsum(n), # Compute cumulative count
         cumulative_fraction = cumulative_species / max(cumulative_species) # Normalize to [0,1]
  )  
# Plot cumulative species descriptions over time
p1 <- ggplot(species_per_year, aes(x = year_description)) +
  geom_line(aes(y = cumulative_species),color = "#D55E00", linewidth = 0.5) +
  labs(#title = "Cumulative number of described species", 
    x = "Year", 
    y = "Cumulative number of described species") +
  theme_minimal()+
  theme(axis.title = element_text(colour = "black", face = "bold", size = 8),
        axis.text = element_text(colour = "black", size = 7),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = 2,linewidth = 0.5)
  )

# Plot annual species descriptions
p2 <- ggplot(species_per_year, aes(x = year_description, y = n)) +
  geom_col(fill = "#D55E00",width = 0.5) +
  #geom_smooth()+
  labs(#title = "Annual Number of Newly Described Species", 
    x = "Year", 
    y = "Annual number of newly-described species") +
  theme_minimal()+
  theme(axis.title = element_text(colour = "black", face = "bold", size = 8),
        axis.text = element_text(colour = "black", size = 7),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linetype = 2,linewidth = 0.5))

# combine
p1+p2
ggsave("figures/supplement/Figure_S16_A.png", width = 20, height = 8, units = "cm", dpi = 300)
################################################################################
df <- data %>%
  mutate(basin = str_split(basin, ";")) %>%        
  mutate(basin = lapply(basin, unique)) %>%        
  mutate(basin = lapply(basin, sort)) %>%
  unnest(basin) %>%
  group_by(basin) %>%
  summarise(count = n())

inland <- readRDS("input/raw/basin/basin_sf_v1.rds")
setdiff(unique(inland$basin),df$basin)
setdiff(df$basin,unique(inland$basin))


dt <- inland

world_map <- ne_countries(scale = 50, type = "countries", returnclass = "sf")


moll_proj <- st_crs("+proj=moll")


lat_points <- data.frame(lon = 0, lat = c(-60, 90)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)


projected_points <- st_transform(lat_points, crs = moll_proj)

y_limits <- st_coordinates(projected_points)[, "Y"]


dt$count_group <- cut(
  dt$n_species,
  breaks = c(1, 10, 25, 50, 100, 500, 1000,1500,Inf),  
  labels = c("<10","10-25","25-50" ,"50-100", "100-500","500-1000","1000-1500",">1500"),  
  include.lowest = TRUE
)

table(dt$count_group)
original_colors <- c("#B4E4D6", "#6DCFB0", "#009E73", "#006D4F", "#003D2C")
custom_colors <- colorRampPalette(original_colors)(8)
dt_fixed <- st_wrap_dateline(dt, options = c("WRAPDATELINE=YES"))

p3 <- ggplot(data = dt_fixed) +
  #geom_sf(data = world_map, fill = "lightgrey", colour = NA) +  # 背景地图
  #geom_sf(aes(fill = count_group), colour = NA) +  # 按分段填充颜色
  ggrastr::rasterise(geom_sf(data = world_map, fill = "lightgrey", colour = NA) ,dpi = 300)+
  ggrastr::rasterise(geom_sf(aes(fill = count_group), colour = NA) ,dpi = 300)+
  scale_fill_manual(
    values = custom_colors,  
    name = "Number of freshwater fish per drainage basin",  
    na.translate = FALSE  
  ) +
  theme_void() +  
  theme(
    panel.border = element_blank(),  
    legend.position = "bottom",  
    #legend.position = c(0.6,0.1),
    legend.text = element_text(size = 5),
    legend.title = element_text(hjust = 0.5,face = "bold",size = 6),  
    legend.box = "horizontal"  
  ) +
  guides(
    fill = guide_legend(
      title.position = "top",  
      title.hjust = 0.5,  
      label.position = "bottom",  
      keywidth = 1.3,  
      keyheight = 0.4,  
      nrow = 1  
    )
  ) +
  coord_sf(crs = "+proj=moll +lon_0=0", ylim = y_limits, expand = FALSE)

ggsave("figures/supplement/Figure_S16_B.png", width = 12, height = 10, units = "cm", dpi = 300)
################################################################################
library(arrow)
tre <- readRDS("input/raw/fish_tree_final.rds")
seq_annotation <- read_parquet("input/raw/seq_annotation_v1.parquet")
length(unique(seq_annotation$valid_name))
# 10756(57%)
tre$tip.label <- gsub("[*]*$", "", tre$tip.label)

meta_data <- data.frame(
  label = gsub("[*]*$", "", tre$tip.label)
)
meta_data$Status <- NA
meta_data$Status <- ifelse(gsub("_"," ",meta_data$label) %in% unique(seq_annotation$valid_name),"exist","grafted")

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("ggtree")
library(ggstar)
library(ggtree)
library(ggtreeExtra)  # geom_fruit() #BiocManager::install("ggtreeExtra")

tree_df <- full_join(as_tibble(tre), meta_data, by = c("label" = "label"))
p <- ggtree(tre,layout = "fan", open.angle = 180,linewidth=0.15,col = "lightgrey") +
  geom_rootedge(rootedge = 40,col = "lightgrey",linewidth = 0.15)

tree_df <- p$data %>% 
  left_join(meta_data, by = c("label" = "label"))

pp <- p + 
  geom_tree(data = tree_df, aes(color = Status), size = 0.15) +
  scale_color_manual(values = c("exist" = "#435792", "grafted" = "#8A2BE2"),na.value = NA)

p4 <- pp + 
  geom_fruit(
    data = meta_data, 
    #geom = geom_point,  
    geom = geom_star,
    mapping = aes(y = label, colour = Status),
    starshape = 28, size = 0.05
  ) +
  scale_colour_manual("Sequence (%)",
                      values = c("exist" = "#435792", "grafted" = "#8A2BE2"),na.value = NA,
                      labels = c("exist" = "Available (57%)", "grafted" = "Unavailable (43%)")
  ) +
  theme(
    legend.position = "none",
    legend.title = element_text(face = "bold", size = 6, hjust = 0, vjust = 0),
    legend.text = element_text(size = 5),
    legend.spacing.y = unit(0, "cm"),
    legend.background = element_blank()
  ) +
  guides(colour = guide_legend(
    override.aes = list(size = 3, shape = 15)
  ))

p4

ggsave("figures/supplement/Figure_S16_C.png", width = 10, height = 12, units = "cm", dpi = 300)

ggplot(data = data.frame(type = c("Available","Unavailable"),
                         n = c(57,43)))+
  geom_col(aes(type,n, fill = type),show.legend = F, width = 0.1)+
  geom_text(aes(type,n+5, label = paste0(n,"%")),size = 5)+
  scale_fill_manual(values = c("Available" = "#435792", "Unavailable" = "#8A2BE2"))+
  coord_flip()+
  scale_x_discrete(expand = c(0.1,0.1))+
  scale_y_continuous(limits = c(0,100),expand = c(0,0))+
  xlab("")+
  ylab("Percentage (%)")+
  theme_classic()+
  theme(axis.line.y = element_blank(),
        axis.line.x = element_line(colour = "black",linewidth = 0.2),
        axis.ticks = element_line(colour = "black",linewidth = 0.2),
        axis.text = element_text(face = "bold",size = 10,colour = "black"),
        axis.title = element_text(face = "bold",size = 10,colour = "black"))
