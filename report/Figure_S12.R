# ------------------------------------------------------------------------------
# Supplementary Figure S10 Prioritization of freshwater fish biodiversity collections at the national scale based on multi-objective trade-offs
# ------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(GGally)
library(rPref)
library(patchwork)
library(legendry)
library(ggh4x)
library(cowplot)

cost <- readRDS("input/processed/cost_country_all.rds")
sdg <- readRDS("input/processed/sdg_country.rds")
shortfall <- readRDS("input/processed/country_shortfall.rds")
country_list <- read.csv("input/raw/country_list.csv")
country_list <- country_list %>% 
  mutate(continent = recode(continent,
                            "Africa"   = "AFR",
                            "Asia"     = "ASI",
                            "Europe"  = "EUR",
                            "North America" = "NAM",
                            "South America" = "SAM",
                            "Oceania" = "OCE"
  ))

sdg.pca <- prcomp(sdg[,2:18], center = TRUE, scale. = TRUE)
sdg.pca$rotation
# PC1(Socioeconomics)# ensure higher value = stronger socioeconomic capacity
# Governance = PC2  # already correct direction
sdg$PC1 <- -sdg.pca$x[, 1] 
sdg$PC2 <- sdg.pca$x[, 2] 

names(sdg)[21:22] <- c("Socioeconomics","Governance")

data <- shortfall %>% select(iso3 = iso3, Shortfall = NSCI_norm) %>%
  left_join(cost %>% select(iso3 = iso3, Investment = cost_scale), by = "iso3") %>%
  left_join(sdg %>% select(iso3,Socioeconomics,Governance), by = "iso3") %>%
  left_join(country_list[,c(1,5)],by = "iso3") 

country_shortfalls <- read.csv("output/tables/country_shortfalls.csv") %>%
  select(iso3,scenario) %>%
  na.omit() %>%
  pull(iso3)

# data <- data %>% filter(basin_id %in% basin_shortfalls)
##################################################################################
# linnaean <- read.csv("output/tables/basin_linnaean_shortfall.csv")
# wallacean <- read.csv("output/tables/basin_wallacean_shortfall.csv")
# darwinian <- read.csv("output/tables/basin_darwinian_shortfall.csv")
# 
# df_basin <- linnaean[,c(1,4)] %>%
#   left_join(wallacean[,1:2], by = "basin_id") %>%
#   left_join(darwinian[,1:2], by = "basin_id")
# 
# df_nsc <- df_basin %>%
#   mutate(
#     SL = min_max_normalize_safe(prob_undesc_mean),
#     SW = min_max_normalize_safe(prob_nongeoloc_mean),
#     SD = min_max_normalize_safe(prob_noseq_mean),
# 
#     # nested adjustment
#     SW_adj = (1 - SL) * SW,
#     SD_adj = (1 - SL) * SD,
# 
#     # rank composite index
#     NSCI = min_max_normalize_safe(SL + SW_adj + SD_adj)
#   )

# df_nsc <- shortfall
# hist(df_nsc$NSCI)
# 
# df_nsc <- df_nsc %>% left_join(country_list[,c(1,5)],by = "iso3")
# 
# rr <- df_nsc %>% select(iso3,NSCI,continent) %>% na.omit()
# # Step 1: Realm-based expectation model
# model <- lm(NSCI ~ continent, data = rr)
# 
# # Step 2: Realm-adjusted residuals
# rr$resid_NSCI <- residuals(model)
# 
# # Step 3: Compute μ and σ of residuals
# mu_resid <- mean(rr$resid_NSCI, na.rm = TRUE)
# sd_resid <- sd(rr$resid_NSCI, na.rm = TRUE)
# 
# # Step 4: Cutoff values consistent with the published method
# z <- 1.96
# upper_cut <- mu_resid + z * sd_resid
# lower_cut <- mu_resid - z * sd_resid
# 
# # Step 5: Classification
# rr$NSCI_high   <- rr$resid_NSCI >= upper_cut
# rr$NSCI_low    <- rr$resid_NSCI <= lower_cut
# rr$NSCI_normal <- rr$resid_NSCI > lower_cut & rr$resid_NSCI < upper_cut
# id <- rr %>% filter(NSCI_high == TRUE) %>% pull(iso3)
# 
# unique(country_shortfalls, id)


data <- data %>% filter(iso3 %in% country_shortfalls)

min_max_normalize_safe <- function(x, epsilon = 1e-6) {
  min_val <- min(x,na.rm = T)
  max_val <- max(x,na.rm = T)
  range_val <- max_val - min_val + epsilon
  return((x - min_val) / range_val)
}

data$Investment <- min_max_normalize_safe(data$Investment)
data$Shortfall <- min_max_normalize_safe(data$Shortfall)
data$Socioeconomics <- min_max_normalize_safe(data$Socioeconomics)
data$Governance <- min_max_normalize_safe(data$Governance)
rm(cost,sdg,sdg.pca,shortfall,country_list,country_shortfalls,min_max_normalize_safe)
################################################################################
# Realm-specific Pareto frontier

#High Shortfall × High Socioeconomics（Best）
#High Shortfall × Low Socioeconomics（Worst）
# Figure X. Trade-off between knowledge shortfall and socioeconomic feasibility, with Pareto-optimal strategic priorities highlighted.


library(dplyr)
library(rPref)
library(emoa) #Pareto dominance(minimization) emoa::nds_rank()
library(ggplot2)
library(ggrepel)
# ------------------------------------------------------------
# Prepare data and specify objective directions
#    Shortfall  : higher is better  → multiply by -1
#    Investment : lower is better  → keep unchanged
# ------------------------------------------------------------
df1 <- data %>%
  select(iso3,continent, Shortfall, Investment) %>%
  filter(!is.na(Investment), !is.na(Shortfall)) %>%
  mutate(
    # Convert objectives to minimization form for nds_rank():
    # High Shortfall (good) → -Shortfall (so smaller = better)
    Shortfall_obj  = -Shortfall,
    Investment_obj = Investment # low Investment (good) → keep as original value
  )

df_plot_realm_1 <- df1 %>%
  group_by(continent) %>%
  group_modify(~ {
    obj_mat <- as.matrix(.x[, c("Shortfall_obj", "Investment_obj")])
    fronts  <- emoa::nds_rank(t(obj_mat))
    stopifnot(length(fronts) == nrow(.x))
    .x %>%
      mutate(front = fronts,
             front_rank = 1 - percent_rank(front),# 1 = best (front=1), 0 = worst (front=max)
             front_scaled = rank(-front, ties.method = "min")
      )   
  }) 

q <- 0.3
df_plot_realm_1 <- df_plot_realm_1 %>% 
  as.data.frame() %>%
  mutate(
    strategy = NULL,
    strategy = case_when(
      front_scaled >= quantile(front_scaled, 1 - q, na.rm = TRUE) ~ "best",
      front_scaled <= quantile(front_scaled, q,     na.rm = TRUE) ~ "worst",
      TRUE ~ "not_significant"
    )
  ) 


p1 <- ggplot() +
  geom_point(data = df_plot_realm_1,aes(x = Investment,y = Shortfall,colour = front_rank),size  = 1,alpha = 0.5) +
  geom_point(data = df_plot_realm_1 %>% group_by(continent) %>%
               group_modify(~ {psel(.x, high(Shortfall) * low(Investment))}) %>% ungroup(),
             aes(x = Investment, y = Shortfall,group = continent),shape  = 21,fill   = "#00B050",colour = "black",size = 1,stroke = 0.3) +
  # geom_path(data = df_plot_realm_1 %>% group_by(continent) %>% 
  #             group_modify(~ {psel(.x, high(Shortfall) * low(Investment)) %>% arrange(Investment)}) %>% ungroup(),
  #           aes(x = Investment, y = Shortfall, lty = continent),linewidth = 0.2,colour = "#00B050",show.legend = FALSE, arrow = arrow(type = "closed", length = unit(0.2, "cm"))
  # )+
  geom_text_repel(data = df_plot_realm_1,
                  aes(x = Investment, y = Shortfall,label = iso3), nudge_x = 0.03,direction = "x", hjust = 0, segment.color = NA, size = 2,colour = "black") +
  ####
  geom_point(data = df_plot_realm_1 %>% group_by(continent) %>%
               group_modify(~ {psel(.x, low(Shortfall) * high(Investment))}) %>% ungroup(),
             aes(x = Investment, y = Shortfall,group = continent),shape  = 21,fill= "#FFC000",colour = "black",size = 1,stroke = 0.3) +
  # geom_path(data = df_plot_realm_1 %>% group_by(continent) %>% 
  #             group_modify(~ {psel(.x, low(Shortfall) * high(Investment)) %>% arrange(Investment)}) %>% ungroup(),
  #           aes(x = Investment, y = Shortfall, lty = continent),linewidth = 0.2,colour = "#FFC000",show.legend = FALSE, arrow = arrow(type = "closed", length = unit(0.2, "cm"))
  # )+
  # geom_text_repel( data = df_plot_realm_1 %>% group_by(continent) %>%
  #                    group_modify(~ {psel(.x, low(Shortfall) * high(Investment))}) %>% ungroup() %>% group_by(continent) %>% 
  #                    filter(Investment == min(Investment)) %>% slice(1) %>% ungroup(),
  #                  aes(x = Investment, y = Shortfall,label = continent), nudge_x = 0.02,nudge_y = 0.03,direction = "y", hjust = 0, segment.color = NA, size = 2,colour = "black") +
  scale_colour_gradientn(name = "Pareto frontier",colours = c("#FFC000", "#BFBFBF", "#00B050"),
                         limits = c(0, 1),breaks = c(0, 1), labels = c("Worst","Best"),
                         guide = guide_colbar(show = c(TRUE, TRUE), oob = "squish"))+
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.1) +
  labs(x = "Investment", y = "Shortfall (NSCI)") +
  theme_classic() +
  theme(
    axis.text  = element_text(colour = "black", size = 7),
    axis.title = element_text(colour = "black", size = 8),
    legend.title = element_text(size = 6, hjust = 0.5),
    legend.text  = element_text(size = 5),
    legend.key.width = unit(0.3,"cm"),
    legend.key.height = unit(0.15,"cm"),
    legend.background = element_blank(),
    legend.position = c(0.18,0.90),
    legend.title.position = "top",
    legend.direction = "horizontal",
    plot.margin = margin(0,0,0,0)
  )
################################################################################
# Prepare data: Shortfall + Socioeconomics
df2 <- data %>%
  select(iso3, continent,Shortfall, Socioeconomics) %>%       
  filter(!is.na(Shortfall), !is.na(Socioeconomics)) %>%
  mutate(
    # Both objectives need to be maximized,
    # so convert them into minimization for nds_rank():
    Shortfall_obj = -Shortfall, # higher shortfall = better priority
    Socio_obj = -Socioeconomics # higher socioeconomic capacity = better
  )

df_plot_realm_2 <- df2 %>%
  group_by(continent) %>%
  group_modify(~ {
    obj_mat <- as.matrix(.x[, c("Shortfall_obj", "Socio_obj")])
    fronts <- emoa::nds_rank(t(obj_mat))
    stopifnot(length(fronts) == nrow(.x))
    .x %>%
      mutate(front= fronts,
             front_rank = 1 - percent_rank(front),# 1 = best (front=1), 0 = worst (front=max)
             front_scaled = rank(-front, ties.method = "min")
      )
  }) 


q <- 0.2
df_plot_realm_2 <- df_plot_realm_2 %>% 
  as.data.frame() %>%
  mutate(
    strategy = NULL,
    strategy = case_when(
      front_scaled >= quantile(front_scaled, 1 - q, na.rm = TRUE) ~ "best",
      front_scaled <= quantile(front_scaled, q,     na.rm = TRUE) ~ "worst",
      TRUE ~ "not_significant"
    )
  ) 




p2 <- ggplot() +
  geom_point(data = df_plot_realm_2,aes(x = Socioeconomics,y = Shortfall,colour = front_rank),size  = 1,alpha = 0.5) +
  geom_point(data = df_plot_realm_2 %>% group_by(continent) %>%
               group_modify(~ {psel(.x, high(Shortfall) * high(Socioeconomics))}) %>% ungroup(),
             aes(x = Socioeconomics, y = Shortfall,group = continent),shape  = 21,fill   = "#00B050",colour = "black",size = 1,stroke = 0.3) +
  # geom_path(data = df_plot_realm_2 %>% group_by(continent) %>% 
  #             group_modify(~ {psel(.x, high(Shortfall) * high(Socioeconomics)) %>% arrange(Socioeconomics)}) %>% ungroup(),
  #           aes(x = Socioeconomics, y = Shortfall, lty = continent),linewidth = 0.2,colour = "#00B050",show.legend = FALSE, arrow = arrow(type = "closed", length = unit(0.2, "cm"))
  # )+
  geom_text_repel( data = df_plot_realm_2,
                   aes(x = Socioeconomics, y = Shortfall,label = iso3), nudge_x = 0.03,direction = "x", hjust = 0, segment.color = NA, size = 2,colour = "black") +
  ####
  geom_point(data = df_plot_realm_2 %>% group_by(continent) %>%
               group_modify(~ {psel(.x, low(Shortfall) * low(Socioeconomics))}) %>% ungroup(),
             aes(x = Socioeconomics, y = Shortfall,group = continent),shape  = 21,fill= "#FFC000",colour = "black",size = 1,stroke = 0.3) +
  # geom_path(data = df_plot_realm_2 %>% group_by(continent) %>% 
  #             group_modify(~ {psel(.x, low(Shortfall) * low(Socioeconomics)) %>% arrange(Socioeconomics)}) %>% ungroup(),
  #           aes(x = Socioeconomics, y = Shortfall, lty = continent),linewidth = 0.2,colour = "#FFC000",show.legend = FALSE, arrow = arrow(type = "closed", length = unit(0.2, "cm"))
  # )+
  # geom_text_repel( data = df_plot_realm_2 %>% group_by(continent) %>%
  #                    group_modify(~ {psel(.x, low(Shortfall) * low(Socioeconomics))}) %>% ungroup() %>% group_by(continent) %>% 
  #                    filter(Socioeconomics == min(Socioeconomics)) %>% slice(1) %>% ungroup(),
  #                  aes(x = Socioeconomics, y = Shortfall,label = continent), nudge_x = 0.03,direction = "x", hjust = 0, segment.color = NA, size = 2,colour = "black") +
  scale_colour_gradientn(name = "Pareto frontier",colours = c("#FFC000", "#BFBFBF", "#00B050"),
                         limits = c(0, 1),breaks = c(0, 1), labels = c("Worst","Best"),
                         guide = guide_colbar(show = c(TRUE, TRUE), oob = "squish"))+
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.1) +
  labs(x = "Socioeconomics", y = "Shortfall (NSCI)") +
  theme_classic() +
  theme(
    axis.text  = element_text(colour = "black", size = 7),
    axis.title = element_text(colour = "black", size = 8),
    legend.title = element_text(size = 8, hjust = 0.5),
    legend.text  = element_text(size = 6),
    legend.key.width = unit(0.4,"cm"),
    legend.key.height = unit(0.2,"cm"),
    legend.background = element_blank(),
    legend.position = "none",
    legend.title.position = "top",
    legend.direction = "horizontal",
    plot.margin = margin(0,0,0,0)
  )


################################################################################

#    Goal: prioritize basins with
#      - high Shortfall   (large knowledge gaps)
#      - high Governance  (strong capacity to fill those gaps)
# ------------------------------------------------------------
df3 <- data %>%
  select(iso3, continent,Governance, Shortfall) %>%
  filter(!is.na(Governance), !is.na(Shortfall)) %>%
  mutate(
    # Convert both objectives to "minimization" form for nds_rank():
    # We want to MAXIMIZE Shortfall → use -Shortfall
    Shortfall_obj  = -Shortfall,
    # We want to MAXIMIZE Governance → use -Governance
    Governance_obj = -Governance
  )

df_plot_realm_3 <- df3 %>%
  group_by(continent) %>%
  group_modify(~ {
    obj_mat <- as.matrix(.x[, c("Shortfall_obj", "Governance_obj")])
    fronts <- emoa::nds_rank(t(obj_mat))
    stopifnot(length(fronts) == nrow(.x))
    .x %>%
      mutate(front = fronts,
             front_rank = 1 - percent_rank(front),# 1 = best (front=1), 0 = worst (front=max)
             front_scaled = rank(-front, ties.method = "min")
      )
  }) 

q <- 0.2
df_plot_realm_3 <- df_plot_realm_3 %>% 
  as.data.frame() %>%
  mutate(
    strategy = NULL,
    strategy = case_when(
      front_scaled >= quantile(front_scaled, 1 - q, na.rm = TRUE) ~ "best",
      front_scaled <= quantile(front_scaled, q,     na.rm = TRUE) ~ "worst",
      TRUE ~ "not_significant"
    )
  ) 



p3 <- ggplot() +
  geom_point(data = df_plot_realm_3,aes(x = Governance,y = Shortfall,colour = front_rank),size  = 1,alpha = 0.5) +
  geom_point(data = df_plot_realm_3 %>% group_by(continent) %>%
               group_modify(~ {psel(.x, high(Shortfall) * high(Governance))}) %>% ungroup(),
             aes(x = Governance, y = Shortfall,group = continent),shape  = 21,fill   = "#00B050",colour = "black",size = 1,stroke = 0.3) +
  # geom_path(data = df_plot_realm_3 %>% group_by(continent) %>% 
  #             group_modify(~ {psel(.x, high(Shortfall) * high(Governance)) %>% arrange(Governance)}) %>% ungroup(),
  #           aes(x = Governance, y = Shortfall, lty = continent),linewidth = 0.2,colour = "#00B050",show.legend = FALSE, arrow = arrow(type = "closed", length = unit(0.2, "cm"))
  # )+
  geom_text_repel( data = df_plot_realm_3,
                   aes(x = Governance, y = Shortfall,label = iso3), nudge_x = 0.03,direction = "x", hjust = 0, segment.color = NA, size = 2,colour = "black") +
  ####
  geom_point(data = df_plot_realm_3 %>% group_by(continent) %>%
               group_modify(~ {psel(.x, low(Shortfall) * low(Governance))}) %>% ungroup(),
             aes(x = Governance, y = Shortfall,group = continent),shape  = 21,fill= "#FFC000",colour = "black",size = 1,stroke = 0.3) +
  # geom_path(data = df_plot_realm_3 %>% group_by(continent) %>% 
  #             group_modify(~ {psel(.x, low(Shortfall) * low(Governance)) %>% arrange(Governance)}) %>% ungroup(),
  #           aes(x = Governance, y = Shortfall, lty = continent),linewidth = 0.2,colour = "#FFC000",show.legend = FALSE, arrow = arrow(type = "closed", length = unit(0.2, "cm"))
  # )+
  # geom_text_repel( data = df_plot_realm_3 %>% group_by(continent) %>%
  #                    group_modify(~ {psel(.x, low(Shortfall) * low(Governance))}) %>% ungroup() %>% group_by(continent) %>% 
  #                    filter(Governance == min(Governance)) %>% slice(1) %>% ungroup(),
  #                  aes(x = Governance, y = Shortfall,label = continent), nudge_x = 0.03,direction = "x", hjust = 0, segment.color = NA, size = 2,colour = "black") +
  scale_colour_gradientn(name = "Pareto frontier",colours = c("#FFC000", "#BFBFBF", "#00B050"),
                         limits = c(0, 1),breaks = c(0, 1), labels = c("Worst","Best"),
                         guide = guide_colbar(show = c(TRUE, TRUE), oob = "squish"))+
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.1) +
  geom_hline(yintercept = 0.5, linetype = "dashed", linewidth = 0.1) +
  labs(x = "Governance", y = "Shortfall (NSCI)") +
  theme_classic() +
  theme(
    axis.text  = element_text(colour = "black", size = 7),
    axis.title = element_text(colour = "black", size = 8),
    legend.title = element_text(size = 8, hjust = 0.5),
    legend.text  = element_text(size = 6),
    legend.key.width = unit(0.4,"cm"),
    legend.key.height = unit(0.2,"cm"),
    legend.background = element_blank(),
    legend.position = "none",
    legend.title.position = "top",
    legend.direction = "horizontal",
    plot.margin = margin(0,0,0,0)
  )

################################################################################
# ------------------------------------------------------------
#  Prepare data for 4-objective trade-off:
#    Goal: High Shortfall × Low Investment × High Socioeconomics × High Governance
#
#    Assumptions:
#    - Shortfall:       higher = larger knowledge gaps (priority ↑)
#    - Investment:      higher = higher total cost (priority ↓)
#    - Socioeconomics:  higher = stronger socioeconomic capacity
#    - Governance:      higher = stronger institutional/governance capacity
# ------------------------------------------------------------
df4 <- data %>%
  select(iso3, continent, Investment, Socioeconomics, Governance, Shortfall) %>%
  na.omit() %>%
  mutate(
    # Convert all objectives to "minimize" form for nds_rank():
    # We want to MAXIMIZE Shortfall → minimize -Shortfall
    Shortfall_obj      = -Shortfall,
    # We want to MINIMIZE Investment → keep as is
    Investment_obj     =  Investment,
    # We want to MAXIMIZE Socioeconomics → minimize -Socioeconomics
    Socioeconomics_obj = -Socioeconomics,
    # We want to MAXIMIZE Governance → minimize -Governance
    Governance_obj     = -Governance
  )

df_plot_realm_4 <- df4 %>%
  group_by(continent) %>%
  group_modify(~ {
    obj_mat <- as.matrix(.x[, c(
      "Shortfall_obj",
      "Investment_obj",
      "Socioeconomics_obj",
      "Governance_obj"
    )])
    fronts <- emoa::nds_rank(t(obj_mat))
    stopifnot(length(fronts) == nrow(.x))
    .x %>%
      mutate(front = fronts,
             front_rank = 1 - percent_rank(front),# 1 = best (front=1), 0 = worst (front=max)
             front_scaled = rank(-front, ties.method = "min")
      )
  }) 


q <- 0.2
df_plot_realm_4 <- df_plot_realm_4 %>% 
  as.data.frame() %>%
  group_by(continent) %>%
  mutate(
    strategy = NULL,
    strategy = case_when(
      front_scaled >= quantile(front_scaled, 1 - q, na.rm = TRUE) ~ "best",
      front_scaled <= quantile(front_scaled, q,     na.rm = TRUE) ~ "worst",
      TRUE ~ "not_significant"
    )
  ) 



p4 <- ggparcoord(df_plot_realm_4,columns = c(7,8,9,10), groupColumn = 12,
                 scale = "uniminmax", showPoints = F,
                 alphaLines = 0.5)+ 
  scale_colour_gradientn(colours = c("#FFC000", "grey90","#00B050"))+
  #geom_point(size = 1)+
  #geom_path(linewidth = 0.02,alpha = 0.5)+
  #annotate(geom = "text",x = 3.6,y = 1.05,label = "BWR = 34", size = 2)+
  theme_minimal()+
  #coord_fixed(ratio = 3) +
  scale_x_discrete(expand = c(0.02,0.05),labels = c("Shortfall \n (NSCI)","Investment","Socioeconomics","Governance"))+
  xlab("")+
  ylab("")+
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey50", linewidth = 1),
        legend.position = "none",
        axis.text.x = element_text(colour = "black",size = 6.5),
        axis.text.y = element_text(colour = "black",size = 7),
        plot.margin = margin(0,0,0,0)
  )
################################################################################
################################################################################
linnaean <- read.csv("output/tables/country_linnaean_shortfall.csv")
wallacean <- read.csv("output/tables/country_wallacean_shortfall.csv")
darwinian <- read.csv("output/tables/country_darwinian_shortfall.csv")


front <- rbind(df_plot_realm_1 %>% select(continent,iso3,front_scaled,strategy) %>% mutate(scenario = "scenario_1"),
               df_plot_realm_2 %>% select(continent,iso3,front_scaled,strategy) %>% mutate(scenario = "scenario_2"),
               df_plot_realm_3 %>% select(continent,iso3,front_scaled,strategy) %>% mutate(scenario = "scenario_3"),
               df_plot_realm_4 %>% select(continent,iso3,front_scaled,strategy) %>% mutate(scenario = "scenario_4")
)
shortfall <- readRDS("input/processed/country_shortfall.rds")
dt <- front %>% left_join(shortfall[,c(1,7:9)],by = "iso3")
names(dt)[6:8] <- c("LS","WS","DS")
min_max_normalize_safe <- function(x, epsilon = 1e-3) {
  min_val <- min(x,na.rm = T)
  max_val <- max(x,na.rm = T)
  range_val <- max_val - min_val + epsilon
  return((x - min_val) / range_val)
}
dt$LS <- min_max_normalize_safe(dt$LS)
dt$WS <- min_max_normalize_safe(dt$WS)
dt$DS <- min_max_normalize_safe(dt$DS)
head(dt)
dt <- dt %>% left_join(data[,1:5],by = "iso3")


library(dplyr)
library(purrr)
library(tidyr)
get_cost <- function(df, scenario) {
  if (scenario == "scenario_1") {
    return(df$Investment)
  } else if (scenario == "scenario_2") {
    return(df$Socioeconomics)
  } else if (scenario == "scenario_3") {
    return(df$Governance)
  } else if (scenario == "scenario_4") {
    return(df$Investment * df$Socioeconomics * df$Governance)
  } else {
    stop("Unknown scenario")
  }
}

#' Compute Best-Worst Ratio of Efficiency (BWR)
#' 
#' This function calculates the ratio of mean efficiency between `best` and `worst` strategy groups,
#' along with key statistical metrics (sample size, standard error, coefficient of variation, 95% CI)
#' required for scientific reporting.
#' @param df A data frame containing strategy labels, scenario identifiers and target variable.
#' @param scenario Character string specifying the scenario to filter.
#' @param S_var Character string of the target variable name for efficiency calculation.
#' @param conf_level Numeric value of confidence level for CI calculation (default = 0.95).
#' @return A tibble with comprehensive statistics of BWR_S.
#' @export
compute_BWR <- function(df, scenario, S_var, conf_level = 0.95) {
  # -------------------------- Input Validation -------------------------- #
  # Stop execution if input df is not a data frame
  if (!is.data.frame(df)) stop("Argument 'df' must be a valid data frame")
  # Stop execution if target scenario is not present in df
  if (!scenario %in% df$scenario) stop(paste0("Scenario '", scenario, "' not found in df$scenario"))
  # Stop execution if target variable is not present in df
  if (!S_var %in% colnames(df)) stop(paste0("Variable '", S_var, "' not found in data frame columns"))
  
  # -------------------------- Data Preprocessing -------------------------- #
  # Filter data for the target scenario and remove rows with missing values in S_var
  df_filtered <- df %>% 
    filter(scenario == !!scenario) %>% 
    drop_na(all_of(S_var))  # Exclude NA values to avoid calculation errors
  
  # Extract cost values (reuse external get_cost function)
  X <- get_cost(df_filtered, scenario)
  
  # Subset data for best and worst strategy groups respectively
  best_grp <- df_filtered %>% filter(strategy == "best")
  worst_grp <- df_filtered %>% filter(strategy == "worst")
  
  # -------------------------- Efficiency Calculation -------------------------- #
  # Calculate efficiency metric E = S / cost for best and worst groups
  E_best <- if (nrow(best_grp) > 0) {
    best_grp[[S_var]] / get_cost(best_grp, scenario)
  } else {
    numeric(0)  # Return empty vector if best group has no data
  }
  
  E_worst <- if (nrow(worst_grp) > 0) {
    worst_grp[[S_var]] / get_cost(worst_grp, scenario)
  } else {
    numeric(0)  # Return empty vector if worst group has no data
  }
  
  # Filter valid efficiency values (finite & positive) to exclude outliers/errors
  E_best <- E_best[is.finite(E_best) & E_best > 0]
  E_worst <- E_worst[is.finite(E_worst) & E_worst > 0]
  
  # -------------------------- Handle Empty Groups -------------------------- #
  # Return tibble with NA values if either group has no valid data
  if (length(E_best) == 0 | length(E_worst) == 0) {
    return(tibble(
      scenario = scenario,
      S_var = S_var,
      n_best = length(E_best),
      n_worst = length(E_worst),
      mean_best = NA_real_,
      mean_worst = NA_real_,
      mean_ratio = NA_real_,
      se_best = NA_real_,
      se_worst = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      cv_best = NA_real_,
      cv_worst = NA_real_
    ))
  }
  
  # -------------------------- Core Statistical Metrics -------------------------- #
  # Calculate mean, standard error (SE) and coefficient of variation (CV, %)
  mean_best <- mean(E_best)
  mean_worst <- mean(E_worst)
  se_best <- sd(E_best) / sqrt(length(E_best))  # SE = SD / sqrt(n)
  se_worst <- sd(E_worst) / sqrt(length(E_worst))
  cv_best <- (sd(E_best) / mean_best) * 100     # CV as percentage
  cv_worst <- (sd(E_worst) / mean_worst) * 100
  
  # Calculate the core metric: ratio of mean efficiency (best / worst)
  mean_ratio <- mean_best / mean_worst
  
  # -------------------------- Confidence Interval Calculation -------------------------- #
  # Calculate SE of log-transformed ratio to avoid bias (log transformation method)
  se_log_ratio <- sqrt((se_best / mean_best)^2 + (se_worst / mean_worst)^2)
  # Degrees of freedom for t-distribution (independent samples approximation)
  df_total <- length(E_best) + length(E_worst) - 2
  # Critical t-value for specified confidence level
  t_crit <- qt((1 + conf_level) / 2, df = df_total)
  
  # Calculate CI in log scale then back-transform to original scale
  log_ratio <- log(mean_ratio)
  log_ci_lower <- log_ratio - t_crit * se_log_ratio
  log_ci_upper <- log_ratio + t_crit * se_log_ratio
  ci_lower <- exp(log_ci_lower)  # Back transformation
  ci_upper <- exp(log_ci_upper)
  
  # -------------------------- Result Compilation -------------------------- #
  # Assemble all statistics into a structured tibble for output
  result <- tibble(
    scenario = scenario,          # Scenario identifier
    S_var = S_var,                # Target variable name
    n_best = length(E_best),      # Sample size of best group
    n_worst = length(E_worst),    # Sample size of worst group
    mean_best = mean_best,        # Mean efficiency of best group
    mean_worst = mean_worst,      # Mean efficiency of worst group
    mean_ratio = mean_ratio,      # Core ratio: mean_best / mean_worst
    se_best = se_best,            # Standard error of best group mean
    se_worst = se_worst,          # Standard error of worst group mean
    ci_lower = ci_lower,          # Lower bound of ratio 95% CI
    ci_upper = ci_upper,          # Upper bound of ratio 95% CI
    cv_best = cv_best,            # Coefficient of variation of best group (%)
    cv_worst = cv_worst           # Coefficient of variation of worst group (%)
  )
  
  return(result)
}


scenarios <- unique(dt$scenario)
shortfalls <- c("LS", "WS", "DS")

BWR_results <- expand_grid(
  scenario = scenarios,
  S_var = shortfalls
) %>%
  mutate(
    stats = map2(scenario, S_var, ~ compute_BWR(dt, .x, .y))
  ) %>%
  select(-scenario, -S_var) %>%
  unnest(cols = stats) %>%
  rename(
    bwr_mean_ratio = mean_ratio,    
    ci_95_lower = ci_lower,         
    ci_95_upper = ci_upper,         
    cv_best_pct = cv_best,          
    cv_worst_pct = cv_worst
  )
BWR_results

BWR_results$S_var <- factor(BWR_results$S_var,levels = c("LS","WS","DS"))
pp1_r <- ggplot(
  BWR_results %>% filter(scenario == "scenario_1"),
  aes(x = S_var, y = bwr_mean_ratio,ymin = ci_95_lower, ymax = ci_95_upper)) +
  geom_col(fill = "#3182bd", alpha = 0.5, width = 0.6) +
  geom_errorbar(width = 0.2, linewidth = 0.2,color = "grey") +
  # geom_text(
  #   aes(label = round(bwr_mean_ratio, 2)),size = 3, vjust = -0.5,  color = "black") +
  labs(x = "Shortfall", y = "BWR (best-to-worst ratio)") +
  theme_classic() +
  scale_y_continuous(limits = c(0, 10),  expand = c(0, 0)) +
  theme(
    axis.text = element_text(size = 5, color = "black"),
    axis.title = element_text(size = 6, color = "black"),
    axis.line = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.ticks.length = unit(0.05, "cm"),
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  )



pp2_r <- ggplot(
  BWR_results %>% filter(scenario == "scenario_2"),
  aes(x = S_var, y = bwr_mean_ratio,ymin = ci_95_lower, ymax = ci_95_upper)) +
  geom_col(fill = "#3182bd", alpha = 0.5, width = 0.6) +
  geom_errorbar(width = 0.2, linewidth = 0.2,color = "grey") +
  # geom_text(
  #   aes(label = round(bwr_mean_ratio, 2)),size = 3, vjust = -0.5,  color = "black") +
  labs(x = "Shortfall", y = "BWR (best-to-worst ratio)") +
  theme_classic() +
  scale_y_continuous(limits = c(0, 5),  expand = c(0, 0)) +
  theme(
    axis.text = element_text(size = 5, color = "black"),
    axis.title = element_text(size = 6, color = "black"),
    axis.line = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.ticks.length = unit(0.05, "cm"),
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  )


pp3_r <- ggplot(
  BWR_results %>% filter(scenario == "scenario_3"),
  aes(x = S_var, y = bwr_mean_ratio,ymin = ci_95_lower, ymax = ci_95_upper)) +
  geom_col(fill = "#3182bd", alpha = 0.5, width = 0.6) +
  geom_errorbar(width = 0.2, linewidth = 0.2,color = "grey") +
  # geom_text(
  #   aes(label = round(bwr_mean_ratio, 2)),size = 3, vjust = -0.5,  color = "black") +
  labs(x = "Shortfall", y = "BWR (best-to-worst ratio)") +
  theme_classic() +
  scale_y_continuous(limits = c(0, 5),  expand = c(0, 0)) +
  theme(
    axis.text = element_text(size = 5, color = "black"),
    axis.title = element_text(size = 6, color = "black"),
    axis.line = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.ticks.length = unit(0.05, "cm"),
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  )


pp4_r <- ggplot(
  BWR_results %>% filter(scenario == "scenario_4"),
  aes(x = S_var, y = bwr_mean_ratio,ymin = ci_95_lower, ymax = ci_95_upper)) +
  geom_col(fill = "#3182bd", alpha = 0.5, width = 0.6) +
  geom_errorbar(width = 0.2, linewidth = 0.2,color = "grey") +
  # geom_text(
  #   aes(label = round(bwr_mean_ratio, 2)),size = 3, vjust = -0.5,  color = "black") +
  labs(x = "Shortfall", y = "BWR (best-to-worst ratio)") +
  theme_classic() +
  scale_y_continuous(limits = c(0, 5),  expand = c(0, 0)) +
  theme(
    axis.text = element_text(size = 5, color = "black"),
    axis.title = element_text(size = 6, color = "black"),
    axis.line = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.ticks.length = unit(0.05, "cm"),
    plot.background = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  )


################################################################################
################################################################################
get_denominator_vector <- function(df_subset, scenario) {
  valid_scenarios <- c("scenario_1", "scenario_2", "scenario_3", "scenario_4")
  if (!scenario %in% valid_scenarios) {
    stop(paste("Unsupported scenario:", scenario, "| Valid scenarios:", paste(valid_scenarios, collapse = ", ")))
  }
  if (scenario == "scenario_1") {
    return(df_subset$Investment)
  } else if (scenario == "scenario_2") {
    return(df_subset$Socioeconomics)
  } else if (scenario == "scenario_3") {
    return(df_subset$Governance)
  } else if (scenario == "scenario_4") {
    return(df_subset$Investment * df_subset$Socioeconomics * df_subset$Governance)
  }
}
calculate_BWR_with_CI <- function(df, scenario) {
  best  <- df %>% filter(strategy == "best")
  worst <- df %>% filter(strategy == "worst")
  if (nrow(best) == 0 | nrow(worst) == 0) {
    result <- tibble(
      scenario = scenario,         
      n_best = 0, 
      n_worst = 0,
      BWR = NA_real_,               
      BWR_CI_95_lower = NA_real_,
      BWR_CI_95_upper = NA_real_
    )
    cat(sprintf("Scenario %s: No valid best/worst group (n_best=%d, n_worst=%d)\n", 
                scenario, result$n_best, result$n_worst))
    return(result)
  }
  
  eps <- 1e-6
  denom_best <- get_denominator_vector(best, scenario)
  E_best_raw <- best$Shortfall / denom_best
  E_best     <- E_best_raw[is.finite(E_best_raw) & !is.na(E_best_raw) & denom_best > eps]
  denom_worst <- get_denominator_vector(worst, scenario)
  E_worst_raw <- worst$Shortfall / denom_worst
  E_worst     <- E_worst_raw[is.finite(E_worst_raw) & !is.na(E_worst_raw) & denom_worst > eps]
  if (length(E_best) == 0 | length(E_worst) == 0) {
    result <- tibble(
      scenario = scenario,
      n_best = length(E_best),
      n_worst = length(E_worst),
      BWR = NA_real_,
      BWR_CI_95_lower = NA_real_,
      BWR_CI_95_upper = NA_real_
    )
    cat(sprintf("Scenario %s: No valid E_best/E_worst (n_best=%d, n_worst=%d)\n", 
                scenario, result$n_best, result$n_worst))
    return(result)
  }
  mean_best  <- mean(E_best, na.rm = TRUE)
  mean_worst <- mean(E_worst, na.rm = TRUE)
  BWR        <- mean_best / mean_worst
  
  se_best  <- sd(E_best, na.rm = TRUE) / sqrt(length(E_best))
  se_worst <- sd(E_worst, na.rm = TRUE) / sqrt(length(E_worst))
  if (mean_best == 0 | mean_worst == 0) {
    result <- tibble(
      scenario = scenario,
      n_best = length(E_best),
      n_worst = length(E_worst),
      BWR = 0,
      BWR_CI_95_lower = 0,
      BWR_CI_95_upper = 0
    )
    cat(sprintf("Scenario %s: BWR=0 (mean_best=%.4f, mean_worst=%.4f)\n", 
                scenario, mean_best, mean_worst))
    return(result)
  }
  se_log_BWR <- sqrt((se_best / mean_best)^2 + (se_worst / mean_worst)^2)
  df_total   <- length(E_best) + length(E_worst) - 2
  t_crit     <- if (df_total > 0) qt(0.975, df = df_total) else 1.96 
  
  log_BWR     <- log(BWR)
  CI_lower    <- exp(log_BWR - t_crit * se_log_BWR)
  CI_upper    <- exp(log_BWR + t_crit * se_log_BWR)
  result <- tibble(
    scenario = scenario,
    n_best = length(E_best),          
    n_worst = length(E_worst),
    BWR = round(BWR, 3),              
    BWR_CI_95_lower = round(CI_lower, 4),  
    BWR_CI_95_upper = round(CI_upper, 3)   
  )
  cat(sprintf(
    "BWR = %.3f (95%% CI: %.3f–%.3f), n_best = %d, n_worst = %d\n",
    result$BWR, result$BWR_CI_95_lower, result$BWR_CI_95_upper,
    result$n_best, result$n_worst
  ))
  #return(result)
}

# ---------------------- scenario_1 ---------------------- #
calculate_BWR_with_CI(df_plot_realm_1 %>% mutate(scenario = "scenario_1") , scenario = "scenario_1")
# BWR = 2.242 (95% CI: 1.425–3.529), n_best = 13, n_worst = 17

# ---------------------- scenario_2 ---------------------- #
calculate_BWR_with_CI(df_plot_realm_2 %>% mutate(scenario = "scenario_2") , scenario = "scenario_2")
# BWR = 0.786 (95% CI: 0.250–2.468), n_best = 12, n_worst = 12

# ---------------------- scenario_3 ---------------------- #
calculate_BWR_with_CI(df_plot_realm_3 %>% mutate(scenario = "scenario_3") , scenario = "scenario_3")
# BWR = 0.089 (95% CI: 0.017–0.455), n_best = 8, n_worst = 9

# ---------------------- scenario_4 ---------------------- #
calculate_BWR_with_CI(df_plot_realm_4 %>% mutate(scenario = "scenario_4") , scenario = "scenario_4")
# BWR = 0.169 (95% CI: 0.038–0.759), n_best = 23, n_worst = 13
################################################################################
library(sf)
library(ggplot2)
sf_use_s2(FALSE)
inland <- readRDS("input/raw/country_20251212.rds")
moll_proj <- st_crs("+proj=moll")
lat_points <- data.frame(lon = 0, lat = c(-60, 90)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)
projected_points <- st_transform(lat_points, crs = moll_proj)
y_limits <- st_coordinates(projected_points)[, "Y"]


case1 <- inland %>% left_join(df_plot_realm_1[,c(2,8)], by = "iso3")
pp1 <- ggplot()+
  ggrastr::rasterise(geom_sf(data =case1 ,aes(fill = front_rank),colour = "white", linewidth = 0.03),dpi = 300)+
  scale_fill_gradientn(name = "Pareto frontier",
                       colours = c("#FFC000", "grey90","#00B050"),
                       limits = c(0, 1),                 
                       breaks = c(0, 1), 
                       na.value = "grey70",
                       labels = c("Worst","Best"),
                       guide = guide_colbar(
                         show = c(TRUE, TRUE), 
                         oob = "squish")
  )+
  theme_void()+
  coord_sf(crs = "+proj=moll +lon_0=0", ylim = y_limits, expand = FALSE)+
  theme(legend.position = "none",
        plot.margin = margin(0,0,0,0))

pp1 <- ggdraw() +
  draw_plot(pp1, 0, 0, 1, 1) +                             
  draw_plot(pp1_r, 0.06, 0.02, 0.2, 0.8) +
  draw_label(label = expression(Delta[NSCI] == 0.02),
             x = 0.58, y = 0.08,hjust = 1, vjust = 1, size = 7.2)+
  draw_label(label = expression("E-test, " * italic(p) * " < 0.05"),
             x = 0.8, y = 0.08,hjust = 1, vjust = 1, size = 7.2)+
  draw_label(label = expression(BWR[NSCI] * ' = 2.42'),
             x = 0.3, y = 0.98,hjust = 1, vjust = 1, size = 7.5)

################################################################################
case2 <- inland %>% left_join(df_plot_realm_2 [,c(2,8)], by = "iso3")
pp2 <- ggplot()+
  ggrastr::rasterise(geom_sf(data =case2 ,aes(fill = front_rank),colour = "white", linewidth = 0.03),dpi = 300)+
  scale_fill_gradientn(name = "Pareto frontier",
                       colours = c("#FFC000", "grey90","#00B050"),
                       limits = c(0, 1),                 
                       breaks = c(0, 1), 
                       na.value = "grey70",
                       labels = c("Worst","Best"),
                       guide = guide_colbar(
                         show = c(TRUE, TRUE), 
                         oob = "squish")
  )+
  theme_void()+
  coord_sf(crs = "+proj=moll +lon_0=0", ylim = y_limits, expand = FALSE)+
  theme(legend.position = "none",
        plot.margin = margin(0,0,0,0))

pp2 <- ggdraw() +
  draw_plot(pp2, 0, 0, 1, 1) +                             
  draw_plot(pp2_r, 0.06, 0.02, 0.2, 0.8)  +
  draw_label(label = expression(Delta[NSCI] == 0.09),
             x = 0.58, y = 0.08,hjust = 1, vjust = 1, size = 7.2)+
  draw_label(label = expression("E-test, " * italic(p) * " = 0.08"),
             x = 0.8, y = 0.08,hjust = 1, vjust = 1, size = 7.2)+
  draw_label(label = expression(BWR[NSCI] * ' = 0.79'),
             x = 0.3, y = 0.98,hjust = 1, vjust = 1, size = 7.5)


################################################################################
case3 <- inland %>% left_join(df_plot_realm_3[,c(2,8)], by = "iso3")
pp3 <- ggplot()+
  ggrastr::rasterise(geom_sf(data =case3 ,aes(fill = front_rank),colour = "white", linewidth = 0.03),dpi = 300)+
  scale_fill_gradientn(name = "Pareto frontier",
                       colours = c("#FFC000", "grey90","#00B050"),
                       limits = c(0, 1),                 
                       breaks = c(0, 1), 
                       na.value = "grey70",
                       labels = c("Worst","Best"),
                       guide = guide_colbar(
                         show = c(TRUE, TRUE), 
                         oob = "squish")
  )+
  theme_void()+
  coord_sf(crs = "+proj=moll +lon_0=0", ylim = y_limits, expand = FALSE)+
  theme(legend.position = "none",
        plot.margin = margin(0,0,0,0))

pp3 <- ggdraw() +
  draw_plot(pp3, 0, 0, 1, 1) +                             
  draw_plot(pp3_r, 0.06, 0.02, 0.2, 0.8)  +
  draw_label(label = expression(Delta[NSCI] == 0.11),
             x = 0.58, y = 0.08,hjust = 1, vjust = 1, size = 7.2)+
  draw_label(label = expression("E-test, " * italic(p) * " = 0.09"),
             x = 0.8, y = 0.08,hjust = 1, vjust = 1, size = 7.2)+
  draw_label(label = expression(BWR[NSCI] * ' = 0.09'),
             x = 0.3, y = 0.98,hjust = 1, vjust = 1, size = 7.5)
################################################################################
case4 <- inland %>% left_join(df_plot_realm_4[,c(2,12)], by = "iso3")
pp4 <- ggplot()+
  ggrastr::rasterise(geom_sf(data =case4 ,aes(fill = front_rank),colour = "white", linewidth = 0.03),dpi = 300)+
  scale_fill_gradientn(name = "Pareto frontier",
                       colours = c("#FFC000", "grey90","#00B050"),
                       limits = c(0, 1),                 
                       breaks = c(0, 1), 
                       na.value = "grey70",
                       labels = c("Worst","Best"),
                       guide = guide_colbar(
                         show = c(TRUE, TRUE), 
                         oob = "squish")
  )+
  theme_void()+
  coord_sf(crs = "+proj=moll +lon_0=0", ylim = y_limits, expand = FALSE)+
  theme(legend.position = "none",
        plot.margin = margin(0,0,0,0))
pp4 <- ggdraw() +
  draw_plot(pp4, 0, 0, 1, 1) +                             
  draw_plot(pp4_r, 0.06, 0.02, 0.2, 0.8)  +
  draw_label(label = expression(Delta[NSCI] == -0.22),
             x = 0.58, y = 0.08,hjust = 1, vjust = 1, size = 7.2)+
  draw_label(label = expression("E-test, " * italic(p) * " < 0.05"),
             x = 0.8, y = 0.08,hjust = 1, vjust = 1, size = 7.2)+
  draw_label(label = expression(BWR[NSCI] * ' = 0.17'),
             x = 0.3, y = 0.98,hjust = 1, vjust = 1, size = 7.5)
################################################################################
row_pair <- function(left_plot, right_plot,labels) {
  plot_grid(
    left_plot  + theme(plot.margin = margin(5, 5, 0, 0)),
    right_plot + theme(plot.margin = margin(5, 0, 0, -5)),
    ncol = 2,
    rel_widths = c(0.4, 0.6),
    label_size = 9,
    label_fontface = "bold",
    align = "h",
    axis  = "tb",
    labels = labels
  )
}

row1 <- row_pair(p1, pp1, labels = c("A","B"))
row2 <- row_pair(p2, pp2, labels = c("C","D"))
row3 <- row_pair(p3, pp3, labels = c("E","F"))
row4 <- row_pair(p4, pp4, labels = c("G","H"))

row1/row2/row3/row4

ggsave("figures/supplement/Figure_S12.png",dpi = 300, units = "cm", width = 18, height = 20)




