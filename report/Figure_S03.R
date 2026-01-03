# ------------------------------------------------------------------------------
# Supplementary Figure S3
# Standardized effects of predictors on species documentation processes in freshwater fishes across continents

rm(list = ls())

library(brms)
library(dplyr)
library(purrr)
library(tidybayes)
library(ggplot2)
library(ggh4x)
library(ggtext)
library(legendry) # This replaces ggh4x for nested guides
library(sf)
library(ggplot2)
library(rnaturalearth) 
library(ggspatial) 
library(ggtext) # 
library(patchwork)
options(sf_use_s2 = F)
options(warn = -1)

compute_continent_slopes <- function(v, draws, cont_levels, baseline_continent = "Africa") {
  # Name of the main-effect coefficient for predictor v
  base_name <- paste0("b_", v)
  
  if (!base_name %in% names(draws)) {
    stop("No coefficient found for ", base_name)
  }
  
  # Baseline slope β_v (effect in the reference continent)
  base <- draws[[base_name]]
  
  # Keep .draw (and optionally .chain, .iteration) to preserve posterior structure
  meta_cols <- draws %>% dplyr::select(.chain, .iteration, .draw)
  
  # For each continent, construct the posterior slope
  out <- map_dfr(cont_levels, function(ct) {
    if (ct == baseline_continent) {
      # Baseline continent: slope = β_v
      slope <- base
    } else {
      # Non-baseline continent: slope = β_v + β_v:continent=ct
      inter_name <- paste0("b_", v, ":continent", ct)
      if (inter_name %in% names(draws)) {
        slope <- base + draws[[inter_name]]
      } else {
        # Safety fallback: no interaction term found (should not occur given your formula)
        slope <- base
      }
    }
    
    tibble(
      continent = ct,
      slope     = slope
    )
  })
  
  # Bind back posterior meta info and add variable name
  bind_cols(
    meta_cols[rep(1:nrow(meta_cols), times = length(cont_levels)), ],
    out
  ) %>%
    mutate(variable = v)
}


# The global marginal slope for predictor $x$ is computed as
# $$
#   \beta_x^{\text{global}} 
# = \beta_x 
# + \sum_{c \ne \text{Africa}} w_c \,\beta_{x:\text{continent}=c},
# $$
#   where $w_c$ denotes the proportion of observations from continent $c$.

compute_global_slope <- function(v, draws, w, cont_levels, baseline_continent = "Africa") {
  # Construct the main-effect coefficient name for variable v
  base_name <- paste0("b_", v)
  
  # Ensure that the baseline coefficient exists in the posterior draws
  if (!base_name %in% names(draws)) {
    stop("No coefficient found for ", base_name)
  }
  
  # Extract the baseline slope β_v (i.e., effect in the reference continent)
  base <- draws[[base_name]]
  
  # Initialize the global slope with the baseline effect
  slope <- base
  
  # Add weighted interaction terms for all non-baseline continent
  for (ct in cont_levels[cont_levels != baseline_continent]) {
    inter_name <- paste0("b_", v, ":continent", ct)
    if (inter_name %in% names(draws)) {
      # Weighted contribution from continent-specific deviation
      slope <- slope + as.numeric(w[ct]) * draws[[inter_name]]
    } else {
      # Fallback: if no interaction term exists (should not happen in this model)
      slope <- slope + 0
    }
  }
  
  # Return a tibble containing variable name and posterior samples of the global slope
  tibble(
    variable = v,
    slope    = slope
  )
}

################################################################################
fit <- readRDS("output/model/country_linnaean_all.rds")$weibull 
# bayes_R2(fit)
# Estimate    Est.Error      Q2.5    Q97.5
# 0.8476622 0.00127659 0.8448291 0.8497559

draws <- as_draws_df(fit)
dat <- fit$data
levels(dat$continent)
w <- prop.table(table(dat$continent))
cont_levels <- names(w)
n_continent <- as.data.frame(table(dat$continent))
names(n_continent) <- c("continent","n")

vars_linnaean <- c(
  "body_size", "taxonmic_effort", "taxonomic_activity", "country_area",
  "range_size", "elevation", "latitude", "discharge", "watertemp",
  "preserved_specimen", "sampling_effort","sequencing_effort", "range_rarity", "population_density"
)

global_slopes_linnaean <- map_dfr(
  vars_linnaean,
  ~ compute_global_slope(.x, draws = draws, w = w, cont_levels = cont_levels)
) %>%
  mutate(slope = as.numeric(slope)) %>%   
  group_by(variable) %>%
  tidybayes::mean_qi(slope) %>%
  ungroup() %>%
  mutate(continent = "Global",
         n = nrow(dat))

continent_slopes_linnaean <- map_dfr(
  vars_linnaean,
  ~ compute_continent_slopes(
    v             = .x,
    draws         = draws,
    cont_levels   = cont_levels,
    baseline_continent = "Africa"
  )
) %>%
  mutate(
    slope    = as.numeric(slope),
    continent = factor(continent, levels = cont_levels),
    variable  = factor(variable, levels = vars_linnaean)
  ) %>%
  group_by(variable, continent) %>%
  tidybayes::mean_qi(slope) %>%
  left_join(n_continent,by = "continent")


linnaean_slopes <- bind_rows(
  continent_slopes_linnaean,   
  global_slopes_linnaean       
) %>%
  mutate(shortfall = "Linnaean")

rm(fit,w,vars_linnaean,cont_levels,continent_slopes_linnaean,global_slopes_linnaean,draws,dat,n_continent)
################################################################################
fit <- readRDS("output/model/country_wallacean_all.rds")$lognormal
# bayes_R2(fit)
# Estimate    Est.Error      Q2.5    Q97.5
# 0.5796673 0.003815058 0.571929 0.5865687
draws <- as_draws_df(fit)
dat <- fit$data
levels(dat$continent)
w <- prop.table(table(dat$continent))
cont_levels <- names(w)
n_continent <- as.data.frame(table(dat$continent))
names(n_continent) <- c("continent","n")

vars_wallacean <- c(
  "body_size", "taxonmic_effort", "taxonomic_activity", "country_area",
  "range_size", "elevation", "latitude", "discharge", "watertemp",
  "preserved_specimen", "sequencing_effort", "range_rarity", "population_density"
)

global_slopes_wallacean <- map_dfr(
  vars_wallacean,
  ~ compute_global_slope(.x, draws = draws, w = w, cont_levels = cont_levels)
) %>%
  mutate(slope = as.numeric(slope)) %>%   
  group_by(variable) %>%
  tidybayes::mean_qi(slope) %>%
  ungroup() %>%
  mutate(continent = "Global",
         n = nrow(dat))

continent_slopes_wallacean <- map_dfr(
  vars_wallacean,
  ~ compute_continent_slopes(
    v             = .x,
    draws         = draws,
    cont_levels   = cont_levels,
    baseline_continent = "Africa"
  )
) %>%
  mutate(
    slope    = as.numeric(slope),
    continent = factor(continent, levels = cont_levels),
    variable  = factor(variable, levels = vars_wallacean)
  ) %>%
  group_by(variable, continent) %>%
  tidybayes::mean_qi(slope) %>%
  left_join(n_continent,by = "continent")


wallacean_slopes <- bind_rows(
  continent_slopes_wallacean,   
  global_slopes_wallacean       
) %>%
  mutate(shortfall = "Wallacean")

rm(fit,w,vars_wallacean,cont_levels,continent_slopes_wallacean,global_slopes_wallacean,draws,dat,n_continent)
################################################################################
fit <- readRDS("output/model/country_darwinian_all.rds")$weibull
# bayes_R2(fit)
# Estimate    Est.Error      Q2.5    Q97.5
# R2 0.5378377 0.004422993 0.5285995 0.5460068
draws <- as_draws_df(fit)
dat <- fit$data
levels(dat$continent)
w <- prop.table(table(dat$continent))
cont_levels <- names(w)
n_continent <- as.data.frame(table(dat$continent))
names(n_continent) <- c("continent","n")

vars_darwinian <- c(
  "body_size", "taxonmic_effort", "taxonomic_activity", "country_area",
  "range_size", "elevation", "latitude", "discharge", "watertemp",
  "preserved_specimen", "sampling_effort","range_rarity", "population_density"
)

global_slopes_darwinian <- map_dfr(
  vars_darwinian,
  ~ compute_global_slope(.x, draws = draws, w = w, cont_levels = cont_levels)
) %>%
  mutate(slope = as.numeric(slope)) %>%   
  group_by(variable) %>%
  tidybayes::mean_qi(slope) %>%
  ungroup() %>%
  mutate(continent = "Global",
         n = nrow(dat))

continent_slopes_darwinian <- map_dfr(
  vars_darwinian,
  ~ compute_continent_slopes(
    v             = .x,
    draws         = draws,
    cont_levels   = cont_levels,
    baseline_continent = "Africa"
  )
) %>%
  mutate(
    slope    = as.numeric(slope),
    continent = factor(continent, levels = cont_levels),
    variable  = factor(variable, levels = vars_darwinian)
  ) %>%
  group_by(variable, continent) %>%
  tidybayes::mean_qi(slope) %>%
  left_join(n_continent,by = "continent")


darwinian_slopes <- bind_rows(
  continent_slopes_darwinian,   
  global_slopes_darwinian       
) %>%
  mutate(shortfall = "Darwinian")

rm(fit,w,vars_darwinian,cont_levels,continent_slopes_darwinian,global_slopes_darwinian,draws,dat,n_continent)
################################################################################
data <- bind_rows(linnaean_slopes, wallacean_slopes, darwinian_slopes) %>%
  mutate(
    shortfall = recode(
      shortfall,
      "Darwinian" = "Sequencing",
      "Wallacean" = "Geolocation",
      "Linnaean"  = "Description"
    )
  )

rename_map <- c(
  body_size          = "Body size",
  preserved_specimen = "N. preserved\nspecimens",
  range_size         = "Range size",
  range_rarity       = "Range rarity",
  latitude           = "Latitude",
  country_area       = "Country area",
  discharge          = "Streamflow",
  watertemp          = "Water temperature",
  elevation          = "Elevation",
  population_density = "Human density",
  taxonmic_effort    = "Taxonomic effort",
  taxonomic_activity = "Taxonomic activity",
  sampling_effort    = "Sampling effort",
  sequencing_effort  = "Sequencing effort"
)

var_levels <- rev(unname(rename_map))

data <- data %>%
  mutate(
    variable  = recode(variable, !!!rename_map),
    variable  = factor(variable, levels = var_levels),
    shortfall = factor(shortfall,
                       levels = c("Sequencing", "Geolocation", "Description")),
    continent = factor(
      continent,
      levels = c(
        "North America",
        "South America",
        "Europe",
        "Africa",
        "Asia",
        "Oceania","Global"
      )
    ),
    group = case_when(
      variable %in% c("Sequencing effort","Sampling effort",
                      "Taxonomic activity","Taxonomic effort") ~ "Activity",
      variable %in% c("Human density","Elevation") ~ "Access",
      variable %in% c("Water temperature","Streamflow") ~ "Environment",
      variable %in% c("Country area","Latitude",
                      "Range rarity","Range size") ~ "Geography",
      variable %in% c("N. preserved\nspecimens","Body size") ~ "Biology",
      TRUE ~ NA_character_
    ),
    group = factor(group,
                   levels = c("Activity","Access",
                              "Environment","Geography","Biology"))
  )

set2_colors <- c(
  "Description" = "#D55E00",
  "Geolocation" = "#009E73",
  "Sequencing"  = "#435792"
)

continent <- c(
  Africa =
    "Africa<br>
        <span style='color:#D55E00;'>n=8,037</span><br>
        <span style='color:#009E73;'>n=10,681</span><br>
        <span style='color:#435792;'>n=10,534</span>",
  Oceania =
    "Oceania<br>
        <span style='color:#D55E00;'>n=597</span><br>
        <span style='color:#009E73;'>n=894</span><br>
        <span style='color:#435792;'>n=844</span>",
  Asia =
    "Asia<br>
        <span style='color:#D55E00;'>n=5,780</span><br>
        <span style='color:#009E73;'>n=10,271</span><br>
        <span style='color:#435792;'>n=9,321</span>",
  `North America` =
    "North America<br>
        <span style='color:#D55E00;'>n=2,395</span><br>
        <span style='color:#009E73;'>n=2,972</span><br>
        <span style='color:#435792;'>n=2,887</span>",
  `South America` =
    "South America<br>
        <span style='color:#D55E00;'>n=10,858</span><br>
        <span style='color:#009E73;'>n=13,871</span><br>
        <span style='color:#435792;'>n=13,633</span>",
  Europe =
    "Europe<br>
        <span style='color:#D55E00;'>n=1,131</span><br>
        <span style='color:#009E73;'>n=2,746</span><br>
        <span style='color:#435792;'>n=2,659</span>",
  Global =
    "Global<br>
        <span style='color:#D55E00;'>bayes R<sup>2</sup>=0.85</span><br>
        <span style='color:#009E73;'>bayes R<sup>2</sup>=0.58</span><br>
        <span style='color:#435792;'>bayes R<sup>2</sup>=0.54</span>"
)


p <- ggplot(data = data, aes(x = slope, y =  interaction(variable, group, sep = ":"),fill = shortfall,colour = shortfall)) +
  annotate("rect", ymin = 0.5, ymax = 4.5,
           xmin = -Inf, xmax = Inf, fill = "grey95") +
  annotate("rect", ymin = 4.5, ymax = 6.5,
           xmin = -Inf, xmax = Inf, fill = "grey90") +
  annotate("rect", ymin = 6.5, ymax = 8.5,
           xmin = -Inf, xmax = Inf, fill = "grey95") +
  annotate("rect", ymin = 8.5, ymax = 12.5,
           xmin = -Inf, xmax = Inf, fill = "grey90") +
  annotate("rect", ymin = 12.5, ymax = 14.5,
           xmin = -Inf, xmax = Inf, fill = "grey95") +
  geom_vline(xintercept = 0, linetype = 2, size = 0.4, color = "gray50") +  # Add a vertical line at zero for reference
  geom_pointrange(aes(xmin = .lower, xmax = .upper),
                  position = position_dodge(0.6), 
                  linewidth = 0.5, shape = 21, stroke = 0.2, size = 0.46) +
  scale_fill_manual(name = "First documentation", values = set2_colors, limits = c("Description", "Geolocation", "Sequencing")) +  # Set custom colors for fill
  scale_color_manual(name = "First documentation", values = set2_colors, limits = c("Description", "Geolocation", "Sequencing")) +  # Set custom colors for color
  theme_classic()+
  facet_wrap2(~ continent, nrow = 1, scales = "free_x",labeller = labeller(continent = continent)) +
  geom_vline(xintercept = 0, linetype = 2)+
  scale_y_discrete(
    expand = c(0, 0),
    guide  = legendry::guide_axis_nested(key = key_range_auto(sep = ":"),
                                         # Customise the appearance of different levels using levels_text
                                         levels_text = list(
                                           element_text(colour = "black",size = 10), # Level 1 text (variable)
                                           element_text(colour = "black", angle = 90, size = 10, hjust = 0.5)  # Level 2 text (group)
                                         )
    )
  ) +
  theme(strip.background = element_blank(),  # Remove background for facet strips
        #panel.background = element_rect(fill = "grey90"),
        legend.title = element_text(colour = "black", size = 9),  # Customize legend title
        legend.text = element_text(colour = "black", size = 8),  # Customize legend text
        legend.key.height = unit(0.5, "cm"),  # Adjust legend key height
        axis.title = element_text(colour = "black", size = 10),  # Customize axis titles
        axis.text.x = element_text(colour = "black", size = 8),  # Customize x-axis text
        strip.text = element_markdown(size = 8),  # Customize strip text
        axis.text.y = element_markdown(size = 8, colour = "black")) + # Use markdown for y-axis text color
  ylab("") +  # Remove y-axis label
  xlab("Standardized coefficient (95% CI)")  # x-axis label for standardized coefficients with 95% CI


################################################################################
inland <- readRDS("input/raw/country.rds")
country <- read.csv("input/raw/country_list.csv")
inland <- inland %>% left_join(country[,c(1,5)], by = "iso3")
inland_grouped <- inland %>%
  group_by(continent) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") 
inland_grouped_fixed <- st_wrap_dateline(inland_grouped, options = c("WRAPDATELINE=YES"))
#plot(inland_grouped_fixed)


world <- ne_countries(scale = "small", returnclass = "sf") %>%
  st_wrap_dateline(options = c("WRAPDATELINE=YES"))


p1 <- ggplot() +
  ggrastr::rasterise(geom_sf(data = world, fill = "gray80", color = NA, size = 0.2),dpi=300) +  
  ggrastr::rasterise(geom_sf(data = filter(inland_grouped_fixed, continent == "North America"), fill = "#F0E442", color = NA, size = 0.5),dpi=300) +  
  coord_sf(crs = "+proj=laea +lon_0=-100 +lat_0=30")+
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#ebebe5", size = 0.01))

p2 <- ggplot() +
  ggrastr::rasterise(geom_sf(data = world, fill = "gray80", color = NA, size = 0.2),dpi=300) +  
  ggrastr::rasterise(geom_sf(data = filter(inland_grouped_fixed, continent == "South America"), fill = "#E69F00", color = NA, size = 0.5),dpi=300) +  
  coord_sf(crs = "+proj=laea +lon_0=-90 +lat_0=0")+
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#ebebe5", size = 0.01))

p3 <- ggplot() +
  #geom_sf(data = graticules, color = "black", size = 0.6) +  
  ggrastr::rasterise(geom_sf(data = world, fill = "gray80", color = NA, size = 0),dpi=300) +  
  ggrastr::rasterise(geom_sf(data = filter(inland_grouped_fixed, continent == "Europe"), fill = "#CC79A7", color = NA, size = 0),dpi=300) +  
  coord_sf(crs = "+proj=laea +lon_0=70 +lat_0=30")+
  theme_void() +
  theme(axis.text = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#ebebe5", size = 0.01))

p4 <- ggplot() +
  ggrastr::rasterise(geom_sf(data = world, fill = "gray80", color = NA, size = 0.2),dpi=300) +  
  ggrastr::rasterise(geom_sf(data = filter(inland_grouped_fixed, continent == "Africa"), fill = "#D55E00", color = NA, size = 0.5),dpi=300) +  
  coord_sf(crs = "+proj=laea +lon_0=15 +lat_0=0")+
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#ebebe5", size = 0.01))

p5 <- ggplot() +
  ggrastr::rasterise(geom_sf(data = world, fill = "gray80", color = NA, size = 0.2),dpi=300) +  
  ggrastr::rasterise(geom_sf(data = filter(inland_grouped_fixed, continent == "Asia"), fill = "#009E73", color = NA, size = 0.5),dpi=300) +  
  coord_sf(crs = "+proj=laea +lon_0=90 +lat_0=10")+
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#ebebe5", size = 0.01))

p6 <- ggplot() +
  ggrastr::rasterise(geom_sf(data = world, fill = "gray80", color = NA, size = 0.2),dpi=300) +  
  ggrastr::rasterise(geom_sf(data = filter(inland_grouped_fixed, continent == "Oceania"), fill = "#0072B2", color = NA, size = 0.5),dpi=300) +  
  coord_sf(crs = "+proj=laea +lon_0=140 +lat_0=0")+
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#ebebe5", size = 0.01))

#https://stackoverflow.com/questions/43207947/whole-earth-polygon-for-world-map-in-ggplot2-and-sf
p8 <- ggplot() +
  ggrastr::rasterise(geom_sf(data = world, fill = "gray80", color = NA, size = 0.2),dpi=300) +  
  ggrastr::rasterise(geom_sf(data = inland_grouped_fixed,aes(fill = continent), color = NA, size = 0.5, show.legend = F),dpi = 300) +  
  scale_fill_manual(values = c(
    "Africa" = "#D55E00",
    "Oceania" = "#0072B2",
    "Asia" = "#009E73",
    "North America" = "#F0E442",
    "South America" = "#E69F00",
    "Europe" = "#CC79A7"
  ))+
  coord_sf(crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +datum=WGS84 +units=m +no_defs")+
  theme_minimal() +
  theme(axis.text = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#ebebe5", size = 0.01),
        )


pp <- p1 +p2 +p3 +p4 +p5+p6+p8 +plot_layout(nrow = 1)


# 1) enlarge the top plot margin of p
ppp <- p + 
  theme(
    plot.margin = unit(c(2, 0.05, 0.05,0), "cm"),   # ↑ leave room for inset
    legend.position = c(-0.08, 1.1)           # ← adjust legend in p
  )

# 2) insert pp into p2, placed in the top margin area
p_inset <- ppp +
  inset_element(
    pp,
    left   = 0.19,   # adjust horizontally
    right  = 1,
    bottom = 0.98,   # allows pp to be above main panel
    top    = 1.17,   # >1 allowed because align_to="plot"
    align_to = "plot"
  ) +
  theme(
    legend.position = c(-0.1, 1.1),  # legend for entire inset layout
    legend.justification = c(1, 1),
    plot.margin = unit(c(2, 0.05, 0.05,0), "cm")
  )


ggsave("figures/supplement/Figure_S3.png",dpi = 600, units = "cm", width = 22, height = 18)

