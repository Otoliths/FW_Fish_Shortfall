library(dplyr)
library(arrow)
library(gtable)
seq_annotation <- read_parquet("input/raw/seq_annotation_v1.parquet")
head(seq_annotation)

seq_summary <- seq_annotation %>%
  filter(!is.na(gene), length_bp > 0) %>%      
  group_by(gene) %>%
  summarise(
    n_sequences   = n(),
    median_length = median(length_bp, na.rm = TRUE),
    min_length    = min(length_bp, na.rm = TRUE),
    max_length    = max(length_bp, na.rm = TRUE),
    n_species     = n_distinct(valid_name),
    n_specimen    = n_distinct(specimen_id),
    n_countries   = n_distinct(iso3),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sequences)) %>%               
  mutate(
    median_length = round(median_length),
    min_length    = round(min_length),
    max_length    = round(max_length)
  )
names(seq_summary)[1] <- "genetic loci"
seq_summary %>%
  slice_head(n = 50) %>%
  write.csv("output/tables/seq_summary_phylo_markers.csv", row.names = F)

################################################################################
#sparkline + summary table
library(dplyr)
library(gt)
library(gtExtras)
library(webshot2)
gene_summary <- seq_annotation %>%
  group_by(gene) %>%
  summarise(
    n_sequences = n(),
    n_species   = n_distinct(valid_name),
    n_basins    = n_distinct(basin_id),
    n_countries = n_distinct(country),
    .groups = "drop"
  ) %>%
  arrange(desc(n_sequences)) %>%
  slice_head(n = 30)


gene_year_counts <- seq_annotation %>%
  count(gene, year, name = "n_seq")

spark_df <- gene_year_counts %>%
  arrange(year) %>%
  group_by(gene) %>%
  summarise(trend = list(n_seq), .groups = "drop")


length_range_df <- seq_annotation %>%
  group_by(gene) %>%
  summarise(
    length_range_text = paste0(median(length_bp, na.rm = TRUE),"(",min(length_bp, na.rm = TRUE), "–", max(length_bp, na.rm = TRUE),")"),
    .groups = "drop"
  )

gene_table <- gene_summary %>%
  left_join(spark_df,        by = "gene") %>%
  left_join(length_range_df, by = "gene") %>%
  select(gene,n_sequences,length_range_text,n_species,n_basins,n_countries,trend)


gt_plot <- gene_table %>%
  gt() %>%
  gt_plt_sparkline(
    trend,
    type = "default",
    same_limit = TRUE
  ) %>%
  gt_fa_rank_change(n_sequences, fa_type = "arrow",font_color = "match",palette = c("#1f77b4", "lightgrey", "#762a83")) %>% 
  gt_fa_rank_change(n_species, fa_type = "angles",font_color = "match") %>% 
  fmt_number(
    columns = c(n_sequences, n_species, n_basins, n_countries),
    sep_mark = ",",
    decimals = 0
  ) %>%
  cols_label(
    gene              = "Genetic loci",
    n_sequences       = "Number of Sequences",
    length_range_text = "Length range (bp)",        
    n_species         = "Number of Species",
    n_basins          = "Number of Basins",
    n_countries       = "Number of Countries",
    trend             = "Annual trend"
  ) %>%
  gt_theme_pff()%>%
  gt_highlight_cols(columns = gene, fill = "#e4e8ec",font_weight = "bold") %>%
  tab_style(
    style = cell_text(transform = "capitalize"),   
    locations = cells_column_labels()
  ) %>%
  tab_source_note(
    source_note = html(
      "<div style='text-align:left;'>
      <em>Note:</em> Loci sorted by number of sequences; data updated to 2024-12.<br>
      <em>Sparkline:</em> annual sequence counts from 1995–2024.
     </div>"
    )
  )


gt_plot

gt::gtsave(
  data = gt_plot,
  filename = "figures/supplement/genetic_loci_summary.html"
)

