 ## **Global Freshwater Fish Knowledge Shortfalls**


A reproducible R-based workflow to quantify global Linnaean, Wallacean, and Darwinian knowledge shortfalls in freshwater fishes and to inform strategic biodiversity collection priorities.

[![R](https://img.shields.io/badge/language-R-blue.svg)](https://www.r-project.org/)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC_BY_4.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)
[![GitHub last commit](https://img.shields.io/github/last-commit/Otoliths/FW_Fish_Shortfall)](https://github.com/Otoliths/FW_Fish_Shortfall)
[![Zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18136683.svg)](https://doi.org/10.5281/zenodo.18136683)
[![figshare DOI](https://img.shields.io/badge/figshare-10.6084%2Fm9.figshare.29262098-orange.svg)](https://doi.org/10.6084/m9.figshare.29262098)


### &#128230; **Code and Data Availability**

All R scripts used for data acquisition, processing, analysis, and model fitting in this study are openly available via [GitHub Repository](https://github.com/Otoliths/FW_Fish_Shortfall) and archived on [Zenodo](https://doi.org/10.5281/zenodo.18136683) for long-term preservation.


### &#128193; **Folder Descriptions**
- **code/**: Contains all R scripts used for data acquisition, preprocessing, statistical analyses, and model fitting.
- **input/**: Stores all datasets used in the analyses, including raw input data (`raw/`), intermediate processed data (`processed/`), and analysis-ready datasets (`data_prep/`).
- **output/**: Contains analysis outputs, including fitted model objects (`model/`), summary tables (`tables/`), and execution logs (`logs/`).
- **figures/**: Includes all figures generated in this study, organized into main-text figures (`main/`) and supplementary figures (`supplement/`).
- **report/**: Contains R scripts and R Markdown files used to generate all manuscript figures, including `Figure_1–4.R` and supplementary figure scripts (`Figure_S01–S19.R`).
- **FW_Fish_Shortfall.Rproj**: RStudio project file used to manage the project structure, file paths, and reproducible workflow.

*Reminder: Download Input & Output Data*

Before running any analysis scripts, download the **input/** and **output/** directories from [Zenodo repository](https://doi.org/10.5281/zenodo.18136683)

This repository is organized as follows:

```text
/project_root
├── /code                                # Core analysis scripts, to be executed sequentially for full reproducibility.
│
│ ├── 00_Get_fishtree.R                  # Retrieve and preprocess the global freshwater fish phylogeny.
│ ├── 00_Get_futurestream.R              # Download and prepare FutureStreams river network data.
│ ├── 00_Get_gbif.R                      # Query and download occurrence records from GBIF.
│ ├── 00_Get_ncbi.R                      # Retrieve sequence metadata and records from NCBI.
│ ├── 00_Get_obis.R                      # Download aquatic occurrence data from OBIS.
│
│ ├── 01_Country_sp.R                    # Compile country-level species checklists.
│ ├── 01_Futurestream_clean.R            # Clean and standardize FutureStreams river attributes.
│ ├── 01_Occurrence_clean.R              # Quality control and harmonization of occurrence records.
│
│ ├── 02_Description.R                   # Quantify taxonomic description coverage.
│ ├── 02_Geolocation.R                   # Assess georeferencing completeness and precision.
│ ├── 02_Sequence.R                      # Summarize molecular sequence availability.
│ ├── 02_Seq_summary_phylo_markers.R     # Summarize sequencing effort across phylogenetic markers.
│
│ ├── 03_Basin_discharge.R               # Basin-level mean river discharge.
│ ├── 03_Basin_elevation.R               # Basin-level elevation statistics.
│ ├── 03_Basin_human_density.R           # Human population density aggregated by basin.
│ ├── 03_Basin_latitude.R                # Latitudinal position of drainage basins.
│ ├── 03_Basin_preserved_specimen.R      # Preserved specimen counts per basin.
│ ├── 03_Basin_range_size.R              # Species geographic range size at the basin scale.
│ ├── 03_Basin_rarity.R                  # Basin-level species rarity metrics.
│ ├── 03_Basin_sampling_effort.R         # Occurrence-based sampling effort per basin.
│ ├── 03_Basin_sequence_effort.R         # Molecular sequencing effort per basin.
│ ├── 03_Basin_watertemp.R               # Basin-level water temperature.
│ ├── 03_Watershed area.R                # Drainage basin area calculations.
│
│ ├── 04_Country_discharge.R             # Country-level river discharge.
│ ├── 04_Country_elevation.R             # Country-level elevation.
│ ├── 04_Country_human_density.R         # Human population density at the country level.
│ ├── 04_Country_latitude.R              # Country centroid latitude.
│ ├── 04_Country_preserved_specimen.R    # Preserved specimen counts per country.
│ ├── 04_Country_range_size.R            # Species range size summarized by country.
│ ├── 04_Country_rarity.R                # Country-level species rarity.
│ ├── 04_Country_sampling_effort.R       # Occurrence sampling effort by country.
│ ├── 04_Country_sequence_effort.R       # Sequencing effort summarized by country.
│ ├── 04_Country_watertemp.R             # Country-level water temperature.
│
│ ├── 05_Body_size.R                     # Compile and standardize species body size traits.
│ ├── 05_Taxonomic_activity.R            # Metrics of taxonomic activity effort.
│ ├── 05_Taxonomic_effort.R              # Metrics of taxonomic research effort.
│
│ ├── 06_Basin_linnaean_model.R          # Statistical model for basin-level Linnaean shortfalls.
│ ├── 06_Basin_wallacean_model.R         # Statistical model for basin-level Wallacean shortfalls.
│ ├── 06_Basin_darwinian_model.R         # Statistical model for basin-level Darwinian shortfalls.
│                             
│ ├── 06_Country_linnaean_model.R        # Country-level Linnaean shortfall model.
│ ├── 06_Country_wallacean_model.R       # Country-level Wallacean shortfall model.
│ ├── 06_Country_darwinian_model.R       # Country-level Darwinian shortfall model.
│
│ ├── 07_Basin_linnaean_shortfall.R      # Quantification of basin-level Linnaean shortfalls.
│ ├── 07_Basin_wallacean_shortfall.R     # Quantification of basin-level Wallacean shortfalls.
│ ├── 07_Basin_darwinian_shortfall.R     # Quantification of basin-level Darwinian shortfalls.
│                             
│ ├── 07_Country_linnaean_shortfall.R    # Country-level Linnaean shortfalls.
│ ├── 07_Country_wallacean_shortfall.R   # Country-level Wallacean shortfalls.
│ ├── 07_Country_darwinian_shortfall.R   # Country-level Darwinian shortfalls.
│
│ ├── Dawinian_deficits_basin.R          # Aggregated Darwinian knowledge deficits at the basin scale.
│ ├── Dawinian_deficits_country.R        # Aggregated Darwinian knowledge deficits at the country scale.
│
│ ├── Model_report.R                     # Automated model diagnostics and summary tables.
│
│ ├── Cost_Taxonomy.R                    # Cost estimation for taxonomic research.
│ ├── Cost_field_sampling.R              # Cost estimation for field sampling campaigns.
│ ├── Cost_sequencing.R                  # Cost estimation for molecular sequencing.
│ ├── Cost_all.R                         # Integrated cost assessment across knowledge gaps.
│ ├── Shortfall_reductions.R             # Scenario-based evaluation of shortfall reduction strategies.
│
│ ├── functions/                         # Custom helper functions used throughout the pipeline.
│ │   ├── clean_species_occ.R            # Species-level occurrence data cleaning functions.
│ │   ├── TAE_function.R                 # Functions related to Taxonomic Effort.
│ │   ├── xxx.R                          # Additional utility functions.
│ │   └── zzz.R                          # Additional utility functions.
│ 
│ └── exec/
│     ├── wget.exe                       # External binary used for automated data downloads.
│
├── /input
│ ├── raw/                               # Raw, unmodified input data.
│ ├── processed/                         # Intermediate processed datasets.
│ ├── data_prep/                         # Data prepared for modeling and visualization.
│
├── /output
│ ├── model/                             # Fitted model objects and predictions.
│ ├── logs/                              # Execution logs and diagnostics.
│ ├── tables/                            # Final tables used in the manuscript.
│
├── /figures
│ ├── main/                              # Figures included in the main text.
│ ├── supplement/                        # Supplementary figures.
│
├── /report
│ ├── Figure_1.R / Figure_1.Rmd          # Script and R Markdown for generating Figure 1.
│ ├── Figure_2.R / Figure_2.Rmd          # Script and R Markdown for generating Figure 2.
│ ├── Figure_3.R / Figure_3.Rmd          # Script and R Markdown for generating Figure 3.
│ ├── Figure_4.R / Figure_4.Rmd          # Script and R Markdown for generating Figure 4.
│ ├── Figure_1.html                      # Rendered HTML output for Figure 1.
│ ├── Figure_2.html                      # Rendered HTML output for Figure 2.
│ ├── Figure_3.html                      # Rendered HTML output for Figure 3.
│ ├── Figure_4.html                      # Rendered HTML output for Figure 4.
│ ├── Figure_S01-S19.R                   # Scripts for supplementary figures.
│
└── FW_Fish_Shortfall.Rproj              # RStudio project file managing paths and workflow.

```

### &#128206; **[Supporting Information](https://doi.org/10.6084/m9.figshare.29262098)**

- **Table_Extended_S1.xlsx**: Inventory of 18,821 freshwater fish species included in the analysis (updated 24 November 2024).
- **Table_Extended_S2.xlsx**: Body-size data flags distinguishing measured and imputed values.
- **Table_Extended_S3.xlsx**: Inventory of freshwater fish biogeographical deficiency units identified in this study.
- **Table_Extended_S4.xlsx**: Drainage basins prioritized for future freshwater biodiversity data collection.


### &#128214; If you use this code or data, please cite:

Ding, L., *et al.* Global freshwater fish biodiversity knowledge shortfalls inform strategic collection priorities. *[Journal]* (Year).

Ding, L. (2026). Global Freshwater Fish Knowledge Shortfalls (Version 1.0). Zenodo. https://doi.org/10.5281/zenodo.18136683




