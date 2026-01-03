# options(repos = c(
#   rtrees = 'https://daijiang.r-universe.dev',
#   CRAN = 'https://cloud.r-project.org'))
# install.packages("rtrees")

# Load required packages (suppress startup messages)
suppressPackageStartupMessages({
  library(rtrees)
  library(ape)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(xfun)  # For download_file function used in original code
})

# --- Load and normalize CAS data once (avoid redundant processing) ---
cas <- openxlsx::read.xlsx("input/raw/cas_freshwater_v1.xlsx") %>%
  select(
    species = valid_name,  # Rename for consistency
    family, Order, Class
  ) %>%
  mutate(
    genus = sub("^(\\S+).*", "\\1", species),  # Extract genus from species name (first word)
    species = gsub(" ", "_", species)  # Match tree tip-label convention (underscores instead of spaces)
  ) %>%
  select(species, genus, family, order, class)  # Column order required by add_root_info()


######################## Ray-finned fishes tree #########################
################################################################################
# Rabosky, D. L., Chang, J., Title, P. O., Cowman, P. F., Sallan, L.,          #
# Friedman, M., ... & Alfaro, M. E. (2018). An inverse latitudinal gradient    #
# in speciation rate for marine fishes. Nature, 559(7714), 392-395.            #
################################################################################

# Download ray-finned fish tree file (with missing species inserted via birth-death process)
fishurl2 <- "https://fishtreeoflife.org/downloads/actinopt_full.trees.xz"
tempf <- file.path(getwd(), "input/raw/fishtree", "actinopt_full.trees.xz")  # Storage path
download.file(fishurl2, tempf)  # Download the file

# Read tree collection (100 phylogenetic trees)
tree_fish_32k_n100 <- ape::read.tree(tempf)
# unlink(tempf)  # Uncomment to clean up temporary file if needed
# str(tree_fish_32k_n100)  # Inspect structure (100 trees)
# tree_fish_32k_n100[[1]]  # Single tree contains 31516 species
# Nnode(tree_fish_32k_n100[[1]])  # Check number of nodes


###################### Chondrichthyan fishes tree (sharks & rays) #######################
################################################################################
# Stein, R. W., Mull, C. G., Kuhn, T. S., Aschliman, N. C., Davidson, L. N.,   #
# Joy, J. B., ... & Mooers, A. O. (2018). Global priorities for conserving     #
# the evolutionary history of sharks, rays and chimaeras. Nature Ecology and   #
# Evolution, 2(2), 288-298.                                                    #
################################################################################

# Download chondrichthyan calibrated trees
xfun::download_file(
  url = "https://data.vertlife.org/sharktree/10.cal.tree.nex",
  output = "input/raw/fishtree/shark_10.cal.tree.nex"
)

# Read tree collection (1000 calibrated trees)
tree_shark_ray_n1000 <- ape::read.nexus("input/raw/fishtree/shark_10.cal.tree.nex")
# Nnode(tree_shark_ray_n1000[[1]])  # Single tree contains 1191 species (fixed variable name typo)


############################ Lungfishes tree ############################
################################################################################
# Brownstein, C. D., Harrington, R. C., & Near, T. J. (2023). The              #
# biogeography of extant lungfishes traces the breakup of Gondwana. Journal    #
# of Biogeography, 50(7), 1191-1198.                                           #
################################################################################

# Read lungfish tree
lungfish <- read.nexus("input/raw/fishtree/Brownstein_2023/tipdate/lungfish12SUM.tree")
lungfish[["tip.label"]] <- gsub("(_\\d+\\.?\\d*)$", "", lungfish[["tip.label"]])  # Clean numeric suffixes from labels

# Filter tips matching Dipneusti class species in CAS data
lungfish_tre <- ape::keep.tip(
  lungfish,
  intersect(
    cas[cas$Class == "Dipneusti", "species"],  # Lungfish species from CAS
    lungfish$tip.label  # Species present in the tree
  )
)
# plot(lungfish_tre)  # Uncomment for visualization


############################# Hagfish tree ##############################
################################################################################
# Miyashita, T., Coates, M. I., Farrar, R., Larson, P., Manning, P. L.,        #
# Wogelius, R. A., ... & Currie, P. J. (2019). Hagfish from the Cretaceous     #
# Tethys Sea and a reconciliation of the morphological–molecular conflict in   #
# early vertebrate phylogeny. Proceedings of the National Academy of           #
# Sciences, 116(6), 2146-2151.                                                 #
################################################################################

# Read hagfish tree
hagfish <- read.nexus("input/raw/fishtree/Miyashita_2019/MCC_median_heights.tre")
hagfish[["tip.label"]] <- gsub("(_\\d+\\.?\\d*)$", "", hagfish[["tip.label"]])  # Clean label format

# Filter tips matching Cladistii/Petromyzonti species in CAS data
hagfish_tre <- ape::keep.tip(
  hagfish,
  intersect(
    cas[cas$Class %in% c("Cladistii", "Petromyzonti"), "species"],
    hagfish$tip.label
  )
)
# plot(hagfish_tre)  # Uncomment for visualization


################################## Combine trees ###################################
# Reference: https://www.mun.ca/biology/scarr/Agntha_Phylogeny.html
# Calculate scaling factor to standardize branch lengths between hagfish and lungfish trees
scale_factor <- max(node.depth.edgelength(hagfish_tre)) / max(node.depth.edgelength(lungfish_tre))
lungfish_tre$edge.length <- lungfish_tre$edge.length * scale_factor  # Rescale lungfish branch lengths

# Combine hagfish and lungfish trees into Agnatha (jawless fishes) tree
Agnatha <- bind.tree(hagfish_tre, lungfish_tre)
saveRDS(Agnatha, "input/raw/fishtree/Agnatha.rds")  # Save combined tree


################################## Batch process trees ###################################
Agnatha <- readRDS("input/raw/fishtree/Agnatha.rds")  # Load jawless fishes tree
set.seed(2024)  # Fix random seed for reproducibility
# Randomly sample 100 trees from 1000 chondrichthyan trees
tree_shark_ray_n100 <- tree_shark_ray_n1000[sample(length(tree_shark_ray_n1000), 100)]
tree_fish_32k_n100 <- ape::read.tree("input/raw/fishtree/actinopt_full.trees.xz")  # Re-read ray-finned fish trees
fish_n32k <- read.csv("input/raw/fishtree/fish_n32k.csv")  # Load species name conversion table
names(tree_shark_ray_n100) <- 1:100  # Name sampled chondrichthyan trees

# Process 100 trees in loop, skipping existing results
for (i in 1:100) {
  output_file <- paste0("input/processed/tree_fish/tree_fish_18k_n100_", i, ".rds")
  
  if (!file.exists(output_file)) {  # Only process ungenerated files
    cat("Processing tree", i, "...\n")
    
    # Rescale Agnatha tree to match current chondrichthyan tree time scale
    Agnatha$edge.length <- Agnatha$edge.length * 
      max(node.depth.edgelength(tree_shark_ray_n100[[i]])) / 
      max(node.depth.edgelength(Agnatha))
    
    # Combine chondrichthyan tree with Agnatha tree
    tre1 <- ape::bind.tree(tree_shark_ray_n100[[i]], Agnatha)
    
    # Rescale combined tree to match current ray-finned fish tree time scale
    tre1$edge.length <- tre1$edge.length * 
      max(node.depth.edgelength(tree_fish_32k_n100[[i]])) / 
      max(node.depth.edgelength(tre1))
    
    # Combine with ray-finned fish tree to get complete phylogeny
    tre2 <- ape::bind.tree(tree_fish_32k_n100[[i]], tre1)
    
    # Convert tip labels to valid_name using conversion table
    y <- tibble::tibble(species = tre2$tip.label) %>% 
      left_join(fish_n32k[, 1:2], by = "species")
    tre2$tip.label <- y$valid_name
    
    # Add taxonomic information and generate final tree
    res <- rtrees::add_root_info(
      tre2, 
      unique(fish_n32k[, c("genus", "family", "order", "class")])
    )
    result <- get_one_tree(
      sp_list      = cas,
      tree         = res,
      taxon        = "none",
      scenario     = c("at_basal_node", "random_below_basal"),
      show_grafted = TRUE,
      tree_by_user = TRUE,
      .progress    = "text",
      dt           = TRUE
    )
    result <- rm_stars(result)  # Remove asterisks from labels
    
    # Filter for species present in both CAS data and tree
    sp_in_tree <- intersect(cas$species, result$tip.label)
    if (length(sp_in_tree) == 0) {
      warning("No overlapping species in tree ", i, "; skipping.\n")
      next
    }
    
    # Keep only overlapping species and save result
    tre_keep <- keep.tip(result, sp_in_tree)
    saveRDS(tre_keep, output_file)
    
  } else {
    cat("Tree", i, "already exists; skipping.\n")
  }
}

# 3427 species added at genus level (*) 
# 
# 411 species added at family level (**) 
# 
# 2 species have no co-family species in the mega-tree, skipped
# (if you know their family, prepare and edit species list with `rtrees::sp_list_df()` may help): 
#   Tarumania_walkerae, Aenigmachanna_gollum
