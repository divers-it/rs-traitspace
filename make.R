# ' @title rs-traitspace: A Research Compendium
# ' 
# ' @description 
# ' A paragraph providing a full description of the project and describing each 
# ' step of the workflow.
# ' 
# ' @author Andrew Helmstetter \email{andrew.j.helmstetter@gmail.com}
# ' 
# ' @date 2024/07/25

# #  Install Dependencies (listed in DESCRIPTION) ----

devtools::install_deps(upgrade = "never")

# #  Load Project Addins (R Functions and Packages) ----

devtools::load_all()

# #  Global Variables ----

#  You can list global variables here (or in a separate R script)

# #  Run Project ----

#  List all R scripts in a sequential order and using the following form:
#  source(here::here("rscripts", "script_X.R"))

# Standard data set preparation
source(here::here("rscripts", "1_proteus_data_preparation_discrete.R"))
source(here::here("rscripts", "2_proteus_data_preparation_quant.R"))
source(here::here("rscripts", "3_recode_quantitative_discrete.R"))
source(here::here("rscripts", "4_merge_subset_data.R"))
source(here::here("rscripts", "5_clean_filter_data.R"))
source(here::here("rscripts", "6_scale_transform.R"))

# One-hot data set preparation
source(here::here("rscripts", "one_hot_1_proteus_data_preparation_discrete.R"))
source(here::here("rscripts", "one_hot_3_recode_quantitative_discrete.R"))
source(here::here("rscripts", "one_hot_4_merge_subset_data.R"))
source(here::here("rscripts", "one_hot_5_clean_filter_data.R"))
source(here::here("rscripts", "one_hot_6_scale_transform.R"))

# utility scripts to build a tree and use this to impute missing data
# NOTE: This can have issues, rerun if errors occur
source(here::here("rscripts/utility", "generate_phylo.R"))
source(here::here("rscripts/utility", "impute_missing_data.R"))

# NOTE: This can have issues, rerun if errors occur
# NOTE: Contains temporary fixes for errors because of synonyms
source(here::here("rscripts/utility", "get_taxonomy.R"))

# look at correlations between traits in data sets
source(here::here("rscripts", "7_correlation.R"))

# Principal Coordinates Analyses
# NOTE: There seems to be an error downloading phylopics
source(here::here("rscripts", "8.0_pcoa.R"))

# Trait space quality
# NOTE: long to run
source(here::here("rscripts", "9_dimensionality_analyses.R"))

# Comparison with Diaz et al. 2022 data set
# NOTE: long to run and needs ~6GB RAM
source(here::here("rscripts", "8.1_pcoa_diaz_data.R"))

# Clustering
source(here::here("rscripts", "10_clustering_PAM.R"))

# Supplementary clustering methods
source(here::here("rscripts", "supplement/10.1_clustering_kprototype.R"))
source(here::here("rscripts", "supplement/10.2_clustering_hierarchical.R"))

# One-hot clustering
source(here::here("rscripts", "one_hot_10_clustering_PAM.R"))
source(here::here("rscripts", "supplement/one_hot_10.1_clustering_kprototype.R"))
source(here::here("rscripts", "supplement/one_hot_10.2_clustering_hierarchical.R"))
 
# Bring clustering together results
source(here::here("rscripts/utility", "collate_clustering_results.R"))

# Calculate functional indices
source(here::here("rscripts", "11_functional_space_mfd.R"))

# Make phylogenetic tree plot
source(here::here("rscripts", "12_phylo_plot.R"))

# Currently unable to source due to memory error
source(here::here("rscripts", "one_hot_13_loadings.R"))

# UMAP analyses and PCoAs with all traits
source(here::here("rscripts", "14_UMAP.R"))

# Simulating neutral data
# NOTE: long to run
source(here::here("rscripts", "15_simulations.R"))

# Predicting flower sex and mating system using random forest
source(here::here("rscripts", "16_prediction.R"))
