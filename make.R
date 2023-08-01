#' @title sseReview: A Research Compendium
#' 
#' @description 
#' A paragraph providing a full description of the project and describing each 
#' step of the workflow.
#' 
#' @author Andrew Helmstetter \email{andrew.j.helmstetter@gmail.com}
#' 
#' @date 2021/04/15

## Install Dependencies (listed in DESCRIPTION) ----

devtools::install_deps(upgrade = "never")

## Load Project Addins (R Functions and Packages) ----

devtools::load_all()

## Global Variables ----

# You can list global variables here (or in a separate R script)

## Run Project ----

# List all R scripts in a sequential order and using the following form:
# source(here::here("rscripts", "script_X.R"))

#Standard data set preparation
source(here::here("rscripts", "1_proteus_data_preparation_discrete.R"))
source(here::here("rscripts", "2_proteus_data_preparation_quant.R"))
source(here::here("rscripts", "3_recode_quantitative_discrete.R"))
source(here::here("rscripts", "4_merge_subset_data.R"))
source(here::here("rscripts", "5_clean_filter_data.R"))
source(here::here("rscripts", "6_scale_transform.R"))
source(here::here("rscripts", "7_correlation.R"))

#One-hot data set preparation
source(here::here("rscripts", "one_hot_1_proteus_data_preparation_discrete.R"))
source(here::here("rscripts", "one_hot_3_recode_quantitative_discrete.R"))
source(here::here("rscripts", "one_hot_4_merge_subset_data.R"))
source(here::here("rscripts", "one_hot_5_clean_filter_data.R"))
source(here::here("rscripts", "one_hot_6_scale_transform.R"))

#utility scripts
source(here::here("rscripts/utility", "generate_phylo.R"))
source(here::here("rscripts/utility", "impute_missing_data.R"))

#NOTE: Contains temporary fixes for errors because of synonyms
source(here::here("rscripts/utility", "get_taxonomy.R"))

#Principal Coordinates Analyses
source(here::here("rscripts", "8_pcoa.R"))
source(here::here("rscripts", "8.1_pcoa_no_reproductive.R"))
source(here::here("rscripts", "8.2_pcoa_no_vegetative.R"))
source(here::here("rscripts", "8.3_pcoa_diaz_data.R"))
source(here::here("rscripts", "8.4_pcoa_diaz_data_shared.R"))

#Clustering
source(here::here("rscripts", "11_clustering_hierarchical.R"))
source(here::here("rscripts", "11.1_clustering_kprototype.R"))
source(here::here("rscripts", "11.2_clustering_PAM.R"))
source(here::here("rscripts", "11.3_clustering_LCM.R"))
source(here::here("rscripts", "11.4_clustering_densityClust.R"))

#One-hot clustering
source(here::here("rscripts", "one_hot_11_clustering_hierarchical.R"))
source(here::here("rscripts", "one_hot_11.1_clustering_kprototype.R"))
source(here::here("rscripts", "one_hot_11.2_clustering_PAM.R"))
source(here::here("rscripts", "one_hot_11.3_clustering_LCM.R"))
source(here::here("rscripts", "one_hot_11.4_clustering_densityClust.R"))

#Bring clustering together results
source(here::here("rscripts/utility", "collate_clustering_results.R"))

#Trait space quality
source(here::here("rscripts", "9_dimensonality_analyses.R"))
source(here::here("rscripts", "10_functional_space_mfd.R"))

#Phylogenetic analyses
source(here::here("rscripts", "phylo_signal.R"))
source(here::here("rscripts", "ASR_corHMM.R"))

#Currently unable to source due to memory error
source(here::here("rscripts", "one_hot_13_loadings.R"))
