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

source(here::here("rscripts", "1_proteus_data_preparation_discrete.R"))
source(here::here("rscripts", "2_proteus_data_preparation_quant.R"))
source(here::here("rscripts", "3_recode_quantitative_discrete.R"))
source(here::here("rscripts", "4_merge_subset_data.R"))
source(here::here("rscripts", "5_clean_filter_df.R"))
source(here::here("rscripts", "6_scale_transform.R"))
source(here::here("rscripts", "7_correlation.R"))
source(here::here("rscripts", "8_pcoa.R"))
source(here::here("rscripts", "8.1_pcoa_no_reproductive.R"))
source(here::here("rscripts", "9_dimensonality_analyses.R"))
source(here::here("rscripts", "10_functional_space_mfd.R"))
source(here::here("rscripts", "11_clustering_hierarchical.R"))
source(here::here("rscripts", "11.1_clustering_kprototype.R"))
source(here::here("rscripts", "11.2_clustering_PAM.R"))
source(here::here("rscripts", "11.3_clustering_LCM.R"))