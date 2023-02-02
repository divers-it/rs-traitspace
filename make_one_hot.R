#' @title one-hot and loadings
#' 
#' @description 
#' A paragraph providing a full description of the project and describing each 
#' step of the workflow.
#' 
#' @author Andrew Helmstetter \email{andrew.j.helmstetter@gmail.com}
#' 
#' @date 2021/04/15

#Modified for one-hot encoding allowing trait loadings by Jos Käfer (jos.kafer@cnrs.fr) on 2023/02/02

## Install Dependencies (listed in DESCRIPTION) ----

devtools::install_deps(upgrade = "never")


## Load Project Addins (R Functions and Packages) ----

devtools::load_all()


## Global Variables ----

# You can list global variables here (or in a separate R script)

## Run Project ----

# List all R scripts in a sequential order and using the following form:
# source(here::here("rscripts", "script_X.R"))

source(here::here("rscripts", "1_proteus_data_preparation_discrete_one_hot.R"))
source(here::here("rscripts", "2_proteus_data_preparation_quant_one_hot.R"))
source(here::here("rscripts", "3_recode_quantitative_discrete.R"))
source(here::here("rscripts", "4_merge_subset_data_one_hot.R"))
source(here::here("rscripts", "5_clean_filter_df_one_hot.R"))
source(here::here("rscripts", "6_scale_transform_one_hot.R"))
source(here::here("rscripts", "8_pcoa_one_hot.R"))
source(here::here("rscripts", "11_clustering_hierarchical_one_hot.R"))

