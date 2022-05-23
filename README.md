
<!-- README.md is generated from README.Rmd. Please edit that file -->

# divers-tree

### Content

This repository is structured as follow:

-   [`data/`](https://github.com/ajhelmstetter/sseReview/tree/master/data):
    contains all raw data required to perform analyses

-   [`rscripts/`](https://github.com/ajhelmstetter/sseReview/tree/master/rscripts/):
    contains R scripts to run each step of the workflow

-   [`outputs/`](https://github.com/ajhelmstetter/sseReview/tree/master/outputs):
    contains all the results created during the workflow

-   [`figures/`](https://github.com/ajhelmstetter/sseReview/tree/master/figures):
    contains all the figures created during the workflow

-   [`paper/`](https://github.com/ajhelmstetter/sseReview/tree/master/paper):
    contains all the manuscript and related content (biblio, templates,
    etc.)

-   [`R/`](https://github.com/ajhelmstetter/sseReview/tree/master/R):
    contains R functions developed especially for this project

-   [`man/`](https://github.com/ajhelmstetter/sseReview/tree/master/man):
    contains help files of R functions

-   [`DESCRIPTION`](https://github.com/ajhelmstetter/sseReview/tree/master/DESCRIPTION):
    contains project metadata (author, date, dependencies, etc.)

-   [`make.R`](https://github.com/ajhelmstetter/sseReview/tree/master/make.R):
    master R script to run the entire project by calling each R script
    stored in the `rscripts/` folder

## Workflow

### 1\_proteus\_data\_preparation\_discrete.R

Read in PROTEUS data and recode discrete data based on trait\_recoding
table

### 2\_proteus\_data\_preparation\_quant.R

Read in PROTEUS data and prepare quantiative data

### 3\_recode\_quantitative\_discrete.R

Convert some quantitative variables to discrete ones

### 4\_merge\_subset\_data.R

Combine discrete, quantitative and discretized data into one data frame

### 5\_clean\_filter\_df.R

Clean up data and filter species / traits based on proportion of missing
data

## Other scripts

### discrete\_state\_freqs.R

The number of species with each discrete state (original) from the
PROTEUS dataset.

## Notes

Do we centre and normalise traits?
What about outliers? There are some very extreme values.

