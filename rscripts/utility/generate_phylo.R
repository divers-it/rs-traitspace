rm(list=ls())

# load/install packages
library(ape)
library(TNRS)
if(!require(V.PhyloMaker2)){
  devtools::install_github("jinyizju/V.PhyloMaker2")
}
library(V.PhyloMaker2)

# load formatted data
df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

# get species list
spec_list <- rownames(df)

####
## U.Taxonstand names ----
####

# NOTE: Not run as results are the same as TNRS
# # load/install package
#
# if(!require(U.Taxonstand)){
#   devtools::install_github("ecoinfor/U.Taxonstand")
# }
# library(U.Taxonstand)
#
# # read in database
# # The Plant List (TPL) taken from:
# # https://github.com/nameMatch/Database/
# # NOTE: THESE WILL NEED TO BE DOWNLOADED IF YOU WANT TO RUN U.Taxonstand. 
# db1<-read.csv("data/Plants_WCVP_database_part1.csv")
# db2<-read.csv("data/Plants_WCVP_database_part2.csv")
# db3<-read.csv("data/Plants_WCVP_database_part3.csv")
# db<-rbind(db1,db2,db3)
# 
# # get standardized list of taxon names
# wcvp_list <- nameMatch(spec_list,spSource=db)
# 
# # NOTE: There could be some species with issues
# wcvp_list[wcvp_list$Fuzzy==1,]
# wcvp_list[wcvp_list$name.dist>1,]
# 
# # reformat standardized list of names
# spec_df <- wcvp_list[,c("Accepted_SPNAME","Genus_in_database","Family")]
# 
# # generate phylogenetic tree from species list
# # Use scenario 3
# maker_tree <- phylo.maker(spec_df)
# maker_tree

####
## TNRS names ----
####

Sys.sleep(2)
# run TNRS to check species (best result only)
check_species <- TNRS::TNRS(spec_list, matches="best", sources="wcvp")
Sys.sleep(2)

# how many name issues?
table(check_species$Name_submitted == check_species$Accepted_species)

# species with issues
issue_species <- check_species[check_species$Name_submitted != check_species$Accepted_species,]

# examine species
issue_species[,c("Name_submitted",
                 "Accepted_name")]
# NOTE: verified on POWO: https://powo.science.kew.org/

# Number of families according to TNRS ----
unique(spec_df$Accepted_family)

# reformat standardized list of names
spec_df <- check_species[,c("Accepted_species","Genus_matched","Accepted_family")]
maker_tree_tnrs <- phylo.maker(spec_df)

# NOT RUN:
# Are TNRS edges the same as U.Taxonstand? (YES)
# table(maker_tree$scenario.3$edge.length == maker_tree_tnrs$scenario.3$edge.length)

# make tree object
maker_tree_s3<-maker_tree_tnrs$scenario.3

### Check if any names need to be resolved ----

# in tree but not DiveRS data
setdiff(gsub("_"," ",maker_tree_s3$tip.label),rownames(df))

# in DiveRS data but not tree
setdiff(rownames(df),gsub("_"," ",maker_tree_s3$tip.label))

# visualize tree
plot(maker_tree_s3,type="fan",cex=0.25)

# check names match
rownames(df)==gsub("_"," ",sort(maker_tree_s3$tip.label))

# NOT RUN: check matching of names
# View(cbind(rownames(df),gsub("_"," ",sort(maker_tree_s3$tip.label))))

# write tree to file
write.tree(maker_tree_s3,file="outputs/pruned_tree.tre")
