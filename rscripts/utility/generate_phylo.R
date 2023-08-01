rm(list=ls())

#load packages
library(ape)
library(V.PhyloMaker2)
library(U.Taxonstand)

#load formatted data
df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

#get species list
spec_list <- rownames(df)

#read in database
#The Plant List (TPL) taken from:
#https://github.com/nameMatch/Database/
db1<-read.csv("data/Plants_TPL_database_part1.csv")
db2<-read.csv("data/Plants_TPL_database_part2.csv")
db3<-read.csv("data/Plants_TPL_database_part3.csv")
db<-rbind(db1,db2,db3)

#get standardized list of taxon names
tpl_list <- nameMatch(spec_list,spSource=db)

#NOTE: There are some species with issues
tpl_list[tpl_list$Fuzzy==1,]
tpl_list[tpl_list$name.dist>1,]

#reformat standardized list of names
spec_df <- tpl_list[,c("Accepted_SPNAME","Genus_in_database","Family")]

#generate phylogenetic tree from species list
#Use scenario 3
maker_tree <- phylo.maker(spec_df)
str(maker_tree)

#make tree object
maker_tree_s3<-maker_tree$scenario.3

#NAMES TO RESOLVE:

#in tree but not DiveRS data
setdiff(gsub("_"," ",maker_tree_s3$tip.label),rownames(df))

#in DiveRS data but not tree
setdiff(rownames(df),gsub("_"," ",maker_tree_s3$tip.label))

#TEMPORARY FIX:
#drop tip
maker_tree_s3<-drop.tip(maker_tree_s3,"Quintinia_apoensis")

#rename species
maker_tree_s3$tip.label[grep("Molineria_latifolia",maker_tree_s3$tip.label)]<-"Curculigo_latifolia"
maker_tree_s3$tip.label[grep("Quintinia_kuborensis",maker_tree_s3$tip.label)]<-"Quintinia_hyehenensis"

#visualize tree
plot(maker_tree_s3,type="fan",cex=0.25)

#check names
rownames(df)==gsub("_"," ",sort(maker_tree_s3$tip.label))

#write tree to file
write.tree(maker_tree_s3,file="outputs/pruned_tree.tre")
