rm(list=ls())

#load libraries
library(tidyverse)
library(gtools)

#import data
df<-read.csv("data/qryDiveRS_Data_2023-05-21.csv")
str(df)

#extract categorical data only
cat_df<-df[grep("(D1)",df$NChrDat),]

#remove extra metadata columns
cat_df<-cat_df[,c(1,2,3,5)]

#pivot to wide format with traits as columns using ID to keep separate multiple records per species
df.wide <- pivot_wider(cat_df, names_from = NChrDat, values_from = c(NCstDat),id_cols = c(NDat,NTaxDat))
head(df.wide)

#get rid of NDat species code column
df.wide<-df.wide[,-1]

#keep only distinct values per trait per species
df.wide<-distinct(df.wide)

#categories that have been chosen for analysis
chosen_cats<-c(
  "1. Habit (D1)",
  "4. Plant sexual system (D1)",
  "5. Life history (D1)",
  "80. Phenotypic mating system (D1)",
  "83. Self-incompatibility system (genetic) (D1)",
  "90. Pollination syndrome (D1)",
  "95. Dispersal syndrome (D1)",
  "100. Floral structural sex (D1)",
  "102. Ovary position (D1)",
  "107. Floral reward (D1)",
  "207. Symmetry of perianth (D1)",
  "239. Main color of perianth at anthesis (D1)",
  "660. Inflorescence attractive organ (D1)",
  "661. Inflorescence attractive color (D1)",
  "15562. Fruit fleshiness (D1)"
)

#subset dataframe to include only those chosen categories
df.wide.chosen<-df.wide[,chosen_cats]

#generate table for all states per trait
for(i in 1:length(colnames(df.wide.chosen))) {
  if (i == 1) {
    trait_codes<-cbind(
      rep(colnames(df.wide.chosen)[i], length(table(df.wide.chosen[, i]))),
      names(table(df.wide.chosen[, i])))
  } else {
    trait_codes<-rbind(trait_codes,cbind(
      rep(colnames(df.wide.chosen)[i], length(table(df.wide.chosen[, i]))),
      names(table(df.wide.chosen[, i]))))
  }
}

#output table of all of the PROTEUS states for each trait (to help construct recoding table)
write.csv(trait_codes,"outputs/one_hot_1_all_states_per_trait.csv")

#subset dataframe to include only those chosen categories and species names
df.wide.chosen<-df.wide[,c("NTaxDat",chosen_cats)]

#rename columns to proteus trait number
colnames(df.wide.chosen)[2:length(colnames(df.wide.chosen))]<-gsub("\\..*","",colnames(df.wide.chosen)[2:length(colnames(df.wide.chosen))])

#change back to data frame from tibble
df.wide.chosen<-as.data.frame(df.wide.chosen)

#remove trait names before trait values
for(i in 2:length(colnames(df.wide.chosen))){
  df.wide.chosen[,i]<-gsub(".*: ","",df.wide.chosen[,i])
}

#read in recoding table
df_recode<-read.csv("data/trait_recoding - Categorical to categorical.csv")

#empty list to store traits
trait_list<-list()

#reset i
i<-1

# loop through traits to build trait list
for(i in 1:length(unique(df_recode$new_trait))){
  
  #get df for trait
  trait_df<-df_recode[df_recode$new_trait==unique(df_recode$new_trait)[i],]
  
  #get trait number for new trait
  trait_no<-unique(trait_df$trait_number)
  trait_no
  
  #combine all old traits into one data frame
  for(j in 1:length(trait_no)){
    if(j == 1){
      new_trait_df<-df.wide.chosen[,c("NTaxDat",as.character(trait_no[j]))]
    } else {
      tmp<-df.wide.chosen[,c("NTaxDat",as.character(trait_no[j]))]
      colnames(tmp)<-colnames(new_trait_df)
      new_trait_df<-rbind(new_trait_df,tmp)
    }
  }
  
  #omit na values
  new_trait_df<-na.omit(new_trait_df)
  
  #select only old states of interest
  new_trait_df<-new_trait_df[new_trait_df[,2]%in%trait_df$old_state,]
  
  #merge, adding new row that recodes old states to new states
  trait_df_merged<-merge(trait_df,new_trait_df,by.x = 'old_state', by.y = 2, all.x = T)
  
  #rename output column based on the new trait to fix names downstream
  colnames(trait_df_merged)[6] <- unique(trait_df_merged$new_trait)
  
  #reduce to species/state columns only
  trait_list[[i]] <- trait_df_merged[,c(6,7)]  
  
  ###
  # duplicate combinations may arise when trait recoding acts on multiple states for the same species
  # e.g. Viburnum rufidulum is coded as a tree and a shrub in PROTEUS leading to two duplicate 'woody' states in the trait list 
  ###
  
  #remove duplicate combinations
  trait_list[[i]]<-unique(trait_list[[i]])
  trait_list[[i]][,2]=as.character(trait_list[[i]][,2])
  
  #one-hot encoding
  trait_list[[i]] <- dcast(data=melt(trait_list[[i]],id.vars="NTaxDat"), NTaxDat ~ variable + value, length)
  
}

str(trait_list)

#Merge dataframes in list
for(i in 1:length(trait_list)){
  if(i == 1){
    disc_df<-trait_list[[i]]
  } else {
    disc_df<-merge(disc_df,trait_list[[i]],by.x = 'NTaxDat', by.y = 'NTaxDat', all = T)
  }
}

#output csv for downstream use
write.csv(disc_df,"outputs/one_hot_1_proteus_discrete_recoded.csv")


