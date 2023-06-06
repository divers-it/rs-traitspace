#load library
library(tidyverse)
library(gtools)

#import data
dg<-read.csv("data/qryDiveRS_Data_2022-12-02.csv")
str(dg)

#filter by quality
#TO DO

#extract categorical data only
cat_dg<-dg[grep("(D1)",dg$NChrDat),]

#remove extra metadata columns
cat_dg<-cat_dg[,c(1,2,3,5)]

#pivot to wide format with traits as columnsm using ID to keep separate multiple records per species
dg.wide <- pivot_wider(cat_dg, names_from = NChrDat, values_from = c(NCstDat),id_cols = c(NDat,NTaxDat))
head(dg.wide)

#get rid of id column
dg.wide<-dg.wide[,-1]

#keep only distinct values per trait per species
dg.wide<-distinct(dg.wide)

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
dg.wide.chosen<-dg.wide[,chosen_cats]

#generate table for all states per trait
for(i in 1:length(colnames(dg.wide.chosen))) {
  if (i == 1) {
    trait_codes<-cbind(
      rep(colnames(dg.wide.chosen)[i], length(table(dg.wide.chosen[, i]))),
      names(table(dg.wide.chosen[, i])))
  } else {
    trait_codes<-rbind(trait_codes,cbind(
      rep(colnames(dg.wide.chosen)[i], length(table(dg.wide.chosen[, i]))),
      names(table(dg.wide.chosen[, i]))))
  }
}

write.csv(trait_codes,"outputs/all_states_per_trait_one_hot.csv")

#subset dataframe to include only those chosen categories and species names
#what about ID?
dg.wide.chosen<-dg.wide[,c("NTaxDat",chosen_cats)]

#rename columns to proteus trait number
colnames(dg.wide.chosen)[2:length(colnames(dg.wide.chosen))]<-gsub("\\..*","",colnames(dg.wide.chosen)[2:length(colnames(dg.wide.chosen))])

#change back to data frame from tibble
dg.wide.chosen<-as.data.frame(dg.wide.chosen)

#remove trait names before trait values
#for(i in 2:length(colnames(dg.wide.chosen))){
#  dg.wide.chosen[,i]<-gsub(".*\\. ","",gsub("\\(.*\\) ","",dg.wide.chosen[,i]))
#}
#remove trait names before trait values
for(i in 2:length(colnames(dg.wide.chosen))){
  dg.wide.chosen[,i]<-gsub(".*: ","",dg.wide.chosen[,i])
}


library(reshape2)

#read in recoding table
#dg_recode<-read.csv("data/trait_recoding - Categorical to more categorical.csv")
dg_recode<-read.csv("data/trait_recoding - Categorical to categorical.csv")

#empty list to store traits
trait_list<-list()

i<-1

for(i in 1:length(unique(dg_recode$new_trait))){

  #get dg for trait
  trait_dg<-dg_recode[dg_recode$new_trait==unique(dg_recode$new_trait)[i],]

  #get trait number for new trait
 trait_no<-unique(trait_dg$trait_number)
  trait_no

  #combine all old traits into one data frame
  for(j in 1:length(trait_no)){
    if(j == 1){
      new_trait_dg<-dg.wide.chosen[,c("NTaxDat",as.character(trait_no[j]))]
    } else {
      tmp<-dg.wide.chosen[,c("NTaxDat",as.character(trait_no[j]))]
      colnames(tmp)<-colnames(new_trait_dg)
      new_trait_dg<-rbind(new_trait_dg,tmp)
    }
  }

  #omit na values
  new_trait_dg<-na.omit(new_trait_dg)

   #select only old states of interest
  new_trait_dg<-new_trait_dg[new_trait_dg[,2]%in%trait_dg$old_state,]

  merge(trait_dg,new_trait_dg,by.x = 'old_state', by.y = 2, all.x = T)

  #merge, adding new row that recodes old states to new states
  trait_dg_merged<-merge(trait_dg,new_trait_dg,by.x = 'old_state', by.y = 2, all.x = T)
  trait_dg_merged

  colnames(trait_dg_merged)[6]=unique(trait_dg_merged$new_trait)

  #reduce to species/state columns only
  trait_dg_merged=trait_dg_merged[,c(6,7)]  

  trait_list[[i]]<-trait_dg_merged
  
  ###
  # Figure out why this is necessary
  ###
  #remove duplicate combinations
  trait_list[[i]]<-unique(trait_list[[i]])
  trait_list[[i]][,2]=as.character(trait_list[[i]][,2])

  #one-hot encoding

  trait_list[[i]]=dcast(data=melt(trait_list[[i]],id.vars="NTaxDat"), NTaxDat ~ variable + value, length)
  
  #name list
  #names(trait_list)[i]<-unique(trait_dg_merged$new_trait)
  
}

#trait_list[[1]]=data.frame(NTaxDat=unique(dg.wide.chosen$NTaxDat))


str(trait_list)

#Merge dataframes in list
for(i in 1:length(trait_list)){
  if(i == 1){
    disc_dg<-trait_list[[i]]
  } else {
    disc_dg<-merge(disc_dg,trait_list[[i]],by.x = 'NTaxDat', by.y = 'NTaxDat', all = T)
  }
}

#fix column names

#colnames(disc_dg)[2:length(colnames(disc_dg))]<-gsub(".*_","",colnames(disc_dg)[2:length(colnames(disc_dg))])



write.csv(disc_dg,"outputs/proteus_discrete_one_hot.csv")


