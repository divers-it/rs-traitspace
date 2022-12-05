#load library
library(tidyverse)
library(gtools)

#import data
df<-read.csv("data/qryDiveRS_Data_2022-12-02.csv")
str(df)

#filter by quality
#TO DO

#extract categorical data only
cat_df<-df[grep("(D1)",df$NChrDat),]

#remove extra metadata columns
cat_df<-cat_df[,c(1,2,3,5)]

#pivot to wide format with traits as columnsm using ID to keep separate multiple records per species
df.wide <- pivot_wider(cat_df, names_from = NChrDat, values_from = c(NCstDat),id_cols = c(NDat,NTaxDat))
head(df.wide)

#get rid of id column
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

write.csv(trait_codes,"outputs/all_states_per_trait.csv")

#subset dataframe to include only those chosen categories and species names
#what about ID?
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

i<-1

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

  merge(trait_df,new_trait_df,by.x = 'old_state', by.y = 2, all.x = T)

  #merge, adding new row that recodes old states to new states
  trait_df_merged<-merge(trait_df,new_trait_df,by.x = 'old_state', by.y = 2, all.x = T)

  #reduce to species/state columns only
  trait_list[[i]]<-trait_df_merged[,c('NTaxDat','new_state')]
  
  ###
  # Figure out why this is necessary
  ###
  #remove duplicate combinations
  trait_list[[i]]<-unique(trait_list[[i]])
  trait_list[[i]][,2]=as.character(trait_list[[i]][,2])
  #get names of polymorphic species
  poly_sp<-names(which(table(trait_list[[i]]$NTaxDat)>1))
  
  #loop through polymorphic species
  for(j in 1:length(poly_sp)){
    
    
    #combine multiple states into one
    poly_state<-paste(sort(trait_list[[i]][trait_list[[i]]$NTaxDat%in%poly_sp[j],]$new_state),collapse = "_")
    
    #replace both states with combined states
    trait_list[[i]]$new_state[grep(poly_sp[j],trait_list[[i]]$NTaxDat)]<-poly_state
    
  }
  
  #remove duplicate combinations
  trait_list[[i]]<-unique(trait_list[[i]])
  
  #name list
  names(trait_list)[i]<-unique(as.character(trait_df_merged$new_trait))
  
}

str(trait_list)

#Merge dataframes in list
for(i in 1:length(trait_list)){
  if(i == 1){
    disc_df<-trait_list[[i]]
  } else {
    disc_df<-merge(disc_df,trait_list[[i]],by.x = 'NTaxDat', by.y = 'NTaxDat', all.x = T)
    colnames(disc_df)[i+1]<-names(trait_list)[i]
  }
}

#fix column name
#colnames(disc_df)[2]<-"Woodiness"

write.csv(disc_df,"outputs/proteus_discrete_recoded.csv")


