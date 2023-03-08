rm(list=ls())
library(dplyr)
library(ggplot2)
library(ggalt)

#load formatted DiveRS data
df<-readRDS(file = here::here("outputs/df_filt_trans.rds"))

#load Diaz data
diaz<-read.csv("data/diaz_species_mean_traits.csv")

#species not in diaz
setdiff(rownames(df),diaz$Species_name_standardized_against_TPL)

#genus of species not in diaz
genera<-sapply(strsplit(setdiff(rownames(df),diaz$Species_name_standardized_against_TPL), " "), "[[", 1)

#diaz only with replacement genera
diaz_genera<-diaz[diaz$Genus%in%genera,]

#find which species of genus has the most information
for(i in 1:length(genera)){
  
  #make tmp df with rows from genus
  diaz_tmp<-diaz_genera[diaz_genera$Genus==genera[i],]
  

  diaz_tmp_max<-diaz_tmp[which.max(diaz_tmp$Number_of_traits_with_values),]
  
  if(i == 1){
    sub_genera<-diaz_tmp_max
  } else {
    sub_genera<-rbind(sub_genera,diaz_tmp_max)
  }
  
}

# Save data frame with substitutes
saveRDS(sub_genera, file = here::here("outputs/sub_genera.rds"))
