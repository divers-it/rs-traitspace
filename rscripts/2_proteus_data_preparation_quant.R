rm(list=ls())

#load libraries
library(tidyverse)
library(gtools)

#import data
df<-read.csv("data/qryDiveRS_Data_2023-05-21.csv")
str(df)

#extract quant data only
q_df<-df[grep("(C1)",df$NChrDat),]

#remove extra metadata columns
q_df<-q_df[,c(1,2,3,6:8)]

#calculate mean value (VAL,MIN,MAX) for each trait for each species
summ_q_df<-q_df %>% group_by(NTaxDat,NChrDat) %>% summarise(meanValDat = mean(ValDat, na.rm=TRUE), meanMinDat = mean(MinDat, na.rm=TRUE), meanMaxDat = mean(MaxDat, na.rm=TRUE), meanMinMaxDat = mean(c(MinDat, MaxDat), na.rm=TRUE) )
#view(summ_q_df)

###
# Note: there are alternative, perhaps better, ways to summarise the quantitative trait
# values per species, particularly when min/max data are available.
###

###
# Exploring mean of means 
###

#make column for mean of means
summ_q_df$meanofmeans<-(summ_q_df$meanMinDat + summ_q_df$meanMaxDat)/2

#column of true or false whether values match
summ_q_df$true_false <- summ_q_df$meanofmeans == summ_q_df$meanMinMaxDat

#what is the relation between mean of all min/max values, and mean of means
plot(log(summ_q_df$meanMinMaxDat),log(summ_q_df$meanofmeans))

#relation between mean Val and mean of means
plot(log(summ_q_df$meanValDat),log(summ_q_df$meanofmeans))

#df where values are not equal (NOTE: they appear to be equal)
summ_q_df[grep("TRUE",!summ_q_df$meanMinMaxDat==summ_q_df$meanofmeans),]

#only data with min values
q_df_min<-q_df[!is.na(q_df$MinDat),]

#counts per species per trait
df2 <- q_df_min[,c(2:3)] %>% group_by(NTaxDat , NChrDat) %>% 
  summarise(total_count=n(),.groups = 'drop') %>% 
  as.data.frame()

#view(q_df_min)
#view(df2)

#replace missing VAL with values calculated from min/max
for(i in 1:length(summ_q_df$meanValDat)){
  
  if(is.na(summ_q_df$meanValDat[i])){
    
    # option #1
    # if val is NA, replace with mean of min and max
    # Note: if one of min/max is NA then this may produce an NA 
    # Note 2: this should be the same as method below, but leaving this in just in case
    # (this doesnt occur in latest dataset after cursory check 07/07/23)
    
    #summ_q_df$meanValDat[i] <- (summ_q_df$meanMinDat[i] + summ_q_df$meanMaxDat[i]) / 2
    
    #option #2
    # if val is NA replace with mean of all min and max values
    
    summ_q_df$meanValDat[i] <- summ_q_df$meanMinMaxDat[i]
    
  }
}

#examine data
#view(summ_q_df)

#are any NA present?
table(is.na(summ_q_df$meanValDat))

#pivot to wide format with traits as columns
df.wide <- pivot_wider(summ_q_df, names_from = NChrDat, values_from = c(meanValDat,meanMinDat,meanMaxDat),id_cols = c(NTaxDat))
#view(df.wide)

#write output csv
write.csv(df.wide,"outputs/2_proteus_quantitative.csv")

