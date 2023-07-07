rm(list=ls())

#load libraries
library(reshape2)

#read in quantitative data
df<-read.csv("outputs/proteus_quantitative.csv")

#get columns for val, min and max outcrossing rate
outc<-df[,grep("84.",colnames(df))]

#add rownames
rownames(outc)<-df$NTaxDat

####
# Notes: columns need to be in order VAL, MIN, MAX
####

#examine raw outcrossing rate data
str(outc)
head(outc)
#view(outc)

#empty vector for output
outc_disc<-vector()

#loop through species, assigning state based on outcrossing rates
for(i in 1:length(outc[,1])){
  
  #if all three values are NA assign NA for species
  if(is.na(outc[i,1]) && is.na(outc[i,2]) && is.na(outc[i,3])){
    
    outc_disc[i]<-NA
    
    #if both min and max have values
  } else if(!is.na(outc[i,2]) && !is.na(outc[i,3])) {
    
    #make all species mixed and those where values are only high or only low outcrossing / selfing
    outc_disc[i] <- "mixed"
    if(outc[i,2] < 0.2 && outc[i,3] < 0.2){outc_disc[i] <- "selfing"}
    if(outc[i,2] > 0.8 && outc[i,3] > 0.8){outc_disc[i] <- "outcrossing"}
    
  } else {
    #for all other cases use the VAL
    #NOTE: Carex has a min value but not a max value, but this does not change classification as both < 0.2
    if(outc[i,1] < 0.2){outc_disc[i] <- "selfing"}
    if(outc[i,1] > 0.2 && outc[i,1] < 0.8){outc_disc[i] <- "mixed"}
    if(outc[i,1] > 0.8){outc_disc[i] <- "outcrossing"}
    
  }
  
}

#make data frame and add rownames and column names
outc_df<-data.frame(rownames(outc),outc_disc)
colnames(outc_df)<-c("species","outcrossing_rate")

#get rid of NA to prepare for merging with discrete data
outc_df <- na.omit(outc_df)
head(outc_df)

#cast to make one column per outcrossing rate state (one-hot)
outc_df_casted <- dcast(data=melt(outc_df,id.vars="species"), species ~ variable + value, length)
head(outc_df_casted)

#write output csv
write.csv(outc_df,"outputs/proteus_quant_recoded_one_hot.csv")
