#load library
library(tidyverse)
library(gtools)

#import data
df<-read.csv("data/qryDiveRS_Data_2022-04-29.csv")
str(df)

#filter by quality
#TO DO

#extract quant data only
q_df<-df[grep("(C1)",df$NChrDat),]

#remove extra metadata columns
q_df<-q_df[,c(1,2,3,6:8)]

#calculate mean value (VAL,MIN,MAX) for each trait for each species
summ_q_df<-q_df %>% group_by(NTaxDat,NChrDat) %>% summarise(meanValDat = mean(ValDat, na.rm=TRUE), meanMinDat = mean(MinDat, na.rm=TRUE), meanMaxDat = mean(MaxDat, na.rm=TRUE))

#view(summ_q_df)

#if val is NA, replace with mean of min and max
for(i in 1:length(summ_q_df$meanValDat)){
  
  if(is.na(summ_q_df$meanValDat[i])){
    
    summ_q_df$meanValDat[i] <- (summ_q_df$meanMinDat[i] + summ_q_df$meanMaxDat[i]) / 2
    
  }
}

view(summ_q_df)


#pivot to wide format with traits as columns
df.wide <- pivot_wider(summ_q_df, names_from = NChrDat, values_from = c(meanValDat,meanMinDat,meanMaxDat),id_cols = c(NTaxDat))
view(df.wide)

write.csv(df.wide,"outputs/proteus_quantitative.csv")

