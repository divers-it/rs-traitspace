df<-read.csv("outputs/proteus_quantitative.csv")

outc<-df[,grep("84.",colnames(df))]

rownames(outc)<-df$NTaxDat

####
# needs to be in order VAL, MIN, MAX
####
str(outc)
head(outc)
#view(outc)

outc_disc<-vector()


for(i in 1:length(outc[,1])){
  
  #if all three values are NA assign NA for species
  if(is.na(outc[i,1]) && is.na(outc[i,2]) && is.na(outc[i,3])){
    
    outc_disc[i]<-NA
    
    #if both min and max have values
  } else if(!is.na(outc[i,2]) && !is.na(outc[i,3])) {
    
    outc_disc[i] <- "mixed"
    if(outc[i,2] < 0.2 && outc[i,3] < 0.2){outc_disc[i] <- "selfing"}
    if(outc[i,2] > 0.8 && outc[i,3] > 0.8){outc_disc[i] <- "outcrossing"}
    
  } else {
    
    if(outc[i,1] < 0.2){outc_disc[i] <- "selfing"}
    if(outc[i,1] > 0.2 && outc[i,1] < 0.8){outc_disc[i] <- "mixed"}
    if(outc[i,1] > 0.8){outc_disc[i] <- "outcrossing"}
    
  }
  
}




outc_df<-data.frame(rownames(outc),outc_disc)
colnames(outc_df)<-c("species","outcrossing_rate")
head(outc_df)

write.csv(outc_df,"outputs/proteus_quant_recoded.csv")
