
#cleanup
disc_df<-read.csv("outputs/proteus_discrete_recoded.csv")
disc_df<-disc_df[,-1]
str(disc_df)

qr_df<-read.csv("outputs/proteus_quant_recoded.csv")
qr_df<-qr_df[,-1]

quant_df<-read.csv("outputs/proteus_quantitative.csv")
quant_df<-quant_df[,-1]

#merge discrete and recoded discrete
disc_qr_df<-merge(disc_df, qr_df, by.x = 'NTaxDat', by.y = 'species', all.x = T)
disc_qr_df[,c('NTaxDat','Mating','outcrossing_rate')]

#fill in gaps with new data

for(i in 1:length(disc_qr_df$NTaxDat)){
  
  #if na value for Mating column
  if(!is.na(disc_qr_df$Mating[i])){
    
    #add discretized outcrossing rate
    disc_qr_df$Mating[i]<-disc_qr_df$outcrossing_rate[i]
    
  }
  
}

#make all outcrossing/mixed/selfing combinations into mixed
disc_qr_df$Mating[grep("_",disc_qr_df$Mating)]<-'mixed'

#remove quant outcrossing
disc_qr_df<-disc_qr_df[,-grep("outcrossing_rate",colnames(disc_qr_df))]
