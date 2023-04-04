
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
  if(is.na(disc_qr_df$Mating[i]) & !is.na(disc_qr_df$outcrossing_rate[i])){
    
    #add discretized outcrossing rate
    disc_qr_df$Mating[i]<-disc_qr_df$outcrossing_rate[i]
    
  }
  
}

#make all outcrossing/mixed/selfing combinations into mixed
disc_qr_df$Mating[grep("_",disc_qr_df$Mating)]<-'mixed'

#remove quant outcrossing
disc_qr_df<-disc_qr_df[,-grep("outcrossing_rate",colnames(disc_qr_df))]

#view(disc_qr_df)

#remove min and max
quant_df<-quant_df[,c(1,grep("meanValDat",colnames(quant_df)))]

#remove quant outcrossing
quant_df<-quant_df[,-grep("Outcrossing",colnames(quant_df))]

#merge discrete+recoded discrete and quantitative
disc_qr_quant_df<-merge(disc_qr_df, quant_df, by.x = 'NTaxDat', by.y = 'NTaxDat', all.x = T)

#view(disc_qr_quant_df)

#clean up names
colnames(disc_qr_quant_df)<-gsub("\\.","",colnames(disc_qr_quant_df))
colnames(disc_qr_quant_df)<-gsub("C1","",colnames(disc_qr_quant_df))
colnames(disc_qr_quant_df)<-gsub("meanValDat_","",colnames(disc_qr_quant_df))
colnames(disc_qr_quant_df)<-gsub('[0-9]+', '', colnames(disc_qr_quant_df))
colnames(disc_qr_quant_df)

#remove first row (numbers) and add to rownames
proteus_combined<-disc_qr_quant_df[,-1]
rownames(proteus_combined)<-disc_qr_quant_df[,"NTaxDat"]

#combine length and diameter
fl<-proteus_combined$Flowerlength
fw<-proteus_combined$Flowerdiameter

#make NA 0 to get max values
fl[is.na(fl)]<-0
fw[is.na(fw)]<-0

#retrieve max values
fs<-pmax(fl,fw)
fs

#convert 0s back to NAs
fs[fs==0]<-NA

#make new column with maximum of length and diameter
proteus_combined$flowerSize <- fs

#remove flowerDiameter (included in flowerSize)
proteus_combined<-subset(proteus_combined, select=-c(Flowerdiameter))

#fix column name
colnames(proteus_combined)[1]<-"Woodiness"


write.csv(proteus_combined,"outputs/proteus_combined.csv")
