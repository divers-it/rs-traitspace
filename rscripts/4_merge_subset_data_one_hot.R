
#cleanup
disc_df<-read.csv("outputs/proteus_discrete_one_hot.csv")
disc_df<-disc_df[,-1]
str(disc_df)

qr_df<-read.csv("outputs/proteus_quant_recoded_one_hot.csv")
qr_df<-qr_df[,-1]

quant_df<-read.csv("outputs/proteus_quantitative.csv")
quant_df<-quant_df[,-1]

#remove min and max
quant_df<-quant_df[,c(1,grep("meanValDat",colnames(quant_df)))]

#merge discrete+recoded discrete
disc_qr_df<-merge(disc_df, qr_df, by.x = 'NTaxDat', by.y = 'species', all.x = T)

#merge mating traits (discrete- and rate-based)
disc_qr_df[is.na(disc_qr_df$Mating_selfing),]$Mating_selfing=disc_qr_df[is.na(disc_qr_df$Mating_selfing),]$outcrossing_rate_selfing
disc_qr_df[!is.na(disc_qr_df$outcrossing_rate_selfing) & disc_qr_df$outcrossing_rate_selfing==1,]$Mating_selfing=disc_qr_df[!is.na(disc_qr_df$outcrossing_rate_selfing) & disc_qr_df$outcrossing_rate_selfing==1,]$outcrossing_rate_selfing
disc_qr_df[is.na(disc_qr_df$Mating_outcrossing),]$Mating_outcrossing=disc_qr_df[is.na(disc_qr_df$Mating_outcrossing),]$outcrossing_rate_outcrossing
disc_qr_df[!is.na(disc_qr_df$outcrossing_rate_outcrossing) & disc_qr_df$outcrossing_rate_outcrossing==1,]$Mating_outcrossing=disc_qr_df[!is.na(disc_qr_df$outcrossing_rate_outcrossing) & disc_qr_df$outcrossing_rate_outcrossing==1,]$outcrossing_rate_outcrossing
disc_qr_df[is.na(disc_qr_df$Mating_mixed),]$Mating_mixed=disc_qr_df[is.na(disc_qr_df$Mating_mixed),]$outcrossing_rate_mixed
disc_qr_df[!is.na(disc_qr_df$outcrossing_rate_mixed) & disc_qr_df$outcrossing_rate_mixed==1,]$Mating_mixed=disc_qr_df[!is.na(disc_qr_df$outcrossing_rate_mixed) & disc_qr_df$outcrossing_rate_mixed==1,]$outcrossing_rate_mixed

#remove rate-based mating system
disc_qr_df<-disc_qr_df[,-c(grep("outcrossing_rate",colnames(disc_qr_df)))]

#merge with quantitative traits
df<-merge(disc_qr_df, quant_df, by.x = 'NTaxDat', by.y = 'NTaxDat', all.x = T)

#clean up names
colnames(df)<-gsub("\\.+",".",colnames(df))
colnames(df)<-gsub("C1","",colnames(df))
colnames(df)<-gsub("meanValDat_","",colnames(df))
colnames(df)<-gsub('[0-9]+', '', colnames(df))
#colnames(df)<-gsub("syndrome","",colnames(df))
#colnames(df)<-gsub("Plant.","",colnames(df))
#colnames(df)<-gsub("sexual.system.","",colnames(df))
#colnames(df)<-gsub("Life.history.","",colnames(df))
colnames(df)<-gsub(".structural","",colnames(df))
colnames(df)<-gsub(".functional","",colnames(df))
#colnames(df)<-gsub("Symmetry.of.perianth.","",colnames(df))
#colnames(df)<-gsub("Self.incompatibility.system.","",colnames(df))
#colnames(df)<-gsub("position.","",colnames(df))
#colnames(df)<-gsub("Main.color.of.perianth.at.anthesis.","Flowercolor.",colnames(df))
colnames(df)<-gsub("\\.+",".",colnames(df))
colnames(df)<-gsub("^\\.","",colnames(df))
colnames(df)<-gsub("\\.$","",colnames(df))
colnames(df)

#remove ovary position
df=df[,-grep("OvaryPosition",colnames(df))]
#remove first column (names) and add to rownames
proteus_combined<-df[,-1]
rownames(proteus_combined)<-df[,"NTaxDat"]

#combine length and diameter
fl<-proteus_combined$Flower.length
fw<-proteus_combined$Flower.diameter

#make NA 0 to get max values
fl[is.na(fl)]<-0
fw[is.na(fw)]<-0

#retrieve max values
fs<-pmax(fl,fw)
fs

#convert 0s back to NAs
fs[fs==0]<-NA

#make new column with maximum of length and diameter
proteus_combined$Flower.size <- fs

#remove flower diameter and length (included in flowerSize)
proteus_combined<-subset(proteus_combined, select=-c(Flower.diameter))
proteus_combined<-subset(proteus_combined, select=-c(Flower.length))


#####
#add seed mass data
#####

#remove one duplicate row Corokia
#synonyms are present e.g. Cleistes -> Cleistesiopsis
seedMass<-read.csv("data/seedWeight.csv")
rownames(seedMass)<-paste(seedMass$Genus,seedMass$Species,sep=" ")
rownames(seedMass)==rownames(proteus_combined)

#check differences between datasets
setdiff(rownames(proteus_combined),rownames(seedMass))
setdiff(rownames(seedMass),rownames(proteus_combined))

#Synonymy issues
rownames(seedMass)[grep("Arctostaphylos uvaursi",rownames(seedMass))]<-"Arctostaphylos uva-ursi"
rownames(seedMass)[grep("Cleistes bifaria",rownames(seedMass))]<-"Cleistesiopsis bifaria"
rownames(seedMass)[grep("Pitcairnia albifilos",rownames(seedMass))]<-"Pitcairnia albiflos"
rownames(seedMass)[grep("Ruellia nudiflora",rownames(seedMass))]<-"Ruellia ciliatiflora"
rownames(seedMass)[grep("Veronica anagallisaquatica",rownames(seedMass))]<-"Veronica anagallis-aquatica"

str(seedMass[rownames(proteus_combined)%in%rownames(seedMass),])

#percentage missing data
table(is.na(seedMass$Seed_weight))[2]/length(seedMass$Seed_weight)

#make df for merging
seedMass_merge<-seedMass[rownames(proteus_combined)%in%rownames(seedMass),]

#check row names
rownames(proteus_combined)==rownames(seedMass_merge)

#merge
proteus_combined<-cbind(proteus_combined,seedMass_merge$Seed_weight)

write.csv(proteus_combined,"outputs/proteus_combined_one_hot.csv")
