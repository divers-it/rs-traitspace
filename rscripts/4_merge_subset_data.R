rm(list=ls())

# load discrete data
disc_df<-read.csv("outputs/1_proteus_discrete_recoded.csv",row.names = 1)
head(disc_df)

# load quantitative data
quant_df<-read.csv("outputs/2_proteus_quantitative.csv",row.names = 1)
str(quant_df)

# load recoded quantitative data
qr_df<-read.csv("outputs/3_proteus_quant_recoded.csv",row.names = 1)
head(qr_df)

# merge discrete and recoded discrete
disc_qr_df<-merge(disc_df, qr_df, by.x = 'NTaxDat', by.y = 'species', all.x = T)
head(disc_qr_df[,c('NTaxDat','Mating','outcrossing_rate')])

# fill in gaps with recoded outcrossing rate data
for(i in 1:length(disc_qr_df$NTaxDat)){
  
  # if na value for Mating column
  if(is.na(disc_qr_df$Mating[i]) && !is.na(disc_qr_df$outcrossing_rate[i])){
    
    # add discretized outcrossing rate
    disc_qr_df$Mating[i]<-disc_qr_df$outcrossing_rate[i]
    
    # if there are states for both sources of outcrossing rate
  } else if (!is.na(disc_qr_df$Mating[i]) && !is.na(disc_qr_df$outcrossing_rate[i])){
    
    # if the states are not the name
    if(disc_qr_df$Mating[i]!=disc_qr_df$outcrossing_rate[i]){
      
      # replace Mating state with polymorphic state
      disc_qr_df$Mating[i]<-paste(disc_qr_df$Mating[i],disc_qr_df$outcrossing_rate[i],sep="_")
      
    }
    
  }
  
}

# make all outcrossing/mixed/selfing combinations into mixed
disc_qr_df$Mating[grep("_",disc_qr_df$Mating)]<-'mixed'

# remove recoded outcrossing
disc_qr_df<-disc_qr_df[,-grep("outcrossing_rate",colnames(disc_qr_df))]
# view(disc_qr_df)

# remove min and max from quantitative df
quant_df<-quant_df[,c(1,grep("meanValDat",colnames(quant_df)))]

# remove quant outcrossing
quant_df<-quant_df[,-grep("Outcrossing",colnames(quant_df))]

# merge discrete+recoded discrete and quantitative
disc_qr_quant_df<-merge(disc_qr_df, quant_df, by.x = 'NTaxDat', by.y = 'NTaxDat', all.x = T)

# view(disc_qr_quant_df)

# clean up names
colnames(disc_qr_quant_df)<-gsub("\\.","",colnames(disc_qr_quant_df))
colnames(disc_qr_quant_df)<-gsub("C1","",colnames(disc_qr_quant_df))
colnames(disc_qr_quant_df)<-gsub("meanValDat_","",colnames(disc_qr_quant_df))
colnames(disc_qr_quant_df)<-gsub('[0-9]+', '', colnames(disc_qr_quant_df))
colnames(disc_qr_quant_df)

# remove first row and add to rownames
proteus_combined<-disc_qr_quant_df[,-1]
rownames(proteus_combined)<-disc_qr_quant_df[,"NTaxDat"]

# combine length and diameter
fl<-proteus_combined$Flowerlength
fw<-proteus_combined$Flowerdiameter

# retrieve max values
fs<-pmax(fl,fw,na.rm = TRUE)

# make new column with maximum of length and diameter
proteus_combined$flowerSize <- fs

# remove flowerDiameter and Length (included in flowerSize)
proteus_combined<-subset(proteus_combined, select=-c(Flowerdiameter))
proteus_combined<-subset(proteus_combined, select=-c(Flowerlength))

####
## Correct species names ----
####

# get species list
spec_list <- rownames(proteus_combined)

Sys.sleep(1)
# run TNRS to check species (best result only)
check_species <- TNRS(spec_list, matches="best", sources="wcvp")
Sys.sleep(1)

# how many name issues?
table(check_species$Name_submitted == check_species$Accepted_species)

# species with issues
issue_species <- check_species[check_species$Name_submitted != check_species$Accepted_species,]

# examine species
issue_species[,c("Name_submitted",
                 "Accepted_name")]
# NOTE: verified on POWO: https://powo.science.kew.org/

# set correct row names
table(rownames(proteus_combined) == check_species$Name_submitted)
rownames(proteus_combined) <- check_species$Accepted_species

# order rows based on new names
proteus_combined<-proteus_combined[order(rownames(proteus_combined)),]

#### 
## Add seed mass data ----
#### 

# manually added Symplocos rhamnifolia row to data set (no data available) as was not present

# Removed one duplicate row Corokia
seedMass<-read.csv("data/seedWeight.csv")
rownames(seedMass) <- paste(seedMass$Genus,seedMass$Species,sep=" ")

####
### Correct seed mass names ----
####

# get species list
spec_list <- rownames(seedMass)

Sys.sleep(1)
# run TNRS to check species (best result only)
check_species <- TNRS(spec_list, matches="best", sources="wcvp")
Sys.sleep(1)

# how many name issues?
table(check_species$Name_submitted == check_species$Accepted_species)

# species with issues
issue_species <- check_species[check_species$Name_submitted != check_species$Accepted_species,]

# examine species
issue_species[,c("Name_submitted",
                 "Accepted_name")]
# NOTE: verified on POWO: https://powo.science.kew.org/

# set correct row names
table(rownames(seedMass) == check_species$Name_submitted)
rownames(seedMass) <- check_species$Accepted_species

# check rownames match
rownames(seedMass)==rownames(proteus_combined)

# check differences between data sets
setdiff(rownames(proteus_combined),rownames(seedMass))
setdiff(rownames(seedMass),rownames(proteus_combined))

# Fix synonymy issues
rownames(seedMass)[grep("Arctostaphylos uvaursi",rownames(seedMass))]<-"Arctostaphylos uva-ursi"
rownames(seedMass)[grep("Cleistes bifaria",rownames(seedMass))]<-"Cleistesiopsis bifaria"
rownames(seedMass)[grep("Hydnocarpus heterophylla",rownames(seedMass))]<-"Hydnocarpus heterophyllus"
rownames(seedMass)[grep("Ruellia nudiflora",rownames(seedMass))]<-"Ruellia ciliatiflora"
rownames(seedMass)[grep("Peritoma arborea",rownames(seedMass))]<-"Cleomella arborea"
rownames(seedMass)[grep("Pitcairnia albifilos",rownames(seedMass))]<-"Pitcairnia albiflos"
rownames(seedMass)[grep("Veronica anagallisaquatica",rownames(seedMass))]<-"Veronica anagallis-aquatica"

str(seedMass[rownames(proteus_combined)%in%rownames(seedMass),])

# reorder seedMass data frame
seedMass <- seedMass[order(rownames(seedMass)),]

# percentage missing data
table(is.na(seedMass$Seed_weight))[2]/length(seedMass$Seed_weight)

# make df for merging
seedMass_merge<-seedMass[rownames(proteus_combined)%in%rownames(seedMass),]

# check row names after fixing synonyms
rownames(proteus_combined)==rownames(seedMass_merge)
cbind(rownames(proteus_combined),rownames(seedMass_merge))


# check differences between data sets
setdiff(rownames(proteus_combined),rownames(seedMass))
setdiff(rownames(seedMass),rownames(proteus_combined))

# merge
proteus_combined<-cbind(proteus_combined,seedMass_merge$Seed_weight)

# rename
colnames(proteus_combined)[ncol(proteus_combined)]<-"seedMass"

# write dataset
write.csv(proteus_combined,"outputs/4_proteus_combined.csv")

## Supplementary data ----

cleaned_data <- proteus_combined

# get family matches
tnrs_out <- TNRS::TNRS(rownames(proteus_combined), matches = "best")

#  add to output
cleaned_data <- cbind(tnrs_out$Accepted_family, cleaned_data)
colnames(cleaned_data)[1] <- "family"

#  write to supplementary folder
write.csv(cleaned_data, "supplementary_data/cleaned_data.csv")
