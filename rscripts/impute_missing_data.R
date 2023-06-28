library(PVR)
library(ape)

library(dplyr)
library(ape)
library(corHMM)

#load formatted data
df<-readRDS(file = here::here("outputs/df_filt_trans.rds"))
#df<-readRDS(file = here::here("outputs/df_filt_trans_one_hot.rds"))

#insert '_' into rownames to match phylo
df$species<-gsub(" ","_",rownames(df))

#remove rownames
rownames(df)<-NULL

#put species column first
df<-df[,c(length(colnames(df)),1:(length(colnames(df)) - 1))]

#Proportion of missing data
table(is.na(df))[2] / (table(is.na(df))[1] + table(is.na(df))[2])

#read in phylogenetic tree
phy<-read.tree("outputs/pruned_tree.tre")
plot(phy,cex=0.5)

#in dataset but not in phylo
setdiff(df$species,phy$tip.label)

#in phylo but not in dataset
setdiff(phy$tip.label,df$species)

#drop tips not in dataset
phy<-drop.tip(phy,setdiff(phy$tip.label,df$species))

#Decomposing phylogenetic distance matrix derived from tree into a set of orthogonal vectors
x <- PVRdecomp(phy)

str(x)

#first ten eigenvectors
pvrs<-x@Eigen$vector[,c(1:10)]

#IS THIS % EXPLAINED?
sum(x@Eigen$values[1:10])/sum(x@Eigen$values)

#Convert all character columns to factor for missForest
df <- as.data.frame(unclass(df),                     
                       stringsAsFactors = TRUE)

# Imputation without Phylo data
imp <- missForest::missForest(df[2:length(colnames(df))], maxiter = 15, ntree = 100, variablewise = FALSE)

# Combine traits and PVRs
traits.pvrs <- cbind(df[2:length(colnames(df))], pvrs)

# Imputation with Phylo data
phy_imp <- missForest::missForest(traits.pvrs, maxiter = 15, ntree = 100, variablewise = FALSE)

#error
phy_imp$OOBerror

#output dataset
phy_imp_df <- phy_imp$ximp[,c(1:(length(colnames(df))-1))]

#add species names
rownames(phy_imp_df)<-df$species

#CHECK THE ORDER IS CORRECT BETWEEN DF AND IMPUTED DF

#write csv
write.csv(phy_imp_df,"outputs/imputed_with_phylo.csv")

