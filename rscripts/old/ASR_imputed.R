rm(list = ls())

#load libraries
library(dplyr)
library(ape)
library(corHMM)
library(RColorBrewer)
library(phytools)
library(plotrix)

#Load formatted data
#NOTE: Uncomment depending on whether you want standard or one-hot
#df <- read.csv(file = here::here("outputs/one_hot_imputed_with_phylo.csv"))
df <- read.csv(file = here::here("outputs/imputed_with_phylo.csv"),stringsAsFactors = TRUE)

#numeric/factor columns only
nums <- unlist(lapply(df, is.numeric))
facts <- unlist(lapply(df, is.factor))

#make only categorical dataset
df2 <- df[, facts]

#insert '_' into rownames to match phylo
colnames(df2)[1] <- "species"

#remove rownames
rownames(df2) <- NULL

#read in phylogenetic tree
phy <- read.tree("outputs/pruned_tree.tre")
plot(phy, cex = 0.5)

#in dataset but not in phylo
setdiff(df2$species, phy$tip.label)

#in phylo but not in dataset
setdiff(phy$tip.label, df2$species)

#empty vectors
traits <- vector()
rates <- list()
mods <- list()
tree_sizes <- vector()
no_states <- vector()

#make states polymorphic
for(i in 2:length(colnames(df2))){
  df2[,i] <- gsub("_","&",df2[,i])
}

#UNCOMMENT:
#load previous run
mods<-readRDS(file = here::here("outputs/imputed_transition_rate_models.rds"))

# #UNCOMMENT:
# #loop through all characters
# for (i in 2:length(colnames(df2))) {
# 
#  #make new dataframe with only trait of interest
#  df3 <- df2[, c(1, i)]
#  
#  #omit missing data
#  df3 <- na.omit(df3)
#  
#  #drop tips not in dataset
#  phy2 <- drop.tip(phy, setdiff(phy$tip.label, df3$species))
#  
#  #sort df to match order of tips in phylo
#  df3 <- df3[match(phy2$tip.label, df3$species), ]
#  
#  #plot phylogeny and example trait
#  plot(phy2, show.tip.label = FALSE, type = 'fan')
#  tiplabels(pch = 16,
#            col = as.factor(df3[, 2]),
#            cex = 1)
#  
#  #fit model with 1 rate category
#  mod <- corHMM(
#    phy = phy2,
#    data = df3,
#    rate.cat = 1,
#    model = "ARD"
#  )
#  
#  mods[[i - 1]] <- mod
#  traits[i - 1] <- colnames(df3)[2]
#  rates[[i - 1]] <- mod$solution
#  tree_sizes[i - 1] <- length(phy2$tip.label)
#  no_states[i - 1] <- length(unique(df3[, 2]))
#  
# }
# 
# names(mods)<-traits
# names(rates)<-traits
# 
# #NOTE: uncomment to save models transition rate matrices
# saveRDS(mods, file = here::here("outputs/imputed_transition_rate_models.rds"))

#list of models
mods

sg_mods<-readRDS(file = here::here("outputs/mk_list.rds"))


sg_mods$woodiness$mat

mods$Woodiness$states
