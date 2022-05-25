library(dplyr)
library(ape)
library(phytools)
#load formatted data
df<-readRDS(file = here::here("outputs/df_filt.rds"))

#numeric/factor columns only
nums <- unlist(lapply(df, is.numeric))
facts <- unlist(lapply(df, is.factor))

#make only categorical dataset
df2<-df[ , facts]

#insert '_' into rownames to match phylo
df2$species<-gsub(" ","_",rownames(df2))

#remove rownames
#rownames(df2)<-NULL

#put species column first
df2<-df2[,c(length(colnames(df2)),1:(length(colnames(df2)) - 1))]

#read in phylogenetic tree
phy<-read.tree("outputs/pruned_tree.tre")
plot(phy,cex=0.5)

#in dataset but not in phylo
setdiff(df2$species,phy$tip.label)

#in phylo but not in dataset
setdiff(phy$tip.label,df2$species)

#drop tips not in dataset
phy2<-drop.tip(phy,setdiff(phy$tip.label,df2$species))

#sort df to match order of tips in phylo
df2<-df2[match(phy2$tip.label, df2$species),]

#make traits binary SOLVE LATER BY MAKING NEW DF FOR ASR / POLYMORPH
df2$Woodiness<-gsub("herbaceous_woody","woody",df2$Woodiness)
df2$Climbing<-gsub("climbing_non-climbing","climbing",df2$Climbing)


#plot phylogeny and example trait
plot(phy2, show.tip.label = FALSE,type='fan')
tiplabels(pch = 16, col = as.factor(df2$Woodiness), cex = 1)

wood<-df2$Woodiness
names(wood)<-phy2$tip.label

#fits an extended Mk model for discrete character evolution (Lewis, 2001).
fit.ER.wood<-fitMk(phy2, wood, model="ER", fixedQ=NULL)
fit.ER.wood

climb<-df2$Climbing
names(climb)<-phy2$tip.label

#fits an extended Mk model
fit.ER.climb<-fitMk(phy2, climb, model="ER", fixedQ=NULL)
fit.ER.climb

#MCMC version BROKEN
#mcmcMk(phy2, climb, model="ER", ngen=100)

