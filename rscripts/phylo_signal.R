rm(list=ls())
# FROM: https://github.com/mrborges23/delta_statistic

#load libraries
library(ape)
library(phylosignal)
library(adephylo)
library(phylobase)

#load functions
source("R/delta.R")

#set margins
par(mar=c(3,3,3,3))

#NOTE: choose based on whether you want standard or one-hot data
#load formatted data
#df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))
df<-readRDS(file = here::here("outputs/one_hot_6_df_filt_trans.rds"))

#numeric/factor columns only
facts <- unlist(lapply(df, is.factor))

#make only categorical dataset
df2<-df[ , facts]

#insert '_' into rownames to match phylo
df2$species<-gsub(" ","_",rownames(df2))

#remove rownames
rownames(df2)<-NULL

#put species column first
df2<-df2[,c(length(colnames(df2)),1:(length(colnames(df2)) - 1))]

#NOTE: for standard data set only
#remove FloralReward, which currently has too many categories
#df2<-subset(df2, select = -c(FloralReward))

#read in phylogenetic tree
phy<-read.tree("outputs/pruned_tree.tre")
plot(phy,cex=0.5)

#in dataset but not in phylo
setdiff(df2$species,phy$tip.label)

#in phylo but not in dataset
setdiff(phy$tip.label,df2$species)

#empty vectors
deltas<-vector()
pvals<-vector()

pdf("figures/asr_phylos.pdf")

#NOTE: Produces NaNs
#loop through all characters
#for(i in 2:length(colnames(df2))){
for(i in 2:length(colnames(df2))){
  #make new dataframe with only trait of interest
  df3<-df2[,c(1,i)]
  
  #omit missing data
  df3<-na.omit(df3)
  
  #drop tips not in dataset
  phy2<-drop.tip(phy,setdiff(phy$tip.label,df3$species))
  
  #sort df to match order of tips in phylo
  df3<-df3[match(phy2$tip.label, df3$species),]
  
  phy2$tip.label==df3$species
  
  trait <- df3[,2]
  tree <- phy2
  
  #plot phylogeny and example trait
  plot(tree, show.tip.label = FALSE,type='fan')
  tiplabels(pch = 16, col = as.factor(trait), cex = 1)
  
  #calculate delta
  deltaA <- delta(trait,tree,0.1,0.0589,10000,10,100)
  
  random_delta <- rep(NA,100)
  for (j in 1:100){
    rtrait <- sample(trait)
    random_delta[j] <- delta(rtrait,tree,0.1,0.0589,10000,10,100)
  }
  
  p_value <- sum(random_delta>deltaA)/length(random_delta)
  boxplot(random_delta,ylim=c(0,(max(c(random_delta,deltaA))+1)),main=(paste(colnames(df2)[i])))
  abline(h=deltaA,col="red")
  
  pvals[i]<-p_value
  deltas[i]<-deltaA
  
  cat(paste("finished: ",colnames(df2)[i],"\n"))
  
}

dev.off()

#clean up
names(deltas)<-colnames(df2)
deltas<-na.omit(deltas)

names(pvals)<-colnames(df2)
pvals<-na.omit(pvals)

results<-data.frame(deltas,pvals)
results

#NOTE: UNCOMMENT TO SAVE
saveRDS(results, file = here::here("outputs/phylo_signal_qualitative.rds"))

###
# ---- Continuous ----
### 
# FROM: http://www.francoiskeck.fr/phylosignal/demo_general.html

#numeric/factor columns only
nums <- unlist(lapply(df, is.numeric))

#make only continuous dataset
df2<-df[ , nums]

#insert '_' into rownames to match phylo
df2$species<-gsub(" ","_",rownames(df2))

#remove rownames
rownames(df2)<-NULL

#put species column first
df2<-df2[,c(length(colnames(df2)),1:(length(colnames(df2)) - 1))]

#read in phylogenetic tree
phy<-read.tree("outputs/pruned_tree.tre")

#in dataset but not in phylo
setdiff(df2$species,phy$tip.label)

#in phylo but not in dataset
setdiff(phy$tip.label,df2$species)

#loop through all characters
for(i in 2:length(colnames(df2))){
  #make new dataframe with only trait of interest
  df3<-df2[,c(1,i)]
  
  #omit missing data
  df3<-na.omit(df3)
  
  #drop tips not in dataset
  phy2<-drop.tip(phy,setdiff(phy$tip.label,df3$species))
  
  #sort df to match order of tips in phylo
  df3<-df3[match(phy2$tip.label, df3$species),]
  
  phy2$tip.label==df3$species
  
  trait <- df3[,2]
  tree <- phy2
  
  #make p4d object
  p4d <- phylo4d(tree, trait)
  
  #plot phylogeny and trait
  barplot.phylo4d(p4d, tree.type = "phylo", tree.ladderize = TRUE)
  
  #calculate signal with all methods
  p4d_s<-phyloSignal(p4d = p4d, method = "all")
  
  if(i == 2){
    pvals <- p4d_s$pvalue
    signals <- p4d_s$stat
  } else {
    pvals <- rbind(pvals,p4d_s$pvalue)
    signals <- rbind(signals,p4d_s$pvalue)
  }
  
  cat(paste("finished: ",colnames(df2)[i],"\n"))
  
}

#row and column names
rownames(pvals)<-colnames(df2)[2:length(colnames(df2))]
rownames(signals)<-colnames(df2)[2:length(colnames(df2))]

#signal values
signals

#p values
colnames(pvals)<-paste("p_",colnames(pvals),sep="")

#NOTE: UNCOMMENT TO SAVE
saveRDS(cbind(pvals,signals), file = here::here("outputs/phylo_signal_quantitative.rds"))

