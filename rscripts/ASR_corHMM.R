library(dplyr)
library(ape)
library(corHMM)
#load formatted data
df<-readRDS(file = here::here("outputs/df_filt_trans.rds"))

#numeric/factor columns only
nums <- unlist(lapply(df, is.numeric))
facts <- unlist(lapply(df, is.factor))

#make only categorical dataset
df2<-df[ , facts]

#insert '_' into rownames to match phylo
df2$species<-gsub(" ","_",rownames(df2))

#remove rownames
rownames(df2)<-NULL

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

#make woodiness binary SOLVE LATER BY MAKING NEW DF FOR ASR
df2$Woodiness<-gsub("herbaceous_woody","woody",df2$Woodiness)

#plot phylogeny and example trait
plot(phy2, show.tip.label = FALSE,type='fan')
tiplabels(pch = 16, col = as.factor(df2$Woodiness), cex = 1)

#fit model with 1 rate category on woodiness trait
MK_2state <- corHMM(phy = phy2, data = df2[,c(1:2)], rate.cat = 1)
MK_2state

#plot transition rates
plotMKmodel(MK_2state)

#retrieve data, model and tree from corHMM object
phy <- MK_2state$phy
data <- MK_2state$data
model <- MK_2state$solution

#not sure what this is doing - making it symmetric?
model[is.na(model)] <- 0
diag(model) <- -rowSums(model)

# run stochastic character mapping
simmap <- makeSimmap(tree = phy, data = data, model = model, rate.cat = 1, nSim = 1,
                     nCores = 1)

#plot results
phytools::plotSimmap(simmap[[1]], fsize = 0.5, type="fan")

#add another rate category
HMM_2state <- corHMM(phy = phy, data = data, rate.cat = 2, model = "SYM", get.tip.states = TRUE)
HMM_2state

