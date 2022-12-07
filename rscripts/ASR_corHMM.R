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

#empty vectors
states<-vector()
rates<-vector()
tree_sizes<-vector()
no_states<-vector()

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

#make woodiness binary SOLVE LATER BY MAKING NEW DF FOR ASR
#df2$Woodiness<-gsub("herbaceous_woody","woody",df2$Woodiness)

#plot phylogeny and example trait
plot(phy2, show.tip.label = FALSE,type='fan')
tiplabels(pch = 16, col = as.factor(df3[,2]), cex = 1)

#fit model with 1 rate category on woodiness trait
MK_2state <- corHMM(phy = phy2, data = df3, rate.cat = 1,model="ER")

states[i-1]<-colnames(df3)[2]
rates[i-1]<-MK_2state$solution[1,2]
tree_sizes[i-1]<-length(phy2$tip.label)
no_states[i-1]<-length(unique(df3[,2]))
}

df_rates<-data.frame(states,rates,tree_sizes,no_states)
df_rates[order(df_rates$rates,decreasing = T),]

#EXTRA FROM TUTORIAL:
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

