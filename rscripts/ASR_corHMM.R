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
#df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))
df <- readRDS(file = here::here("outputs/one_hot_6_df_filt_trans.rds"))

#numeric/factor columns only
nums <- unlist(lapply(df, is.numeric))
facts <- unlist(lapply(df, is.factor))

#make only categorical dataset
df2 <- df[, facts]

#insert '_' into rownames to match phylo
df2$species <- gsub(" ", "_", rownames(df2))

#remove rownames
rownames(df2) <- NULL

#put species column first
df2 <- df2[, c(length(colnames(df2)), 1:(length(colnames(df2)) - 1))]

#read in phylogenetic tree
phy <- read.tree("outputs/pruned_tree.tre")
plot(phy, cex = 0.5)

#in dataset but not in phylo
setdiff(df2$species, phy$tip.label)

#in phylo but not in dataset
setdiff(phy$tip.label, df2$species)

#empty vectors
states <- vector()
rates <- vector()
tree_sizes <- vector()
no_states <- vector()

#NOTE: Uncommen to read in transition rates from previous run
#df_rates<-readRDS(file = here::here("outputs/one_hot_transition_rates.rds"))

#loop through all characters
for (i in 2:length(colnames(df2))) {
  #make new dataframe with only trait of interest
  df3 <- df2[, c(1, i)]
  
  #omit missing data
  df3 <- na.omit(df3)
  
  #drop tips not in dataset
  phy2 <- drop.tip(phy, setdiff(phy$tip.label, df3$species))
  
  #sort df to match order of tips in phylo
  df3 <- df3[match(phy2$tip.label, df3$species), ]
  
  #plot phylogeny and example trait
  plot(phy2, show.tip.label = FALSE, type = 'fan')
  tiplabels(pch = 16,
            col = as.factor(df3[, 2]),
            cex = 1)
  
  #fit model with 1 rate category on woodiness trait
  MK_2state <- corHMM(
    phy = phy2,
    data = df3,
    rate.cat = 1,
    model = "ER"
  )
  
  states[i - 1] <- colnames(df3)[2]
  rates[i - 1] <- MK_2state$solution[1, 2]
  tree_sizes[i - 1] <- length(phy2$tip.label)
  no_states[i - 1] <- length(unique(df3[, 2]))
}

df_rates <- data.frame(states, rates, tree_sizes, no_states)
df_rates[order(df_rates$rates, decreasing = T), ]

#NOTE: uncomment to save transition rates
saveRDS(df_rates, file = here::here("outputs/one_hot_transition_rates.rds"))

###
# ---- Compare phylogenetic signal and transition rates ----
###

#read in phylogenetic signal results for qualitiative traits
psq<-readRDS(file = here::here("outputs/phylo_signal_qualitative.rds"))

#reorder
d1 <- psq[order(row.names(psq)), ]

#reorder transition rates
d2 <- df_rates[order(df_rates$states), ]

#check order matches
d2$states == row.names(d1)

#combine data frames
d1d2<-cbind(d1,d2)

#remove outliers
d1d2<-d1d2[d1d2$deltas<100,]

#plot signal vs transition rates
plot(d1d2$deltas ~ d1d2$rates)

###
# ---- Transitions between strategies ----
###

#read in clustering results
clust_df <-
  read.csv("outputs/collate_clustering_results.csv", row.names = 1)

#make new species column and insert '_' into rownames to match phylo
clust_df$species <- gsub(" ", "_", rownames(clust_df))

#remove rownames
rownames(clust_df) <- NULL

#put species column first
clust_df <-
  clust_df[, c(length(colnames(clust_df)), 1:(length(colnames(clust_df)) - 1))]

#read in phylogenetic tree
phy <- read.tree("outputs/pruned_tree.tre")

#in dataset but not in phylo
setdiff(clust_df$species, phy$tip.label)

#in phylo but not in dataset
setdiff(phy$tip.label, clust_df$species)

#drop tips
phy <- drop.tip(phy, setdiff(phy$tip.label, clust_df$species))

#check order
clust_df <- clust_df[match(phy$tip.label, clust_df$species), ]
phy$tip.label == clust_df$species

#NOTE: change column to choose different classification
#make associated data frame by choosing clusters/robust column
dat <- clust_df[, c("species", "ward")]

#add another rate category
HMM_ARD <-
  corHMM(
    phy = phy,
    data = dat,
    rate.cat = 1,
    model = "ARD",
    node.states = 'marginal',
    get.tip.states = TRUE
  )

#examine model
HMM_ARD

#plot transition rates
plotMKmodel(HMM_ARD)

#NOTE: linked to chosen clustering method
#make transition rates table
trans_mat <- data.frame(round(HMM_ARD$solution, 4))
colnames(trans_mat)<-c("->1","->2","->3")
rownames(trans_mat)<-c("1->","2->","3->")

#plot tree with marginal probabilites on nodes
png(
  "figures/phylo_asr_clusters.png",
  width = 1500,
  height = 1500,
  res = 100
)

par(mar = c(1, 1, 1, 1))

#plot phylo
plot(
  phy,
  show.tip.label = FALSE,
  type = 'fan',
  x.lim = c(-160, 160),
  y.lim = c(-160, 160)
)

#NOTE: linked to chosen clustering method
#add tip and node labels based on tip states and reconstructions
tiplabels(pch = 16,
          col = as.factor(clust_df$ward),
          cex = 2)
nodelabels(pie = HMM_ARD$states,
           piecol = brewer.pal(3, "Set1"),
           cex = 0.3)

#NOTE: linked to chosen clustering method
#add legend
legend(
  x = -160,
  y = 160,
  legend = c(
    "1. Monomorphic herbs",
    "2. Dimorphic herbs & trees",
    "3. Monomorphic trees"
  ),
  fill = brewer.pal(3, "Set1")
)

#NOTE: linked to chosen clustering method
#add transition table
addtable2plot(
  -160,
  120,
  trans_mat,
  bty = "o",
  display.rownames = TRUE,
  hlines = TRUE,
  vlines = TRUE
)

dev.off()
