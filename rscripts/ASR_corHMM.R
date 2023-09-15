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

#NOTE: Uncomment to read in transition rates from previous run
df_rates<-readRDS(file = here::here("outputs/one_hot_transition_rates.rds"))

# #loop through all characters
# for (i in 2:length(colnames(df2))) {
#   #make new dataframe with only trait of interest
#   df3 <- df2[, c(1, i)]
#   
#   #omit missing data
#   df3 <- na.omit(df3)
#   
#   #drop tips not in dataset
#   phy2 <- drop.tip(phy, setdiff(phy$tip.label, df3$species))
#   
#   #sort df to match order of tips in phylo
#   df3 <- df3[match(phy2$tip.label, df3$species), ]
#   
#   #plot phylogeny and example trait
#   plot(phy2, show.tip.label = FALSE, type = 'fan')
#   tiplabels(pch = 16,
#             col = as.factor(df3[, 2]),
#             cex = 1)
#   
#   #fit model with 1 rate category on woodiness trait
#   MK_2state <- corHMM(
#     phy = phy2,
#     data = df3,
#     rate.cat = 1,
#     model = "ER"
#   )
#   
#   states[i - 1] <- colnames(df3)[2]
#   rates[i - 1] <- MK_2state$solution[1, 2]
#   tree_sizes[i - 1] <- length(phy2$tip.label)
#   no_states[i - 1] <- length(unique(df3[, 2]))
# }
# 
# df_rates <- data.frame(states, rates, tree_sizes, no_states)
 

#NOTE: uncomment to save/load transition rates
#write.csv(df_rates, file = here::here("outputs/one_hot_transition_rates.csv"))
df_rates<-read.csv(file = here::here("outputs/one_hot_transition_rates.csv"),row.names = 1)

#ordered results
df_rates[order(df_rates$rates, decreasing = T), ]

###
# ---- Compare phylogenetic signal and transition rates ----
###

#read in phylogenetic signal results for qualitiative traits
psq<-read.csv(file = here::here("outputs/phylo_signal_qualitative.csv"),row.names = 1)

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

####
# ---- Figure 5: Reproductive strategy ASR phylogenetic tree ----
####

library(wesanderson)

#plot tree with marginal probabilites on nodes
png(
  "figures/phylo_asr_clusters.png",
  width = 1500,
  height = 1500,
  res = 100
)

par(mar = c(1, 1, 1, 1))
palette(wes_palette("IsleofDogs1", 3))

#plot phylo
plot(
  phy,
  show.tip.label = FALSE,
  type = 'fan',
  x.lim = c(-150, 150),
  y.lim = c(-150, 150)
)

#NOTE: linked to chosen clustering method
#add tip and node labels based on tip states and reconstructions
tiplabels(pch = 21,
          col="black",
          bg=as.factor(clust_df$ward),
          piecol = as.factor(clust_df$ward),
          cex = 1.25)
nodelabels(pie = HMM_ARD$states,
           piecol = wes_palette("IsleofDogs1", 3),
           cex = 0.25)

#NOTE: linked to chosen clustering method
#add legend
legend(
  x = -151,
  y = 150,
  legend = c(
    "1. Monomorphic herbaceous",
    "2. Dimorphic",
    "3. Monomorphic woody"
  ),
  fill = wes_palette("IsleofDogs1", 3),
  bty = "n",
  cex=1.3
)

#NOTE: linked to chosen clustering method
#add transition table
addtable2plot(
  -147,
  107,
  trans_mat,
  bty = "o",
  display.rownames = TRUE,
  hlines = TRUE,
  vlines = TRUE,
  cex=1.2
)

dev.off()

library(ggtree)
ggtree(phy, layout="fan") + 
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + 
  geom_cladelab(node=379, label="test label", angle=0, fontsize=8, offset=.5, vjust=.5)



ancstats <- as.data.frame(HMM_ARD$states)
ancstats$node <- 1:phy$Nnode+Ntip(phy)

## cols parameter indicate which columns store stats
bars <- nodebar(ancstats, cols=1:3)
bars <- lapply(bars, function(g) g+scale_fill_manual(values=wes_palette("IsleofDogs1", 3)))

tree2 <- full_join(phy, data.frame(label = rownames(df), stat = as.character(clust_df$ward) ), by = 'label')

ggtree(tree2)

p <- ggtree(tree2,layout = "fan") + geom_tiplab() +
  geom_tippoint(aes(color = stat)) +
  scale_color_manual(values=wes_palette("IsleofDogs1", 3)[1:3]) +
  theme(legend.position = "right") + 
  xlim(NA, 8)
p

p1 <- p + geom_inset(bars, width = .08, height = .05, x = "branch")   
p1

library(phytools)
data(anoletree)
x <- getStates(anoletree,"tips")
tree <- as.phylo(anoletree)

cols <- setNames(palette()[1:length(unique(x))],sort(unique(x)))
fitER <- ape::ace(x,tree,model="ER",type="discrete")
ancstats <- as.data.frame(fitER$lik.anc)
ancstats$node <- 1:tree$Nnode+Ntip(tree)

## cols parameter indicate which columns store stats
bars <- nodebar(ancstats, cols=1:6)
bars <- lapply(bars, function(g) g+scale_fill_manual(values = cols))

tree2 <- full_join(tree, data.frame(label = names(x), stat = x ), by = 'label')
p <- ggtree(tree2) + geom_tiplab() +
  geom_tippoint(aes(color = stat)) + 
  scale_color_manual(values = cols) +
  theme(legend.position = "right") + 
  xlim(NA, 8)
p1 <- p + geom_inset(bars, width = .08, height = .05, x = "branch")   
p1
