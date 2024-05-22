rm(list = ls())

#load libraries
library(dplyr)
library(ape)
library(corHMM)
library(RColorBrewer)
library(phytools)
library(plotrix)
library(ggplot2)
library(ggtree)
library(ggpp)
library(rphylopic)

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

####
# ---- MK 2-state ----
####

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
 


####
# ---- Mk + hidden states (1 rate cat) ----
####

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
#   MK_hidd <- corHMM(
#     phy = phy2,
#     data = df3,
#     rate.cat = 2,
#     model = "ARD"
#   )
#   
#   states[i - 1] <- colnames(df3)[2]
#   rates[i - 1] <- MK_hidd$solution[1, 2]
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


#make labels removing low values 
d1d2$labels <- d1d2$states
d1d2$labels[d1d2$deltas<7 & d1d2$rates<0.005]<-NA


#plot signal vs transition rates
ggplot(d1d2, aes(x = deltas, y = rates)) +
  geom_point(
    fill=wesanderson::wes_palette("Royal1",2)[2],
    shape=21,
    alpha=0.4,
    size=5,
    stroke = 0.5
  ) + theme_bw() +
  ggrepel::geom_text_repel(aes(label = d1d2$labels),size = 5,max.overlaps = Inf, nudge_x = 3, nudge_y = 0.0001) +
  xlab("Phylogenetic signal (Delta)") +
  ylab("Transition rates")

ggsave("figures/scatterplot_signal_transition_rate.png",height=10, width=10)

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

#make model
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

####
# Hidden rates model 
####

# #Add another rate category
# HMM_ARD2 <-
#   corHMM(
#     phy = phy,
#     data = dat,
#     rate.cat = 2,
#     model = "ARD",
#     node.states = 'marginal',
#     get.tip.states = TRUE
#   )
# 
# #examine model
# HMM_ARD2
# 
# #plot transition rates
# plotMKmodel(HMM_ARD2)
# 
# #compare AIC
# print("AIC of simple Mk model")
# HMM_ARD$AIC
# print("AIC of model with hidden rates")
# HMM_ARD2$AIC

#NOTE: linked to chosen clustering method
#make transition rates table
trans_mat <- data.frame(round(HMM_ARD$solution, 4))
colnames(trans_mat)<-c("->1","->2","->3")
rownames(trans_mat)<-c("1->","2->","3->")

####
# ---- Figure 5: Reproductive strategy ASR phylogenetic tree ----
####

library(wesanderson)

#get taxonomy
tax<-readRDS("outputs/taxonomy.rds")

str(tax)

tax<-tax[match(phy$tip.label,gsub(" ","_",rownames(tax))),]

#make new tree for plotting
phy2<-phy

#make names unique
phy2$tip.label<-make.unique(tax$order)

#plot tree with marginal probabilites on nodes
png(
  "figures/phylo_asr_clusters.png",
  width = 2500,
  height = 2500,
  res = 180
)

par(bg=NA)
par(mar = c(1, 1, 1, 1))
palette(wes_palette("IsleofDogs1", 3))

#plot phylo
plot(
  phy2,
  show.tip.label = FALSE,
  #cex=0.5,
  type = 'fan',
  x.lim = c(-160, 160),
  y.lim = c(-160, 160)
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
  x = -166,
  y = 165,
  legend = c(
    "1. Bisexual herbaceous",
    "2. Unisexual",
    "3. Bisexual woody"
  ),
  pt.bg = wes_palette("IsleofDogs1", 3),
  bty = "n",
  cex=1.3,
  pch=21,
  pt.cex=2,
  #col = wes_palette("IsleofDogs1", 3)
  col="black"
)

#NOTE: linked to chosen clustering method
#add transition table
addtable2plot(
  -162,
  117,
  trans_mat,
  bty = "o",
  display.rownames = TRUE,
  hlines = TRUE,
  vlines = TRUE,
  cex=1.2
)

dev.off()

####
# ---- Phylopics ----
####

# Get a single image uuid
uuid <- get_uuid(name = "Primula", n = 1)
# Get the image for that uuid
ericales_pp <- get_phylopic(uuid = uuid[1])

uuid <- get_uuid(name = "Olea", n = 1)
lamiales_pp <- get_phylopic(uuid = uuid[1])

uuid <- get_uuid(name = "Crepis", n = 1)
asterales_pp <- get_phylopic(uuid = uuid[1])

uuid <- get_uuid(name = "Magnolia", n = 1)
magnoliales_pp <- get_phylopic(uuid = uuid[1])

uuid <- get_uuid(name = "Panicum miliaceum", n = 1)
poales_pp <- get_phylopic(uuid = uuid[1])

uuid <- get_uuid(name = "Cypripedium", n = 1)
asparagales_pp <- get_phylopic(uuid = uuid[1])

uuid <- get_uuid(name = "Cynomorium", n = 1)
saxifragales_pp <- get_phylopic(uuid = uuid[1])

uuid <- get_uuid(name = "Acer platanoides", n = 1)
sapindales_pp <- get_phylopic(uuid = uuid[1])

uuid <- get_uuid(name = "Arabidopsis thaliana", n = 1)
brassicales_pp <- get_phylopic(uuid = uuid[1])

uuid <- get_uuid(name = "Passiflora incarnata", n = 1)
malpighiales_pp <- get_phylopic(uuid = uuid[1])

uuid <- get_uuid(name = "Ficus carica", n = 1)
rosales_pp <- get_phylopic(uuid = uuid[1])

uuid <- get_uuid(name = "Nepenthes distillatoria", n = 1)
caryophyllales_pp <- get_phylopic(uuid = uuid[1])

uuid <- get_uuid(name = "Quercus robur", n = 2)
fagales_pp <- get_phylopic(uuid = uuid[2])

uuid <- get_uuid(name = "Capsicum annuum", n = 1)
solanales_pp <- get_phylopic(uuid = uuid[1])

uuid <- get_uuid(name = "Anthurium", n = 2)
alismatales_pp <- get_phylopic(uuid = uuid[2])

####
# ---- GGTREE ----
####

x2<-paste("Cluster",clust_df$ward,sep="")
names(x2)<-gsub(" ","_",rownames(df))

x2[phy$tip.label]

ancstats2 <- as.data.frame(HMM_ARD$states)
ancstats2$node <- 1:phy$Nnode+Ntip(phy)
colnames(ancstats2) <- c("Cluster1","Cluster2","Cluster3","node")

ggtree(phy, linewidth=0.5) +
  geom_tiplab(size=2)

#get taxonomy
tax <- readRDS("outputs/taxonomy.rds")

#read in phylogenetic tree
phy <- read.tree("outputs/pruned_tree.tre")

rownames(tax)<-gsub(" ","_",rownames(tax))

clust_df <- clust_df[match(phy$tip.label,clust_df$species),]
tax <- tax[match(phy$tip.label,rownames(tax)),]

metadata<-data.frame(rownames(tax),clust_df$ward,tax$order)
colnames(metadata)<-c("label","cluster","order")

phy$tip.label==metadata$label


#Make dataframe for clade nodes
clades.df <- data.frame(
  clade=unique(metadata$order),
  node=NA
)

#Find the most recent common ancestor for each clade
for (i in 1:length(clades.df$clade)) {
  
  clades.df$node[i] <- MRCA(
    phy,
    metadata$label[metadata$order == clades.df$clade[i]]
  )
  
}


#Add highlights
g1 <- ggtree(phy, linetype=NA) %<+% metadata +
  geom_highlight(data=clades.df, 
                 aes(node=node, fill=clade),
                 alpha=1,
                 align="right",
                 extend=0.1,
                 show.legend=FALSE) +
  geom_tree(linewidth=0.5) +
  geom_tiplab(aes(label=label), size=2)

#Order the clades dataframe to match the tree
clades.df <- clades.df[match(g1$data %>%
                               filter(isTip == "TRUE") %>%
                               arrange(y) %>%
                               pull(order) %>%
                               unique(),
                             clades.df$clade),]

#Add column with alternating binary value
clades.df$highlight <- rep(c(0,1),
                           length.out=length(clades.df$clade))
head(clades.df)

#Add clade labels
g1 <- ggtree(phy, layout="circular", linetype=NA) %<+% metadata +
  geom_highlight(data=clades.df, 
                 aes(node=node, fill=as.factor(highlight)),
                 alpha=1,
                 align="right",
                 extend=0.04,
                 show.legend=FALSE) +
  geom_cladelab(data=clades.df,
                mapping=aes(node=node, label=clade),
                fontsize=3,
                align="TRUE",
                angle="auto",
                offset=0.04,
                offset.text=0.01) +
  geom_tree(linewidth=0.3) +
  geom_tippoint(mapping=aes(color=as.factor(cluster)), 
                size=1.5) +
  xlim(-40,190) +
  scale_fill_manual(values=c("#ECECEC", "#FCFCFC")) +
  theme(legend.position = c(0.49, 0.50),
        legend.text = element_text(size=11),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key=element_blank(),
        legend.spacing.y = unit(-0.2, "cm"),
        ) +
  guides(colour = guide_legend(override.aes = list(size=5),
                               byrow = TRUE))

g1

colours<-harrypotter::hp(3,option="ronweasley2")

pies <- nodepie(ancstats2, cols=1:3, col=colours)

g2 <- g1 + geom_plot(data=td_filter(!isTip), 
               mapping=aes(x=x,y=y, label=pies),
               vp.width=0.0125, 
               vp.height=0.0125, 
               hjust=0.6,
               vjust=0.6
               ) + scale_color_manual(values=colours,
                                      labels=c('G1: Bisexual herbaceous', 'G2: Unisexual', 'G3: Bisexual woody'))

#no. species = 360 for plotting y coord
g2 +
  add_phylopic(img=magnoliales_pp, x=195, y=14, ysize = 9,col = "grey30") +
  add_phylopic(img=alismatales_pp, x=195, y=29, ysize = 10,col = "grey30") +
  add_phylopic(img=asparagales_pp, x=195, y=55, ysize = 10,col = "grey30") +
  add_phylopic(img=poales_pp, x=195, y=78, ysize = 10,col = "grey30") +
  add_phylopic(img=saxifragales_pp, x=195, y=106, ysize = 10,col = "grey30") +
  add_phylopic(img=sapindales_pp, x=195, y=128, ysize = 10,col = "grey30") +
  add_phylopic(img=brassicales_pp, x=195, y=151, ysize = 10,col = "grey30") +
  add_phylopic(img=rosales_pp, x=195, y=170, ysize = 7,col = "grey30") +
  add_phylopic(img=fagales_pp, x=195, y=185, ysize = 10,col = "grey30") +
  add_phylopic(img=malpighiales_pp, x=195, y=208, ysize = 13,col = "grey30") +
  add_phylopic(img=caryophyllales_pp, x=195, y=244, ysize = 10,col = "grey30") +
  add_phylopic(img=ericales_pp, x=195, y=272, ysize = 10,col = "grey30") + 
  add_phylopic(img=asterales_pp, x=195, y=307, ysize = 10,col = "grey30") +
  add_phylopic(img=solanales_pp, x=195, y=330, ysize = 10,col = "grey30") +
  add_phylopic(img=lamiales_pp, x=195, y=348, ysize = 12,col = "grey30")

ggsave("figures/figure3_tree2.png", width = 15, height=15)
