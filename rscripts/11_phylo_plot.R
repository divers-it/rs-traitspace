rm(list = ls())

# load libraries
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

# Load formatted data
df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

# numeric/factor columns only
nums <- unlist(lapply(df, is.numeric))
facts <- unlist(lapply(df, is.factor))

# make only categorical dataset
df2 <- df[, facts]

# insert '_' into rownames to match phylo
df2$species <- gsub(" ", "_", rownames(df2))

# remove rownames
rownames(df2) <- NULL

# put species column first
df2 <- df2[, c(length(colnames(df2)), 1:(length(colnames(df2)) - 1))]

# read in phylogenetic tree
phy <- read.tree("outputs/pruned_tree.tre")
plot(phy, cex = 0.5)

# in dataset but not in phylo
setdiff(df2$species, phy$tip.label)

# in phylo but not in dataset
setdiff(phy$tip.label, df2$species)

### 
##  Read and format trait groupings ----
### 

# get rid of polymorphic states
wood <- df2$Woodiness
wood[wood == "herbaceous_woody"] <- NA
wood <- droplevels(wood)
table(wood)

# get rid of polymorphic states
fs <- df2$FlowerSex
fs[fs == "bisexual_unisexual"] <- NA
fs <- droplevels(fs)
table(fs)

# make data frame
trait_df <- data.frame(df2$species, paste(wood,fs))
colnames(trait_df) <- c("species","trait")

# NA species
trait_df$trait[grep("NA",trait_df$trait)] <- NA

# remove rownames
rownames(trait_df) <- NULL

# read in phylogenetic tree
phy <- read.tree("outputs/pruned_tree.tre")

# in dataset but not in phylo
setdiff(trait_df$species, phy$tip.label)

# in phylo but not in dataset
setdiff(phy$tip.label, trait_df$species)

# drop tips
phy <- drop.tip(phy, setdiff(phy$tip.label, trait_df$species))

# check order
trait_df <- trait_df[match(phy$tip.label, trait_df$species), ]
phy$tip.label == trait_df$species

# NOTE: change column to choose different classification
# make associated data frame by choosing clusters/robust column
dat <- trait_df


#### 
## Figure 5: Reproductive strategy on phylogenetic tree ----
#### 

# get taxonomy
tax<-readRDS("outputs/taxonomy.rds")

str(tax)

tax<-tax[match(phy$tip.label,gsub(" ","_",rownames(tax))),]

# make new tree for plotting
phy2<-phy

# make names unique
phy2$tip.label<-make.unique(tax$order)

#### 
### Phylopics ----
#### 

#  Get a single image uuid
uuid <- get_uuid(name = "Primula", n = 1)
#  Get the image for that uuid
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
### GGTREE ----
#### 

x2<-trait_df$trait
names(x2)<-trait_df$species

x2[phy$tip.label]

ggtree(phy, linewidth=0.5) +
  geom_tiplab(size=2)

# get taxonomy
tax <- readRDS("outputs/taxonomy.rds")

# read in phylogenetic tree
phy <- read.tree("outputs/pruned_tree.tre")

rownames(tax)<-gsub(" ","_",rownames(tax))

trait_df <- trait_df[match(phy$tip.label,trait_df$species),]
tax <- tax[match(phy$tip.label,rownames(tax)),]

metadata<-data.frame(rownames(tax),trait_df$trait,tax$order)
colnames(metadata)<-c("label","trait","order")

phy$tip.label==metadata$label

# Make dataframe for clade nodes
clades.df <- data.frame(
  clade=unique(metadata$order),
  node=NA
)

# Find the most recent common ancestor for each clade
for (i in 1:length(clades.df$clade)) {
  
  clades.df$node[i] <- tidytree::MRCA(
    phy,
    metadata$label[metadata$order == clades.df$clade[i]]
  )
  
}


# Add highlights
g1 <- ggtree(phy, linetype=NA) %<+% metadata +
  geom_highlight(data=clades.df, 
                 aes(node=node, fill=clade),
                 alpha=1,
                 align="right",
                 extend=0.1,
                 show.legend=FALSE) +
  geom_tree(linewidth=0.5) +
  geom_tiplab(aes(label=label), size=2)

# Order the clades dataframe to match the tree
clades.df <- clades.df[match(g1$data %>%
                               filter(isTip == "TRUE") %>%
                               arrange(y) %>%
                               pull(order) %>%
                               unique(),
                             clades.df$clade),]

# Add column with alternating binary value
clades.df$highlight <- rep(c(0,1),
                           length.out=length(clades.df$clade))
head(clades.df)

# set colours
colours<-harrypotter::hp(3,option="ronweasley2")
colours <- c(colorRampPalette(c("lightblue","navyblue"))(3)[1:2],"darkred","#E7D889")

clades.df$highlight <- as.factor(clades.df$highlight)

# Add clade labels
g1 <- ggtree(phy, layout="circular", linetype=NA) %<+% metadata +
  geom_highlight(data=clades.df,
                 mapping=aes(node=node, fill=highlight),
                  alpha=0.75,
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
  geom_tippoint(mapping=aes(color=as.factor(trait)), 
                size=1.5) +
  xlim(-50,200) +
  scale_fill_manual(values=c("#ECECEC", "#FCFCFC")) +
  theme(legend.position = c(0.53, 0.455),
        legend.text = element_text(size=11),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key=element_blank(),
        legend.spacing.y = unit(-0.2, "cm"),
        ) + scale_color_manual(values=colours,
                               na.value = "white",
                               labels=c('Herbaceous + bisexual flowers', 'Herbaceous + unisexual flowers', 'Woody + bisexual flowers','Woody + unisexual flowers','')) +
    guides(colour = guide_legend(override.aes = list(size=5),
                               byrow = TRUE))

g1

#  no. species = 361 for plotting y coord
g1 +
  add_phylopic(img=magnoliales_pp, x=195, y=14, ysize = 9,col = "black") +
  add_phylopic(img=alismatales_pp, x=195, y=29, ysize = 10,col = "black") +
  add_phylopic(img=asparagales_pp, x=195, y=55, ysize = 10,col = "black") +
  add_phylopic(img=poales_pp, x=195, y=78, ysize = 10,col = "black") +
  add_phylopic(img=saxifragales_pp, x=195, y=106, ysize = 10,col = "black") +
  add_phylopic(img=sapindales_pp, x=195, y=128, ysize = 10,col = "black") +
  add_phylopic(img=brassicales_pp, x=195, y=151, ysize = 10,col = "black") +
  add_phylopic(img=rosales_pp, x=195, y=170, ysize = 7,col = "black") +
  add_phylopic(img=fagales_pp, x=195, y=185, ysize = 10,col = "black") +
  add_phylopic(img=malpighiales_pp, x=195, y=208, ysize = 13,col = "black") +
  add_phylopic(img=caryophyllales_pp, x=195, y=244, ysize = 10,col = "black") +
  add_phylopic(img=ericales_pp, x=195, y=272, ysize = 10,col = "black") + 
  add_phylopic(img=asterales_pp, x=195, y=307, ysize = 10,col = "black") +
  add_phylopic(img=solanales_pp, x=195, y=330, ysize = 10,col = "black") +
  add_phylopic(img=lamiales_pp, x=195, y=348, ysize = 12,col = "black")

ggsave("figures/figure_1_phylo.png", width = 15, height=15)
