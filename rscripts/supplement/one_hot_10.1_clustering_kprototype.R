rm(list=ls())

# load packages
library(dplyr)
library(clustMixType)
library(wesanderson)
library(ggplot2)
library(cluster)
library(networkD3)
library(gridExtra)
library(patchwork)

## load formatted data
df<-readRDS(file = here::here("outputs/one_hot_6_df_filt_trans.rds"))

#### 
## K-prototypes clustering ----
#### 
# 
## UNCOMMENT TO RERUN CLUSTERING
# set up empty vectors
ss <- vector()
clust_memb <- vector()
kproto_list <- list()

# run clustering with different values of K up to 10
for (i in 2:10) {
  kproto_out <-
    kproto(
      df,
      k = i,
      lambda = NULL,
      iter.max = 1000,
      nstart = 10,
      na.rm = F
    )
  
  kproto_list[[i]] <- kproto_out
  
  ss[i] <- kproto_out$tot.withinss
  
  if (i == 2) {
    clust_memb <- kproto_out$cluster
  } else {
    clust_memb <- cbind(clust_memb, kproto_out$cluster)
  }
  
}

# save image as takes long to run
save.image(file = "outputs/one_hot_10.1_kpro.RData")

# load previous Kproto run
load("outputs/one_hot_10.1_kpro.RData")

### Select number of clusters (k) ----

# look at clustering output
head(clust_memb)

# check alignment of names
names(kproto_list[[2]]$cluster)==rownames(clust_memb)

# plot total ss to choose 
plot(ss,type='b')

# chosen value of k
kproto_out<-kproto_list[[5]]

####
## PCOA scatterplot with cluster annotation ----
####

# dissimilarity matrix calculation
gower_df <- daisy(df,
                  metric = "gower" )

# convert to distance matrix
dataset_dist <- stats::as.dist(gower_df)

# run PCOA
dataset_pcoa <- ape::pcoa(dataset_dist)

# Recalculate relative eigenvalues by removing negative eigenvalues as in Mouillot et al.  
ev_pcoa <- dataset_pcoa$values$Eigenvalues
ev_pcoa_g0 <- ev_pcoa[ev_pcoa>0]
rel_ev_pcoa_g0 <- ev_pcoa_g0/sum(ev_pcoa_g0)

# plot points on first two axes, coloured by cluster
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2, fill = as.factor(kproto_out$cluster))) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.5,
    size=3,
    stroke = 0.5
  ) + 
  stat_ellipse(geom = "polygon",
               aes(fill = as.factor(kproto_out$cluster)), 
               alpha = 0.25) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2))) + 
  theme(legend.position = "none")

# save plot
# ggsave("figures/10.1_scatterplot_pcoa_kpro_k5_coloured_by_cluster.png",width = 10,height=10)

### 
# ---- Sankey plot ----
### 

# from: https://r-graph-gallery.com/321-introduction-to-interactive-sankey-diagram-2.html
# table of different k values (2-7)

for (i in 1:6) {
  if (i == 1) {
    clust.num.k.2.7 <- paste("k",i+1,"cluster",as.character(clust_memb[,i]),sep="_")
  } else {
    clust.num.k.2.7 <-
      cbind(clust.num.k.2.7, paste("k",i+1,"cluster",as.character(clust_memb[,i]),sep="_"))
  }
}

colnames(clust.num.k.2.7)<-c("2clusters",
                             "3clusters",
                             "4clusters",
                             "5clusters",
                             "6clusters",
                             "7clusters")
clust.num.k.2.7.df <-as.data.frame(clust.num.k.2.7)

rownames(clust.num.k.2.7.df)<-rownames(clust_memb)

saveRDS(clust.num.k.2.7.df, file = here::here("outputs/one_hot_10.1_clust_num_k_2_7_kpro.rds"))

# A connection data frame is a list of flows with intensity for each flow

for(i in 1:(length(colnames(clust.num.k.2.7.df))-1)){
  if(i == 1){
    links<-as.data.frame(table(clust.num.k.2.7.df[,c(i,(i+1))]))
    colnames(links)<-c("source","target","value")
  } else {
    
    tmp<-as.data.frame(table(clust.num.k.2.7.df[,c(i,(i+1))]))
    colnames(tmp)<-c("source","target","value")
    links <-
      rbind(links, tmp)
  }
  
}

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1
links

# remove rows where values are 0
links<-links[links$value>0,]


# Make and plot the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", fontSize = 14,
                   sinksRight=FALSE)
p

# save as html
# saveNetwork(p, "figures/one_hot_10.1_sankey_kpro.html")

####
## Identify robust groups ----
####

# make data frame of combo frequencies
combos <- as.data.frame(table(clust.num.k.2.7.df))

# remove no existant combos
combos <- combos[combos$Freq > 0, ]

# order
combos <- combos[order(combos$Freq, decreasing = T), ]

# change to strings
combos <-data.frame(lapply(combos, as.character), stringsAsFactors = FALSE)
plot(combos$Freq)

# empty list
robust<-list()

# empty vector
robust_vect_kpro<-rep(NA,length(rownames(df)))
names(robust_vect_kpro)<-rownames(df)

# loop through ordered table to extract robust groups
for(i in 1:length(combos$Freq[as.numeric(combos$Freq)>0])){
  foo<-rownames(clust.num.k.2.7.df[clust.num.k.2.7.df[, 1] == combos[i, 1] & 
                                     clust.num.k.2.7.df[, 2] == combos[i, 2] &
                                     clust.num.k.2.7.df[, 3] == combos[i, 3] &
                                     clust.num.k.2.7.df[, 4] == combos[i, 4] &
                                     clust.num.k.2.7.df[, 5] == combos[i, 5] &
                                     clust.num.k.2.7.df[, 6] == combos[i, 6],])
  
  robust[[i]]<-foo
  
  robust_vect_kpro[foo]<-i
  
}

# robust groups
robust

# keep 80% of the species in robust clusters ; others are NA
sum = 0

for (i in 1:max(robust_vect_kpro)) {
  
  if (sum > 0.80 * length(robust_vect_kpro)) {
    robust_vect_kpro[robust_vect_kpro == i] <- NA
    
  }
  
  sum = sum + length(robust[[i]])
  
}

table(robust_vect_kpro)

# complete vector of robust groups and non-robust 
robust_vect_kpro_full<-robust_vect_kpro

saveRDS(robust_vect_kpro_full, file = here::here("outputs/one_hot_10.1_robust_vect_kpro_full.rds"))

#remove species not in robust groups
robust_vect_kpro<-na.omit(robust_vect_kpro)

# check order
rownames(dataset_pcoa$vectors)==rownames(clust.num.k.2.7.df)

# Plot robust groups on PCoA scatterplot
# plot points on first two axes, coloured by robust group, shaped by cluster
ggplot(
  data.frame(dataset_pcoa$vectors),
  aes(
    x = Axis.1,
    y = Axis.2,
    col = as.factor(robust_vect_kpro_full)
  )
) +
  geom_point(
    aes(shape = as.factor(clust.num.k.2.7.df$`5clusters`)),
    alpha = 0.5,
    size = 3,
    stroke = 0.5
  )

# save plot
# ggsave("figures/10.1_scatterplot_pcoa_kpro_k5_coloured_by_robust.png",width=12,height=10)

