rm(list = ls())

# load packages
library(dplyr)
library(gridExtra)
library(Rtsne)
library(ggplot2)
library(cluster)
library(RColorBrewer)
library(ggrepel)
library(networkD3)

# load formatted data
df<-readRDS(file = here::here("outputs/one_hot_6_df_filt_trans.rds"))

# build distance matrix
gower_df <- daisy(df,
                  metric = "gower" )
summary(gower_df)

## Select number of clusters (k) ----

# Silhouette Width to select the optimal number of clusters
# The silhouette width is one of the very popular choices when it comes to selecting the optimal number of clusters. 
# It measures the similarity of each point to its cluster, and compares that to the similarity of the point with the closest neighboring cluster. 
# This metric ranges between -1 to 1, where a higher value implies better similarity of the points to their clusters. 
# Therefore, a higher value of the Silhouette Width is desirable. 
# We calculate this metric for a range of cluster numbers and find where it is maximized. 
# The following code shows the implementation in R:

#empty vector
silhouette <- c()

#run PAM with different values of K from 2-10 and calculate silhouette width
for(i in 2:10){
  pam_clusters <- pam(as.matrix(gower_df),
                      diss = TRUE,
                      k = i)
  silhouette <- c(silhouette ,pam_clusters$silinfo$avg.width)
}

### Plot silhouette width ----
plot(2:10, silhouette,
     xlab = "Clusters",
     ylab = "Silhouette Width")
lines(2:10, silhouette)

# make df with cluster membership for each value of k from 2-7
for(i in 2:7){
  if(i == 2){
    pam.gower = pam(gower_df, diss = TRUE, k = i)
    pam_df<-pam.gower$clustering
  } else {
    pam.gower = pam(gower_df, diss = TRUE, k = i)
    pam_df<-cbind(pam_df,pam.gower$clustering)
  }
}

colnames(pam_df)<-c("2clusters",
                    "3clusters",
                    "4clusters",
                    "5clusters",
                    "6clusters",
                    "7clusters")
clust.num.k.2.7.df <-as.data.frame(pam_df)
rownames(clust.num.k.2.7.df)<-names(pam.gower$clustering)

####
## Sankey plot ----
####

#from: https://r-graph-gallery.com/321-introduction-to-interactive-sankey-diagram-2.html

#table of different k values (2-7)
for (i in 1:6) {
  if (i == 1) {
    clust.num.k.2.7 <- paste("k",i+1,"cluster",as.character(clust.num.k.2.7.df[,i]),sep="_")
  } else {
    clust.num.k.2.7 <-
      cbind(clust.num.k.2.7, paste("k",i+1,"cluster",as.character(clust.num.k.2.7.df[,i]),sep="_"))
  }
}

colnames(clust.num.k.2.7)<-c("2clusters",
                             "3clusters",
                             "4clusters",
                             "5clusters",
                             "6clusters",
                             "7clusters")
clust.num.k.2.7.df <-as.data.frame(clust.num.k.2.7)

# fix rownames again
rownames(clust.num.k.2.7.df)<-names(pam.gower$clustering)

# save RDS
saveRDS(clust.num.k.2.7.df, file = here::here("outputs/one_hot_10_clust_num_k_2_7_pam.rds"))

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

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", fontSize = 14,
                   sinksRight=FALSE)

p

# save as HTML
saveNetwork(p, "figures/10_sankey_pam.html")

####
## Identify robust groups ----
####
# Robust groups are those that consistently stay together as values of k change

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
robust_vect_pam<-rep(NA,length(rownames(df)))
names(robust_vect_pam)<-rownames(df)

# loop through ordered table to extract robust groups
for(i in 1:length(combos$Freq[as.numeric(combos$Freq)>0])){
  foo<-rownames(clust.num.k.2.7.df[clust.num.k.2.7.df[, 1] == combos[i, 1] & 
                                     clust.num.k.2.7.df[, 2] == combos[i, 2] &
                                     clust.num.k.2.7.df[, 3] == combos[i, 3] &
                                     clust.num.k.2.7.df[, 4] == combos[i, 4] &
                                     clust.num.k.2.7.df[, 5] == combos[i, 5] &
                                     clust.num.k.2.7.df[, 6] == combos[i, 6],])
  
  robust[[i]]<-foo
  
  robust_vect_pam[foo]<-i
  
}

# robust groups
robust

# keep 80% of the species in robust clusters ; others are NA
sum = 0

for (i in 1:max(robust_vect_pam)) {
  
  if (sum > 0.8 * length(robust_vect_pam)) {
    robust_vect_pam[robust_vect_pam == i] <- NA
    
  }
  
  sum = sum + length(robust[[i]])
  
}

table(robust_vect_pam)

# complete vector of robust groups and non-robust 
robust_vect_pam_full<-robust_vect_pam

saveRDS(robust_vect_pam_full, file = here::here("outputs/one_hot_10_robust_vect_pam_full.rds"))
