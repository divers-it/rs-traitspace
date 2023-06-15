library(dplyr)
library(clustMixType)
library(wesanderson)

#load formatted data
df<-readRDS(file = here::here("outputs/df_filt_trans_one_hot.rds"))

#numeric columns only
nums <- unlist(lapply(df, is.numeric))
facts <- unlist(lapply(df, is.factor))

df2<-cbind(df[ , nums],df[ , facts])

str(df2)

ss<-vector()

clust_memb<-vector()

for(i in 2:10){
  kproto_out<-kproto(df2, k=i, lambda = NULL, iter.max = 1000, nstart = 10, na.rm = F)
  ss[i]<-kproto_out$tot.withinss
  
  if(i == 2){
    clust_memb<-kproto_out$cluster
  } else {
    clust_memb<-cbind(clust_memb,kproto_out$cluster)
  }
}

#check alignment of names
kproto_out2<-kproto(df2, k=2, lambda = NULL, iter.max = 1000, nstart = 10, na.rm = F)
kproto_out3<-kproto(df2, k=3, lambda = NULL, iter.max = 1000, nstart = 10, na.rm = F)
names(kproto_out2$cluster)==names(kproto_out3$cluster)
names(kproto_out2$cluster)==rownames(clust_memb)

#plot total ss
plot(ss,type='b')
# k = 4 or 6

#rerun with chosen value of k
kproto_out<-kproto(df2, k=3, lambda = NULL, iter.max = 1000, nstart = 10, na.rm = F)

#relationships of dataset properties to clusters


#removes asking for plot
source("R/myclprofiles.R")

n=ncol(df2[,nums])
png("figures/kproto_cluster_characteristics_quant_one_hot.png",height = 2000, width=1500,res=200)
par(mfrow=c(3,ceiling(n/3)))
myclprofiles(kproto_out, df2[,1:6], col = wes_palette("Royal1", 4, type = "continuous"))
dev.off()

m=ncol(df2[,facts])
png("figures/kproto_cluster_characteristics_qual1_one_hot.png",height = 2000, width=1500,res=200)
par(mfrow=c(3,ceiling(m/6)))
myclprofiles(kproto_out, df2[,(n+1):(n+ceiling(m/2))], col = wes_palette("Royal1", 4, type = "discrete"))
dev.off()

png("figures/kproto_cluster_characteristics_qual2_one_hot.png",height = 2000, width=1500,res=200)
par(mfrow=c(3,ceiling(m/6)))
myclprofiles(kproto_out, df2[,(n+ceiling(m/2)+1):(n+m)], col = wes_palette("Royal1", 4, type = "discrete"))
dev.off()

## --------- PCOA scatterplot with cluster annotation ---------

par(mfrow=c(1,1))

#convert to distance matrix
dataset_dist <- stats::as.dist(gower_df)

#run PCOA
dataset_pcoa <- ape::pcoa(dataset_dist)

#plot points on first two axes, coloured by cluster
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
  xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))

ggsave("figures/pcoa_kproto_k3_one_hot.png",width = 12,height=10)


####
# Sankey plot
####

#from: https://r-graph-gallery.com/321-introduction-to-interactive-sankey-diagram-2.html
#table of different k values (2-7)

for (i in 1:6) {
  if (i == 1) {
    clust.num.k.2.7 <- paste("k",i+1,"cluster",as.character(clust_memb[,i]),sep="_")
  } else {
    clust.num.k.2.7 <-
      cbind(clust.num.k.2.7, paste("k",i+1,"cluster",as.character(clust_memb[,i]),sep="_"))
  }
}

clust_memb

colnames(clust.num.k.2.7)<-c("2clusters",
                             "3clusters",
                             "4clusters",
                             "5clusters",
                             "6clusters",
                             "7clusters")
clust.num.k.2.7.df <-as.data.frame(clust.num.k.2.7)

rownames(clust.num.k.2.7.df)<-rownames(clust_memb)

saveRDS(clust.num.k.2.7.df, file = here::here("outputs/clust_num_k_2_7_kpro_one_hot.rds"))

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

#remove rows where values are 0
links<-links[links$value>0,]

# Library
library(networkD3)
library(dplyr)
# Make the Network

p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE)

p

saveNetwork(p, "figures/sankey_kpro_one_hot.html")


###
# Robust combinations
###

#make data frame of combo frequencies
combos <- as.data.frame(table(clust.num.k.2.7.df))

#remove no existant combos
combos <- combos[combos$Freq > 0, ]

#order
combos <- combos[order(combos$Freq, decreasing = T), ]

#change to strings
combos <-data.frame(lapply(combos, as.character), stringsAsFactors = FALSE)

#table of different combinations and their frequencies
combos

#empty list
robust<-list()

#loop through ordered table to extract robust groups
#change value in loop for threshold
for(i in 1:length(combos$Freq[as.numeric(combos$Freq)>20])){
  robust[[i]]<-rownames(clust.num.k.2.7.df[clust.num.k.2.7.df[, 1] == combos[i, 1] & 
                                                clust.num.k.2.7.df[, 2] == combos[i, 2] &
                                                clust.num.k.2.7.df[, 3] == combos[i, 3] &
                                                clust.num.k.2.7.df[, 4] == combos[i, 4] &
                                                clust.num.k.2.7.df[, 5] == combos[i, 5] &
                                                clust.num.k.2.7.df[, 6] == combos[i, 6],])
}

robust



