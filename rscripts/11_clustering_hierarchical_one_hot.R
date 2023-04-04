rm(list = ls())
library(dplyr)

#load formatted data
df2<-readRDS(file = here::here("outputs/df_filt_trans_one_hot.rds"))
df2[sapply(df2, is.factor)] <- lapply(df2[sapply(df2, is.factor)],
                                     as.integer)

#dissimilarity matrix
library(cluster)
gower_df2 <- daisy(df2, metric = "gower" )

summary(gower_df2)

# The main input for the code below is dissimilarity (distance matrix)

## ------------ DIVISIVE CLUSTERING ------------

library(cluster)
divisive.clust <- diana(as.matrix(gower_df2),
                        diss = TRUE, keep.diss = TRUE)

pltree(divisive.clust, main = "Divisive", cex = 0.25)

## ------------ AGGLOMERATIVE CLUSTERING ------------

#running different methods to compare results

aggl.clust.c <- hclust(gower_df2, method = "complete")
plot(aggl.clust.c, main = "Agglomerative, complete linkages", cex = 0.25)

aggl.clust.s <- hclust(gower_df2, method = "single")
plot(aggl.clust.s, main = "Agglomerative, single", cex = 0.25)

aggl.clust.a <- hclust(gower_df2, method = "average")
plot(aggl.clust.a, main = "Agglomerative, average", cex = 0.25)

aggl.clust.w <- hclust(gower_df2, method = "ward.D2")

png(
  "figures/agglomerative_hclust_ward_dendrogram_one_hot.png",
  width = 3000,
  height = 3000,
  res = 300
)
plot(aggl.clust.w, main = "Agglomerative, Ward", cex = 0.25)
dev.off()

## ------------ CLUSTERING STATS ------------

#source function
source("R/cstats.table.R")

# produce stats table, capping the maximum number of clusters at 10
stats.df.divisive <- cstats.table(gower_df2, divisive.clust, 10)
stats.df.divisive

stats.df.aggl.c <- cstats.table(gower_df2, aggl.clust.c, 10)
stats.df.aggl.c

#only 1 per cluster, not useful
stats.df.aggl.s <- cstats.table(gower_df2, aggl.clust.s, 10)
stats.df.aggl.s

stats.df.aggl.a <- cstats.table(gower_df2, aggl.clust.a, 10)
stats.df.aggl.a

stats.df.aggl.w <- cstats.table(gower_df2, aggl.clust.w, 10)
stats.df.aggl.w
write.csv(stats.df.aggl.w, "outputs/stats_hclust_ward_one_hot.csv")


## --------- Choosing the number of clusters - elbow withiness ---------

library(ggplot2)


# Agglomerative ward
ggplot(data = data.frame(t(cstats.table(
  gower_df2, aggl.clust.w, 10
))),
aes(x = cluster.number, y = within.cluster.ss)) +
  geom_point() +
  geom_line() +
  ggtitle("Agglomerative clustering, Ward") +
  labs(x = "Num.of clusters", y = "Within clusters sum of squares (SS)") +
  theme(plot.title = element_text(hjust = 0.5))
# k = 5
ggsave("figures/within_ss_ward_one_hot.png",
       width = 5,
       height = 5)

## --------- Choosing the number of clusters - silhouette ---------

#Agglomerative ward
ggplot(data = data.frame(t(cstats.table(
  gower_df2, aggl.clust.w, 15
))),
aes(x = cluster.number, y = avg.silwidth)) +
  geom_point() +
  geom_line() +
  ggtitle("Agglomerative clustering, Ward") +
  labs(x = "Num.of clusters", y = "Average silhouette width") +
  theme(plot.title = element_text(hjust = 0.5))
# k = 2
ggsave("figures/silwidth_ward_one_hot.png",
       width = 5,
       height = 5)

## --------- Plotting clusters on dendrogram  ---------

#load packages
library(reshape2)
library(purrr)
library(dplyr)
library(dendextend)
library(RColorBrewer)

brewer.pal(7, "Dark2")

# dendrogram with clusters divisive
#dendro <- as.dendrogram(divisive.clust)
#dendro.col <- dendro %>%
#  set("branches_k_color", k = 3, value = brewer.pal(7,"Dark2")) %>%
#  set("branches_lwd", 0.6) %>%
#  set("labels_colors",
#      value = c("darkslategray")) %>%
#  set("labels_cex", 0.5)
#ggd1 <- as.ggdend(dendro.col)
#ggplot(ggd1, theme = theme_minimal()) +
#  labs(x = "Num. observations", y = "Height", title = "Dendrogram divisive, k = 3")
#
## dendrogram with clusters agglomerative complete
#dendro <- as.dendrogram(aggl.clust.c)
#dendro.col <- dendro %>%
#  set("branches_k_color", k = 4, value = brewer.pal(7,"Dark2")) %>%
#  set("branches_lwd", 0.6) %>%
#  set("labels_colors",
#      value = c("darkslategray")) %>%
#  set("labels_cex", 0.5)
#ggd1 <- as.ggdend(dendro.col)
#ggplot(ggd1, theme = theme_minimal()) +
#  labs(x = "Num. observations", y = "Height", title = "Dendrogram agglomerative complete, k = 4")
#
## dendrogram with clusters agglomerative average
#dendro <- as.dendrogram(aggl.clust.a)
#dendro.col <- dendro %>%
#  set("branches_k_color", k = 6, value = brewer.pal(7,"Dark2")) %>%
#  set("branches_lwd", 0.6) %>%
#  set("labels_colors",
#      value = c("darkslategray")) %>%
#  set("labels_cex", 0.5)
#ggd1 <- as.ggdend(dendro.col)
#ggplot(ggd1, theme = theme_minimal()) +
#  labs(x = "Num. observations", y = "Height", title = "Dendrogram agglomerative average, k = 6")

# dendrogram with clusters agglomerative ward k = 3
dendro <- as.dendrogram(aggl.clust.w)
dendro.col <- dendro %>%
  set("branches_k_color",
      k = 3,
      value = brewer.pal(7, "Dark2")) %>%
  set("branches_lwd", 0.6) %>%
  set("labels_colors",
      value = c("darkslategray")) %>%
  set("labels_cex", 0.5)
ggd1 <- as.ggdend(dendro.col)
wk3 <- ggplot(ggd1, theme = theme_minimal()) +
  labs(x = "Num. observations", y = "Height", title = "Dendrogram agglomerative ward, k = 3")

# dendrogram with clusters agglomerative ward k = 6
dendro <- as.dendrogram(aggl.clust.w)
dendro.col <- dendro %>%
  set("branches_k_color",
      k = 6,
      value = brewer.pal(7, "Dark2")) %>%
  set("branches_lwd", 0.6) %>%
  set("labels_colors",
      value = c("darkslategray")) %>%
  set("labels_cex", 0.5)
ggd1 <- as.ggdend(dendro.col)
wk6 <- ggplot(ggd1, theme = theme_minimal()) +
  labs(x = "Num. observations", y = "Height", title = "Dendrogram agglomerative ward, k = 6")

library(patchwork)
wk3 + wk6
ggsave("figures/dendro_ward_k3_k6_one_hot.png",
       width = 20,
       height = 10)


## --------- PCOA scatterplot with cluster annotation ---------

#convert gower df to distance matrix
dataset_dist2 <- stats::as.dist(gower_df2)

library(vegan)
#run PCOA
#dataset_pcoa <- ape::pcoa(dataset_dist)
dataset_pcoa2 <- wcmdscale(d = dataset_dist2, eig = TRUE, add = "lingoes")
traitvectors12=scores(envfit(dataset_pcoa2, df2, na.rm = T, add = T),display="vectors")
traitvectors34=scores(envfit(dataset_pcoa2, df2, na.rm = T, add = T,choices=c(3,4)),display="vectors")
traitd=as.data.frame(traitvectors12)
traitd$Dim3=as.data.frame(traitvectors34)$Dim3
traitd$Dim4=as.data.frame(traitvectors34)$Dim4
traitd$trait=rownames(traitd)

speciesv=scores(dataset_pcoa2,display="species")
speciesd=as.data.frame(speciesv)

# use original traits (multistate) for visualization
dft<-readRDS(file = here::here("outputs/df_filt_trans_one_hot.rds"))
df_ord=as.data.frame(c(dft,speciesd))

#select number of clusters
# divisive.clust = 3
# aggl.clust.c (complete) = 4
# aggl.clust.a (average) = 2
# aggl.clust.w (ward) = 5
clust.num <- cutree(aggl.clust.w, k = 4)
df_clust=merge(clust.num,dft,by=0)
names(df_clust)[c(1,2)]=c("species","clust")
write.csv(file="outputs/species_hclust_k4_one_hot.csv",df_clust)

clust.num <- cutree(aggl.clust.w, k = 6)
df_clust=merge(clust.num,dft,by=0)
names(df_clust)[c(1,2)]=c("species","clust")
write.csv(file="outputs/species_hclust_k6_one_hot.csv",df_clust)

#clust.num[clust.num==1]="showy perennials"
#clust.num[clust.num==2]="dioecious trees"
#clust.num[clust.num==3]="small/short-lived herbs"
#clust.num[clust.num==4]="showy trees"
#clust.num[clust.num==5]="showy trees"
#clust.num[clust.num==6]="wind-pollinated sexually monomorphic"

library(ggrepel)

#minimal value for arrow length
minarrow=0.2

#plot points on first two axes, coloured by cluster
p1=ggplot() + geom_point(data=df_ord,aes(x=Dim1,y=Dim2,color = as.factor(clust.num)))+geom_segment(data=traitd[sqrt(traitd$Dim1^2+traitd$Dim2^2)>minarrow,],aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd[sqrt(traitd$Dim1^2+traitd$Dim2^2)>minarrow,],aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord, geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(clust.num)), alpha = 0.25)+theme_bw() + theme(legend.position="top",legend.title=element_blank()) 

#plot points on second and third axis, coloured by cluster
p2=ggplot() + geom_point(data=df_ord,aes(x=Dim2,y=Dim3,color = as.factor(clust.num)))+geom_segment(data=traitd[sqrt(traitd$Dim2^2+traitd$Dim3^2)>minarrow,],aes(x=0,y=0,xend=Dim2/2,yend=Dim3/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd[sqrt(traitd$Dim2^2+traitd$Dim3^2)>minarrow,],aes(x=Dim2/2,y=Dim3/2,label=trait))+  stat_ellipse(data=df_ord, geom = "polygon", aes(x=Dim2,y=Dim3,fill = as.factor(clust.num)), alpha = 0.25)+theme_bw() + theme(legend.position="none")

#plot points on third and fourth axis, coloured by cluster
p3=ggplot() + geom_point(data=df_ord,aes(x=Dim3,y=Dim4,color = as.factor(clust.num)))+geom_segment(data=traitd[sqrt(traitd$Dim3^2+traitd$Dim4^2)>minarrow,],aes(x=0,y=0,xend=Dim3/2,yend=Dim4/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd[sqrt(traitd$Dim3^2+traitd$Dim4^2)>minarrow,],aes(x=Dim3/2,y=Dim4/2,label=trait))+  stat_ellipse(data=df_ord, geom = "polygon", aes(x=Dim3,y=Dim4,fill = as.factor(clust.num)), alpha = 0.25)+theme_bw() + theme(legend.position="none")

library(patchwork)
p1 / p2 / p3

ggsave("figures/pcoa_hclust_k6_loadings_one_hot.png",
       width = 20,
       height = 40,
       units = 'cm')

clust.num <- cutree(aggl.clust.w, k = 3)

library(ggrepel)

#plot points on first two axes, coloured by cluster
p1=ggplot() + geom_point(data=df_ord,aes(x=Dim1,y=Dim2,color = as.factor(clust.num)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord, geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(clust.num)), alpha = 0.25) + theme(legend.position="none")

#plot points on second and third axis, coloured by cluster
p2=ggplot() + geom_point(data=df_ord,aes(x=Dim2,y=Dim3,color = as.factor(clust.num)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim2/2,yend=Dim3/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim2/2,y=Dim3/2,label=trait))+  stat_ellipse(data=df_ord, geom = "polygon", aes(x=Dim2,y=Dim3,fill = as.factor(clust.num)), alpha = 0.25) + theme(legend.position="none")

#plot points on third and fourth axis, coloured by cluster
p3=ggplot() + geom_point(data=df_ord,aes(x=Dim3,y=Dim4,color = as.factor(clust.num)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim3/2,yend=Dim4/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim3/2,y=Dim4/2,label=trait))+  stat_ellipse(data=df_ord, geom = "polygon", aes(x=Dim3,y=Dim4,fill = as.factor(clust.num)), alpha = 0.25) + theme(legend.position="none")

library(patchwork)
p1 / p2 / p3

ggsave("figures/pcoa_hclust_k3_loadings_one_hot.png",
       width = 20,
       height = 40,
       units = 'cm')


#3d

clust.num <- cutree(aggl.clust.w, k = 6)

library(rgl)
open3d()
get_colors <- function(groups, group.col = palette()){
  groups <- as.factor(groups)
  ngrps <- length(levels(groups))
  if(ngrps > length(group.col)) 
    group.col <- rep(group.col, ngrps)
  color <- group.col[as.numeric(groups)]
  names(color) <- as.vector(groups)
  return(color)
}

plot_loadings3d <- function(loadings,scale=0.5){
for(i in 1:nrow(loadings)){
arrow3d(p0=c(0,0,0),p1=c(loadings[i,"Dim1"]*scale,loadings[i,"Dim2"]*scale,loadings[i,"Dim3"]*scale),type="lines",barblen=0.01)
text3d(x=loadings[i,"Dim1"]*scale,y=loadings[i,"Dim2"]*scale,z=loadings[i,"Dim3"]*scale,loadings$trait[i],cex=0.5)
}
}

rgl.spheres(df_ord[,"Dim1"], df_ord[,"Dim2"], df_ord[,"Dim3"], r = 0.005,color = get_colors(clust.num)) 
plot_loadings3d(traitd)
clear3d()

#end 3d


