library(dplyr)
library(gridExtra)
library(Rtsne)
library(ggplot2)

#load formatted data
df<-readRDS(file = here::here("outputs/df_filt_trans.rds"))

library(cluster)
gower_df <- daisy(df,
                  metric = "gower" )

summary(gower_df)

#introduced NAs - need to solve better
#gower_df[is.na(gower_df)]<-0

#Silhouette Width to select the optimal number of clusters
#The silhouette width is one of the very popular choices when it comes to selecting the optimal number of clusters. It measures the similarity of each point to its cluster, and compares that to the similarity of the point with the closest neighboring cluster. This metric ranges between -1 to 1, where a higher value implies better similarity of the points to their clusters. Therefore, a higher value of the Silhouette Width is desirable. We calculate this metric for a range of cluster numbers and find where it is maximized. The following code shows the implementation in R:
silhouette <- c()
silhouette = c(silhouette, NA)
for(i in 2:10){
  pam_clusters = pam(as.matrix(gower_df),
                     diss = TRUE,
                     k = i)
  silhouette = c(silhouette ,pam_clusters$silinfo$avg.width)
}

plot(1:10, silhouette,
     xlab = "Clusters",
     ylab = "Silhouette Width")
lines(1:10, silhouette)

#construct a PAM model with 2 clusters, and try to interpret the behavior of these clusters with the help of the medoids.
pam.gower = pam(gower_df, diss = TRUE, k = 2)

pdf("figures/pam_medoids_table.pdf",height=3,width=17)
grid.table(df[pam.gower$medoids, c(1:10)])
dev.off()

#To dig deeper into the characteristics of each cluster, we find the summary stats.
pam_summary <- df %>%
  mutate(cluster = pam.gower$clustering) %>%
  group_by(cluster) %>%
  do(cluster_summary = summary(.))

#summary stats of cluster 1
pam_summary$cluster_summary[[1]]

#set palette
library(RColorBrewer)
palette(brewer.pal(6,"Dark2"))

#the t-SNE or the t-Distributed Stochastic Neighbor Embedding technique
#alternative to PCOA
tsne_object <- Rtsne(gower_df, is_distance = TRUE)
tsne_df <- tsne_object$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam.gower$clustering))

#add rownames (NOT SURE IF MATCH)
tsne_df$names<-rownames(df2)

#prune down species name to make readable
inds <- round ( runif(280, 1, length(tsne_df$names)) )
tsne_df$names[inds]<-NA
tsne_df$names[sample(seq_along(tsne_df$names), 280, replace = FALSE)] <- NA

#put in medoids
for(i in 1:length(pam.gower$medoids)){
  ind<-grep(pam.gower$medoids[i],rownames(df2))
  tsne_df$names[ind]<-rownames(df2)[ind]
}

#plot points on first two axes, coloured by cluster
ggplot(tsne_df, aes(x = X, y = Y, fill = as.factor(cluster))) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.5,
    size=3,
    stroke = 0.5
  ) + 
  stat_ellipse(geom = "polygon",
               aes(fill = tsne_df$cluster), 
               alpha = 0.25) +
  xlab("t-SNE Axis 1") +
  ylab("t-SNE Axis 2")

ggsave("figures/scatterplot_pam_clusters.pdf")

clust_df<-data.frame(names(tsne_df$cluster),tsne_df$cluster,row.names=NULL)
colnames(clust_df)<-c("Species","Cluster")
clust_df<-clust_df[order(clust_df$Species),]
head(clust_df)

clust_df$Species<-gsub(" ","_",clust_df$Species)

write.csv(clust_df,"outputs/pam_clustering.csv")

#save.image("outputs/pam_clustering.Rdata")

#pcoa
dataset_dist <- stats::as.dist(gower_df)
dataset_pcoa <- ape::pcoa(dataset_dist)

#plot points on first two axes, coloured by cluster
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2, fill = as.factor(pam.gower$clustering))) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.5,
    size=3,
    stroke = 0.5
  ) + 
  stat_ellipse(geom = "polygon",
               aes(fill =  as.factor(pam.gower$clustering)), 
               alpha = 0.25) +
  xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))
