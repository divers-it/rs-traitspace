library(dplyr)

#load formatted data
df<-readRDS(file = here::here("outputs/df_filt.rds"))

#remove mating system
#df<-subset(df, select=-c(sexmorphs))

#numeric columns only
nums <- unlist(lapply(df, is.numeric))
facts <- unlist(lapply(df, is.factor))

df2<-cbind(df[ , nums],df[ , facts])

str(df2)

#centre and scale
#df_nums<-scale(df_nums)

#transform?

pdf("figures/proteus_trait_hists.pdf")
par(mfrow=c(3,3))
#look at hists
for(i in 1:6){
  hist(df2[,i],main=colnames(df2)[i])
}
dev.off()


pdf("figures/proteus_trait_hists_transformed.pdf")
par(mfrow=c(3,3))
#look at log10 hists
for(i in 1:6){
  hist(log(df2[,i]),main=colnames(df2)[i])
}
dev.off()

#do log transformations
for(i in 1:6){
  df2[,i]<-log(df2[,i])
  #df2[,i]<-scale(df2[,i]) #if we want to scale
}


par(mfrow=c(1,1))
#dissimilarity matrix calc - weights?
library(cluster)
gower_df <- daisy(df2,
                  metric = "gower" )

summary(gower_df)

#introduced NAs - need to solve better
#gower_df[is.na(gower_df)]<-0

#
library(factoextra)

jpeg("figures/proteus_fviz_dist.jpeg",width=1000,height=1000)

fviz_dist(dist.obj = gower_df,
          order = TRUE, show_labels = F)

dev.off()

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
pam_palm = pam(gower_df, diss = TRUE, k = 3)

library("gridExtra")

pdf("figures/proteus_medoids_table.pdf",height=3,width=17)
grid.table(df[pam_palm$medoids, c(1:10)])
dev.off()

#To dig deeper into the characteristics of each cluster, we find the summary stats.
pam_summary <- df %>%
  mutate(cluster = pam_palm$clustering) %>%
  group_by(cluster) %>%
  do(cluster_summary = summary(.))

#summary stats of cluster 1
pam_summary$cluster_summary[[1]]

#set palette
library(RColorBrewer)
palette(brewer.pal(6,"Dark2"))

#the t-SNE or the t-Distributed Stochastic Neighbor Embedding technique
library(Rtsne)
library(ggplot2)
tsne_object <- Rtsne(gower_df, is_distance = TRUE)
tsne_df <- tsne_object$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam_palm$clustering))


#add rownames (NOT SURE IF MATCH)
tsne_df$names<-rownames(df2)

#prune down species name to make readable
inds <- round ( runif(280, 1, length(tsne_df$names)) )
tsne_df$names[inds]<-NA
tsne_df$names[sample(seq_along(tsne_df$names), 280, replace = FALSE)] <- NA


#put in medoids
for(i in 1:length(pam_palm$medoids)){
  ind<-grep(pam_palm$medoids[i],rownames(df2))
  tsne_df$names[ind]<-rownames(df2)[ind]
}

ggplot(aes(x = X, y = Y), data = tsne_df) +
  geom_point(aes(color = tsne_df$cluster)) +
  geom_text(aes(label=names,color = as.factor(tsne_df$cluster)),hjust=0, vjust=0)

ggsave("figures/scatterplot_proteus_clusters.pdf")

clust_df<-data.frame(names(tsne_df$cluster),tsne_df$cluster,row.names=NULL)
colnames(clust_df)<-c("Species","Cluster")
clust_df<-clust_df[order(clust_df$Species),]
head(clust_df)

clust_df$Species<-gsub(" ","_",clust_df$Species)

write.csv(clust_df,"outputs/proteus_trait_clustering.csv")

save.image("outputs/proteus_clustering.Rdata")
