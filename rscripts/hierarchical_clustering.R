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
#not logging ovaries
for(i in c(1,2,3,4,6)){
  df2[,i]<-log(df2[,i])
  #df2[,i]<-scale(df2[,i]) #if we want to scale
}


par(mfrow=c(1,1))
#dissimilarity matrix calc - weights?
library(cluster)
gower_df <- daisy(df2,
                  metric = "gower" )

summary(gower_df)

# The main input for the code below is dissimilarity (distance matrix)

## ------------ DIVISIVE CLUSTERING ------------

library(cluster)
divisive.clust <- diana(as.matrix(gower_df), 
                        diss = TRUE, keep.diss = TRUE)

pdf("figures/dendrogram.pdf",width=15,height=10)
pltree(divisive.clust, main = "Divisive",cex=0.25)
dev.off()

## ------------ AGGLOMERATIVE CLUSTERING ------------

#running different methods to compare results

aggl.clust.c <- hclust(gower_df, method = "complete")
plot(aggl.clust.c, main = "Agglomerative, complete linkages",cex=0.25)

aggl.clust.s <- hclust(gower_df, method = "single")
plot(aggl.clust.s, main = "Agglomerative, single",cex=0.25)

aggl.clust.a <- hclust(gower_df, method = "average")
plot(aggl.clust.a, main = "Agglomerative, average",cex=0.25)

aggl.clust.m <- hclust(gower_df, method = "median")
plot(aggl.clust.m, main = "Agglomerative, median",cex=0.25)

aggl.clust.w <- hclust(gower_df, method = "ward.D2")
plot(aggl.clust.w, main = "Agglomerative, Ward",cex=0.25)

## ------------ CLUSTERING STATS ------------

#source function
source("R/cstats.table.R")

# produce stats table, capping the maximum number of clusters at 7
stats.df.divisive <- cstats.table(gower_df, divisive.clust, 7)
stats.df.divisive

stats.df.aggl.c <-cstats.table(gower_df, aggl.clust.c, 7)
stats.df.aggl.c

#only 1 per cluster, not useful
#stats.df.aggl.s <-cstats.table(gower_df, aggl.clust.s, 7) 
#stats.df.aggl.s

stats.df.aggl.a <-cstats.table(gower_df, aggl.clust.a, 7) 
stats.df.aggl.a

#very large first cluster, not useful
#stats.df.aggl.m <-cstats.table(gower_df, aggl.clust.m, 7) 
#stats.df.aggl.m

stats.df.aggl.w <-cstats.table(gower_df, aggl.clust.w, 7) 
stats.df.aggl.w

## --------- Choosing the number of clusters - elbow withiness ---------

library(ggplot2)

#It shows how the within sum of squares — as a measure of closeness of observations : the lower it is the closer the observations within the clusters are — changes for the different number of clusters. Ideally, we should see a distinctive “bend” in the elbow where splitting clusters further gives only minor decrease in the SS.

# Divisive clustering
ggplot(data = data.frame(t(cstats.table(gower_df, divisive.clust, 15))), 
       aes(x=cluster.number, y=within.cluster.ss)) + 
  geom_point()+
  geom_line()+
  ggtitle("Divisive clustering") +
  labs(x = "Num.of clusters", y = "Within clusters sum of squares (SS)") +
  theme(plot.title = element_text(hjust = 0.5))
# k = 3

# Agglomerative complete
ggplot(data = data.frame(t(cstats.table(gower_df, aggl.clust.c, 15))), 
       aes(x=cluster.number, y=within.cluster.ss)) + 
  geom_point()+
  geom_line()+
  ggtitle("Agglomerative clustering, complete") +
  labs(x = "Num.of clusters", y = "Within clusters sum of squares (SS)") +
  theme(plot.title = element_text(hjust = 0.5))
# k = 4

# Agglomerative average
ggplot(data = data.frame(t(cstats.table(gower_df, aggl.clust.a, 15))), 
       aes(x=cluster.number, y=within.cluster.ss)) + 
  geom_point()+
  geom_line()+
  ggtitle("Agglomerative clustering, average") +
  labs(x = "Num.of clusters", y = "Within clusters sum of squares (SS)") +
  theme(plot.title = element_text(hjust = 0.5))
# k = 2

# Agglomerative ward
ggplot(data = data.frame(t(cstats.table(gower_df, aggl.clust.w, 15))), 
       aes(x=cluster.number, y=within.cluster.ss)) + 
  geom_point()+
  geom_line()+
  ggtitle("Agglomerative clustering, Ward") +
  labs(x = "Num.of clusters", y = "Within clusters sum of squares (SS)") +
  theme(plot.title = element_text(hjust = 0.5))
# k = 5

## --------- Choosing the number of clusters - silhouette ---------

#When it comes to silhouette assessment, the rule is you should choose the number that maximizes the silhouette coefficient because you want clusters that are distinctive (far) enough to be considered separate.
#The silhouette coefficient ranges between -1 and 1, with 1 indicating good consistency within clusters, -1 — not so good.

# Seems to be biased towards two clusters, I prefer elbow

# Divisive clustering
#ggplot(data = data.frame(t(cstats.table(gower_df, divisive.clust, 15))), 
#       aes(x=cluster.number, y=avg.silwidth)) + 
#  geom_point()+
#  geom_line()+
#  ggtitle("Divisive clustering") +
#  labs(x = "Num.of clusters", y = "Average silhouette width") +
#  theme(plot.title = element_text(hjust = 0.5))
## k = 2
#
## Agglomerative complete
#ggplot(data = data.frame(t(cstats.table(gower_df, aggl.clust.c, 15))), 
#       aes(x=cluster.number, y=avg.silwidth)) + 
#  geom_point()+
#  geom_line()+
#  ggtitle("Agglomerative clustering, complete") +
#  labs(x = "Num.of clusters", y = "Average silhouette width") +
#  theme(plot.title = element_text(hjust = 0.5))
## k = 2
#
#
## Agglomerative average
#ggplot(data = data.frame(t(cstats.table(gower_df, aggl.clust.a, 15))), 
#       aes(x=cluster.number, y=avg.silwidth)) + 
#  geom_point()+
#  geom_line()+
#  ggtitle("Agglomerative clustering, average") +
#  labs(x = "Num.of clusters", y = "Average silhouette width") +
#  theme(plot.title = element_text(hjust = 0.5))
## k = 2
#
## Agglomerative ward
#ggplot(data = data.frame(t(cstats.table(gower_df, aggl.clust.w, 15))), 
#       aes(x=cluster.number, y=avg.silwidth)) + 
#  geom_point()+
#  geom_line()+
#  ggtitle("Agglomerative clustering, Ward") +
#  labs(x = "Num.of clusters", y = "Average silhouette width") +
#  theme(plot.title = element_text(hjust = 0.5))
## k = 2

## --------- Plotting clusters on dendrogram  ---------

#load packages
library(reshape2)
library(purrr)
library(dplyr)
library(dendextend)
library(RColorBrewer)

brewer.pal(7,"Dark2")

# dendrogram with clusters divisive
dendro <- as.dendrogram(divisive.clust)
dendro.col <- dendro %>%
  set("branches_k_color", k = 3, value = brewer.pal(7,"Dark2")) %>%
  set("branches_lwd", 0.6) %>%
  set("labels_colors", 
      value = c("darkslategray")) %>% 
  set("labels_cex", 0.5)
ggd1 <- as.ggdend(dendro.col)
ggplot(ggd1, theme = theme_minimal()) +
  labs(x = "Num. observations", y = "Height", title = "Dendrogram divisive, k = 3")

# dendrogram with clusters agglomerative complete
dendro <- as.dendrogram(aggl.clust.c)
dendro.col <- dendro %>%
  set("branches_k_color", k = 4, value = brewer.pal(7,"Dark2")) %>%
  set("branches_lwd", 0.6) %>%
  set("labels_colors", 
      value = c("darkslategray")) %>% 
  set("labels_cex", 0.5)
ggd1 <- as.ggdend(dendro.col)
ggplot(ggd1, theme = theme_minimal()) +
  labs(x = "Num. observations", y = "Height", title = "Dendrogram agglomerative complete, k = 4")

# dendrogram with clusters agglomerative average
dendro <- as.dendrogram(aggl.clust.a)
dendro.col <- dendro %>%
  set("branches_k_color", k = 6, value = brewer.pal(7,"Dark2")) %>%
  set("branches_lwd", 0.6) %>%
  set("labels_colors", 
      value = c("darkslategray")) %>% 
  set("labels_cex", 0.5)
ggd1 <- as.ggdend(dendro.col)
ggplot(ggd1, theme = theme_minimal()) +
  labs(x = "Num. observations", y = "Height", title = "Dendrogram agglomerative average, k = 6")

# dendrogram with clusters agglomerative ward
dendro <- as.dendrogram(aggl.clust.w)
dendro.col <- dendro %>%
  set("branches_k_color", k = 5, value = brewer.pal(7,"Dark2")) %>%
  set("branches_lwd", 0.6) %>%
  set("labels_colors", 
      value = c("darkslategray")) %>% 
  set("labels_cex", 0.5)
ggd1 <- as.ggdend(dendro.col)
ggplot(ggd1, theme = theme_minimal()) +
  labs(x = "Num. observations", y = "Height", title = "Dendrogram agglomerative ward, k = 5")

## --------- Heatmap of cluster properties (categorical only) ---------

#factors only
facts <- unlist(lapply(df2, is.factor))
df3<-df2[ , facts]

#add column with species names
df3$species<-rownames(df2)

#cut tree to get cluster numbers, add to df
clust.num <- cutree(aggl.clust.w, k = 5)
df3.cl <- cbind(df3, clust.num)

# factors have to be converted to characters in order not to be dropped
clust.long <- melt(data.frame(lapply(df3.cl, as.character), stringsAsFactors=FALSE), 
                   id = c("species", "clust.num"), factorsAsStrings=T)

#get counts per state per cluster - can do heatmap with this
clust.long.q <- clust.long %>%
  group_by(clust.num, variable, value) %>%
  mutate(count = n_distinct(species)) %>%
  distinct(clust.num, variable, value, count)

#remove NAs
clust.long.q<-na.omit(clust.long.q)

# calculating the percent of each factor level in the absolute count of cluster members
clust.long.p <- clust.long.q %>%
  group_by(clust.num, variable) %>%
  mutate(perc = count / sum(count)) %>%
  arrange(clust.num)

#plot heatmap
#deeper blue corresponds to a higher relative number of observations within a cluster

heatmap.p <- ggplot(clust.long.p, aes(x = clust.num, y = factor(value, levels = c("woody",
                                                                                  "herbaceous",
                                                                                  "non-climbing",
                                                                                  "climbing",
                                                                                  "monomorphic",
                                                                                  "dimorphic",
                                                                                  "polymorphic",
                                                                                  "short",
                                                                                  "long",
                                                                                  "selfing",
                                                                                  "outcrossing",
                                                                                  "mixed",
                                                                                  "abiotic",
                                                                                  "biotic",
                                                                                  "autonomous",
                                                                                  "bisexual",
                                                                                  "unisexual",
                                                                                  "superior",
                                                                                  "inferior",
                                                                                  "intermediate",
                                                                                  "none",
                                                                                  "pollen",
                                                                                  "nectar",
                                                                                  "oil",
                                                                                  "nursery",
                                                                                  "foodbodies",
                                                                                  "perfume",
                                                                                  "heat",
                                                                                  "actinomorphic",
                                                                                  "zygomorphic",
                                                                                  "other",
                                                                                  "bright",
                                                                                  "green/brown",
                                                                                  "whitish"), ordered = T))) +
  geom_tile(aes(fill = perc), alpha = 0.85)+
  labs(title = "Distribution of characteristics across clusters (percentage)", x = "Cluster number", y = NULL) +
  geom_hline(yintercept = 2.5) + 
  geom_hline(yintercept = 4.5) + 
  geom_hline(yintercept = 7.5) + 
  geom_hline(yintercept = 9.5) + 
  geom_hline(yintercept = 12.5) + 
  geom_hline(yintercept = 14.5) + 
  geom_hline(yintercept = 17.5) + 
  geom_hline(yintercept = 20.5) + 
  geom_hline(yintercept = 26.5) + 
  geom_hline(yintercept = 29.5) + 
  scale_fill_gradient2(low = "darkslategray1", mid = "yellow", high =  "turquoise4")

heatmap.p

## --------- Values of continuous variables per cluster ---------
library(viridis)

clust_df<-cbind(df2,clust.num)

colnames(clust_df)[length(colnames(clust_df))]<-"cluster"

#change to character for plotting
clust_df$cluster<-as.character(clust_df$cluster)

# fertile stamens
ggplot(clust_df, aes(
  x = cluster,
  y = Numberoffertilestamens,
  fill = cluster
)) + geom_boxplot(outlier.shape = NA, width = 0.2) +
  geom_violin(width=0.75) +
  geom_jitter(
    aes(color = cluster),
    size = 1,
    alpha = 0.5,
    width = 0.05,
    height = 0
  ) +
  theme(legend.position = "none",
        axis.title.x = element_text(vjust = 0.3)) +
  scale_y_continuous(name="Number of fertile stamens") +
  scale_x_discrete(name = "cluster") +
  scale_fill_viridis(discrete = TRUE, alpha = 0.3) +
  scale_color_viridis(discrete = TRUE, alpha = 0.8)

# ovules per functional carpel
ggplot(clust_df, aes(
  x = cluster,
  y = Numberofovulesperfunctionalcarpel,
  fill = cluster
)) + geom_boxplot(outlier.shape = NA, width = 0.2) +
  geom_violin(width=0.75) +
  geom_jitter(
    aes(color = cluster),
    size = 1,
    alpha = 0.5,
    width = 0.05,
    height = 0
  ) +
  theme(legend.position = "none",
        axis.title.x = element_text(vjust = 0.3)) +
  scale_y_continuous(name="Number of ovules per functional carpel") +
  scale_x_discrete(name = "cluster") +
  scale_fill_viridis(discrete = TRUE, alpha = 0.3) +
  scale_color_viridis(discrete = TRUE, alpha = 0.8)

# Flower diameter
ggplot(clust_df, aes(
  x = cluster,
  y = Flowerdiameter,
  fill = cluster
)) + geom_boxplot(outlier.shape = NA, width = 0.2) +
  #geom_violin(width=0.75) +
  geom_jitter(
    aes(color = cluster),
    size = 1,
    alpha = 0.5,
    width = 0.05,
    height = 0
  ) +
  theme(legend.position = "none",
        axis.title.x = element_text(vjust = 0.3)) +
  scale_y_continuous(name="Flower diameter") +
  scale_x_discrete(name = "cluster") +
  scale_fill_viridis(discrete = TRUE, alpha = 0.3) +
  scale_color_viridis(discrete = TRUE, alpha = 0.8)

# Structural carpels
ggplot(clust_df, aes(
  x = cluster,
  y = Numberofstructuralcarpels,
  fill = cluster
)) + geom_boxplot(outlier.shape = NA, width = 0.2) +
  #geom_violin(width=0.75) +
  geom_jitter(
    aes(color = cluster),
    size = 1,
    alpha = 0.5,
    width = 0.05,
    height = 0
  ) +
  theme(legend.position = "none",
        axis.title.x = element_text(vjust = 0.3)) +
  scale_y_continuous(name="Number of structural carpels") +
  scale_x_discrete(name = "cluster") +
  scale_fill_viridis(discrete = TRUE, alpha = 0.3) +
  scale_color_viridis(discrete = TRUE, alpha = 0.8)

# Fusion of ovaries
ggplot(clust_df, aes(
  x = cluster,
  y = Fusionofovaries,
  fill = cluster
)) + geom_boxplot(outlier.shape = NA, width = 0.2) +
  #geom_violin(width=0.75) +
  geom_jitter(
    aes(color = cluster),
    size = 1,
    alpha = 0.5,
    width = 0.05,
    height = 0
  ) +
  theme(legend.position = "none",
        axis.title.x = element_text(vjust = 0.3)) +
  scale_y_continuous(name="Fusion of ovaries") +
  scale_x_discrete(name = "cluster") +
  scale_fill_viridis(discrete = TRUE, alpha = 0.3) +
  scale_color_viridis(discrete = TRUE, alpha = 0.8)

# Maximum vertical height
ggplot(clust_df, aes(
  x = cluster,
  y = Maximumverticalheight,
  fill = cluster
)) + geom_boxplot(outlier.shape = NA, width = 0.2) +
  #geom_violin(width=0.75) +
  geom_jitter(
    aes(color = cluster),
    size = 1,
    alpha = 0.5,
    width = 0.05,
    height = 0
  ) +
  theme(legend.position = "none",
        axis.title.x = element_text(vjust = 0.3)) +
  scale_y_continuous(name="Maximum vertical height") +
  scale_x_discrete(name = "cluster") +
  scale_fill_viridis(discrete = TRUE, alpha = 0.3) +
  scale_color_viridis(discrete = TRUE, alpha = 0.8)


## --------- PCOA scatterplot with cluster annotation ---------

#select number of clusters
# divisive.clust = 3 
# aggl.clust.c (complete) = 4
# aggl.clust.a (average) = 2
# aggl.clust.w (ward) = 5
clust.num <- cutree(aggl.clust.a, k = 6)

#convert to distance matrix
dataset_dist <- stats::as.dist(gower_df)

#run PCOA
dataset_pcoa <- ape::pcoa(dataset_dist)

#plot points on first two axes, coloured by cluster
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2, color = as.factor(clust.num))) +
  geom_point() +
  stat_ellipse(geom = "polygon",
               aes(fill = as.factor(clust.num)), 
               alpha = 0.25) +
  xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))

