rm(list = ls())

#load packages
library(dplyr)
library(cluster)
library(ggplot2)
library(reshape2)
library(purrr)
library(dplyr)
library(dendextend)
library(RColorBrewer)
library(patchwork)
library(ggrepel)
library(networkD3)
library(gridExtra)

#load formatted data
df <- readRDS(file = here::here("outputs/one_hot_6_df_filt_trans.rds"))

# dissimilarity matrix
gower_df <- daisy(df,
                  metric = "gower")
summary(gower_df)

### Divisive clustering ----
divisive.clust <- diana(as.matrix(gower_df),
                        diss = TRUE, keep.diss = TRUE)

pltree(divisive.clust, main = "Divisive", cex = 0.25)

### Agglomerative clustering ----

# running different methods to compare results
aggl.clust.c <- hclust(gower_df, method = "complete")
plot(aggl.clust.c, main = "Agglomerative, complete linkages", cex = 0.25)

aggl.clust.a <- hclust(gower_df, method = "average")
plot(aggl.clust.a, main = "Agglomerative, average", cex = 0.25)

aggl.clust.w <- hclust(gower_df, method = "ward.D2")
plot(aggl.clust.w, main = "Agglomerative, Ward", cex = 0.25)

### Clustering stats ----

# source function
source("R/cstats.table.R")

# produce stats table, capping the maximum number of clusters at 10

# divisive
stats.df.d <- cstats.table(gower_df, divisive.clust, 10)
stats.df.d

# complete
stats.df.aggl.c <- cstats.table(gower_df, aggl.clust.c, 10)
stats.df.aggl.c

# average
stats.df.aggl.a <- cstats.table(gower_df, aggl.clust.a, 10)
stats.df.aggl.a

# ward
stats.df.aggl.w <- cstats.table(gower_df, aggl.clust.w, 10)
stats.df.aggl.w

# stats of one clustering approach / K value
cluster.stats(d = gower_df, clustering = cutree(aggl.clust.w, 4))

### Choosing k - Elbow withinness ----

# It shows how the within sum of squares — as a measure of closeness of observations : the lower it is the closer the observations within the clusters are — changes for the different number of clusters. Ideally, we should see a distinctive “bend” in the elbow where splitting clusters further gives only minor decrease in the SS.

# Agglomerative ward
ggplot(data = data.frame(t(cstats.table(
  gower_df, aggl.clust.w, 10
))),
aes(x = cluster.number, y = within.cluster.ss)) +
  geom_point() +
  geom_line() +
  ggtitle("Agglomerative clustering, Ward") +
  labs(x = "Num.of clusters", y = "Within clusters sum of squares (SS)") +
  theme(plot.title = element_text(hjust = 0.5))
# k = 5

### Choosing k - Silhouette ----

# When it comes to silhouette assessment, the rule is you should choose the number that maximizes the silhouette coefficient because you want clusters that are distinctive (far) enough to be considered separate. The silhouette coefficient ranges between -1 and 1, with 1 indicating good consistency within clusters, -1 — not so good.
# Seems to be biased towards two clusters, I prefer elbow

# Agglomerative ward
ggplot(data = data.frame(t(cstats.table(
  gower_df, aggl.clust.w, 15
))),
aes(x = cluster.number, y = avg.silwidth)) +
  geom_point() +
  geom_line() +
  ggtitle("Agglomerative clustering, Ward") +
  labs(x = "Num.of clusters", y = "Average silhouette width") +
  theme(plot.title = element_text(hjust = 0.5))
# k = 3

### Plotting clusters on dendrogram ----

# set palette
brewer.pal(7, "Dark2")

# dendrogram with clusters agglomerative ward k = 3
dendro <- as.dendrogram(aggl.clust.w)
dendro.col <- dendro %>%
  dendextend::set("branches_k_color",
                  k = 3,
                  value = c("#4DAF4A","#377EB8","#E41A1C")) %>%
  dendextend::set("branches_lwd", 0.6) %>%
  dendextend::set("labels_colors",
                  value = c("darkslategray")) %>%
  dendextend::set("labels_cex", 0.5)
ggd1 <- as.ggdend(dendro.col)
wk3 <- ggplot(ggd1, theme = theme_minimal()) +
  labs(x = "Num. observations", y = "Height", title = "Dendrogram agglomerative ward, k = 3")

# dendrogram with clusters agglomerative ward k = 5
dendro <- as.dendrogram(aggl.clust.w)
dendro.col <- dendro %>%
  dendextend::set("branches_k_color",
                  k = 5,
                  value = brewer.pal(7, "Dark2")) %>%
  dendextend::set("branches_lwd", 0.6) %>%
  dendextend::set("labels_colors",
                  value = c("darkslategray")) %>%
  dendextend::set("labels_cex", 0.5)
ggd1 <- as.ggdend(dendro.col)
wk5 <- ggplot(ggd1, theme = theme_minimal()) +
  labs(x = "Num. observations", y = "Height", title = "Dendrogram agglomerative ward, k = 5")

# plot coloured dendrograms
wk3 + wk5

# ggsave("figures/one_hot_10.2_dendrograms_ward2_coloured_by_cluster.png",
#        width = 10,
#        height = 10)

### PCOA scatterplot with cluster annotation ----

#convert gower df to distance matrix
dataset_dist <- stats::as.dist(gower_df)

#run PCOA
dataset_pcoa <- ape::pcoa(dataset_dist)

#Recalculate relative eigenvalues by removing negative eigenvalues as in Mouillot et al.  
ev_pcoa <- dataset_pcoa$values$Eigenvalues
ev_pcoa_g0 <- ev_pcoa[ev_pcoa>0]
rel_ev_pcoa_g0 <- ev_pcoa_g0/sum(ev_pcoa_g0)

# cut dendrogram to get three clusters
clust.num <- cutree(aggl.clust.w, k = 5)

# colors
cols<-harrypotter::hp(5,option="ronweasley2")

# plot points on first two axes, coloured by cluster with species names
ggplot(data.frame(dataset_pcoa$vectors),
       aes(
         x = Axis.1,
         y = Axis.2,
         fill = as.factor(clust.num)
       )) +
  geom_point(
    color = "black",
    shape = 21,
    alpha = 0.5,
    size = 3,
    stroke = 0.5
  ) +  stat_ellipse(geom = "polygon", aes(fill = as.factor(clust.num)), alpha = 0.25) +
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = cols) +
  xlab(paste("Axis 1: relative eigenvalue =", round(rel_ev_pcoa_g0[1], 2))) +
  ylab(paste("Axis 2: relative eigenvalue =", round(rel_ev_pcoa_g0[2], 2))) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  )

# table of different k values (2-7)
for (i in 2:7) {
  if (i == 2) {
    clust.num.k.2.7 <- paste("k",i,"cluster",as.character(cutree(aggl.clust.w, k = i)),sep="_")
  } else {
    clust.num.k.2.7 <-
      cbind(clust.num.k.2.7, paste("k",i,"cluster",as.character(cutree(aggl.clust.w, k = i)),sep="_"))
  }
}

colnames(clust.num.k.2.7)<-c("2clusters",
                             "3clusters",
                             "4clusters",
                             "5clusters",
                             "6clusters",
                             "7clusters")
clust.num.k.2.7.df <-as.data.frame(clust.num.k.2.7)
rownames(clust.num.k.2.7.df)<-names(cutree(aggl.clust.w, k = 2))

# save output for downstream use
saveRDS(clust.num.k.2.7.df, file = here::here("outputs/one_hot_10.2_clust_num_k_2_7_ward.rds"))
