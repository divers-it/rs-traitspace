rm(list = ls())
library(dplyr)

#load formatted data
df <- readRDS(file = here::here("outputs/df_filt_trans.rds"))

#dissimilarity matrix
library(cluster)
gower_df <- daisy(df,
                  metric = "gower")

summary(gower_df)

# The main input for the code below is dissimilarity (distance matrix)

## ------------ DIVISIVE CLUSTERING ------------

library(cluster)
divisive.clust <- diana(as.matrix(gower_df),
                        diss = TRUE, keep.diss = TRUE)

pltree(divisive.clust, main = "Divisive", cex = 0.25)

## ------------ AGGLOMERATIVE CLUSTERING ------------

#running different methods to compare results

aggl.clust.c <- hclust(gower_df, method = "complete")
plot(aggl.clust.c, main = "Agglomerative, complete linkages", cex = 0.25)

aggl.clust.s <- hclust(gower_df, method = "single")
plot(aggl.clust.s, main = "Agglomerative, single", cex = 0.25)

aggl.clust.a <- hclust(gower_df, method = "average")
plot(aggl.clust.a, main = "Agglomerative, average", cex = 0.25)

aggl.clust.w <- hclust(gower_df, method = "ward.D2")

png(
  "figures/agglomerative_hclust_ward_dendrogram.png",
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
stats.df.divisive <- cstats.table(gower_df, divisive.clust, 10)
stats.df.divisive

stats.df.aggl.c <- cstats.table(gower_df, aggl.clust.c, 10)
stats.df.aggl.c

#only 1 per cluster, not useful
stats.df.aggl.s <- cstats.table(gower_df, aggl.clust.s, 10)
stats.df.aggl.s

stats.df.aggl.a <- cstats.table(gower_df, aggl.clust.a, 10)
stats.df.aggl.a

stats.df.aggl.w <- cstats.table(gower_df, aggl.clust.w, 10)
stats.df.aggl.w
write.csv(stats.df.aggl.w, "outputs/stats_hclust_ward.csv")

#stats of one clustering approach / K value
cluster.stats(d = gower_df, clustering = cutree(aggl.clust.w, 4))


## --------- Choosing the number of clusters - elbow withiness ---------

library(ggplot2)

#It shows how the within sum of squares — as a measure of closeness of observations : the lower it is the closer the observations within the clusters are — changes for the different number of clusters. Ideally, we should see a distinctive “bend” in the elbow where splitting clusters further gives only minor decrease in the SS.

# Divisive clustering
#ggplot(data = data.frame(t(cstats.table(gower_df, divisive.clust, 10))),
#       aes(x=cluster.number, y=within.cluster.ss)) +
#  geom_point()+
#  geom_line()+
#  ggtitle("Divisive clustering") +
#  labs(x = "Num.of clusters", y = "Within clusters sum of squares (SS)") +
#  theme(plot.title = element_text(hjust = 0.5))
## k = 3
#
## Agglomerative complete
#ggplot(data = data.frame(t(cstats.table(gower_df, aggl.clust.c, 10))),
#       aes(x=cluster.number, y=within.cluster.ss)) +
#  geom_point()+
#  geom_line()+
#  ggtitle("Agglomerative clustering, complete") +
#  labs(x = "Num.of clusters", y = "Within clusters sum of squares (SS)") +
#  theme(plot.title = element_text(hjust = 0.5))
## k = 4
#
## Agglomerative average
#ggplot(data = data.frame(t(cstats.table(gower_df, aggl.clust.a, 10))),
#       aes(x=cluster.number, y=within.cluster.ss)) +
#  geom_point()+
#  geom_line()+
#  ggtitle("Agglomerative clustering, average") +
#  labs(x = "Num.of clusters", y = "Within clusters sum of squares (SS)") +
#  theme(plot.title = element_text(hjust = 0.5))
## k = 2

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
ggsave("figures/within_ss_ward.png",
       width = 5,
       height = 5)

## --------- Choosing the number of clusters - silhouette ---------

#When it comes to silhouette assessment, the rule is you should choose the number that maximizes the silhouette coefficient because you want clusters that are distinctive (far) enough to be considered separate. The silhouette coefficient ranges between -1 and 1, with 1 indicating good consistency within clusters, -1 — not so good.

#Seems to be biased towards two clusters, I prefer elbow

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

#Agglomerative ward
ggplot(data = data.frame(t(cstats.table(
  gower_df, aggl.clust.w, 15
))),
aes(x = cluster.number, y = avg.silwidth)) +
  geom_point() +
  geom_line() +
  ggtitle("Agglomerative clustering, Ward") +
  labs(x = "Num.of clusters", y = "Average silhouette width") +
  theme(plot.title = element_text(hjust = 0.5))
# k = 2
ggsave("figures/silwidth_ward.png",
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
ggsave("figures/dendro_ward_k3_k6.png",
       width = 20,
       height = 10)

## --------- Heatmap of cluster properties (categorical only) ---------

#factors only
facts <- unlist(lapply(df, is.factor))
df2 <- df[, facts]

#add column with species names
df2$species <- rownames(df)

#cut tree to get cluster numbers, add to df
clust.num <- cutree(aggl.clust.w, k = 3)
df2.cl <- cbind(df2, clust.num)

# factors have to be converted to characters in order not to be dropped
clust.long <-
  melt(
    data.frame(lapply(df2.cl, as.character), stringsAsFactors = FALSE),
    id = c("species", "clust.num"),
    factorsAsStrings = T
  )

#get counts per state per cluster - can do heatmap with this
clust.long.q <- clust.long %>%
  group_by(clust.num, variable, value) %>%
  mutate(count = n_distinct(species)) %>%
  distinct(clust.num, variable, value, count)

#remove NAs
clust.long.q <- na.omit(clust.long.q)

# calculating the percent of each factor level in the absolute count of cluster members
clust.long.p <- clust.long.q %>%
  group_by(clust.num, variable) %>%
  mutate(perc = count / sum(count)) %>%
  arrange(clust.num)

#rename dispersal and pollination values as they are the same, as well as "other"

clust.long.p[clust.long.p$variable == "Pollination", ]$value = paste(clust.long.p[clust.long.p$variable ==
                                                                                    "Pollination", ]$variable, clust.long.p[clust.long.p$variable == "Pollination", ]$value, sep =
                                                                       ".")
clust.long.p[clust.long.p$variable == "Dispersal", ]$value = paste(clust.long.p[clust.long.p$variable ==
                                                                                  "Dispersal", ]$variable, clust.long.p[clust.long.p$variable == "Dispersal", ]$value, sep =
                                                                     ".")
clust.long.p[clust.long.p$variable == "FloralReward", ]$value = paste(clust.long.p[clust.long.p$variable ==
                                                                                     "FloralReward", ]$variable, clust.long.p[clust.long.p$variable == "FloralReward", ]$value, sep =
                                                                        ".")
clust.long.p[clust.long.p$variable == "FlowerSymmetry", ]$value = paste(clust.long.p[clust.long.p$variable ==
                                                                                       "FlowerSymmetry", ]$variable, clust.long.p[clust.long.p$variable == "FlowerSymmetry", ]$value, sep =
                                                                          ".")

#plot heatmap
#deeper blue corresponds to a higher relative number of observations within a cluster

heatmap.p <-
  ggplot(clust.long.p, aes(x = clust.num, y = factor(
    value,
    levels = c(
      "woody",
      "herbaceous",
      "non-climbing",
      "climbing",
      "short",
      "long",
      "selfing",
      "outcrossing",
      "mixed",
      "Pollination.abiotic",
      "Pollination.biotic",
      "Pollination.autonomous",
      "Dispersal.abiotic",
      "Dispersal.biotic",
      "Dispersal.autonomous",
      "bisexual",
      "unisexual",
      "monomorphic",
      "dimorphic",
      "polymorphic",
      "superior",
      "inferior",
      "intermediate",
      "FloralReward.none",
      "FloralReward.pollen",
      "FloralReward.nectar",
      "FloralReward.oil",
      "FloralReward.nursery",
      "FloralReward.foodbodies",
      "FloralReward.perfume",
      "FloralReward.heat",
      "FloralReward.other",
      "FlowerSymmetry.actinomorphic",
      "FlowerSymmetry.zygomorphic",
      "FlowerSymmetry.other",
      "bright",
      "dull",
      "whitish"
    ),
    ordered = T
  ))) +
  geom_tile(aes(fill = perc), alpha = 0.85) +
  labs(title = "Distribution of characteristics across clusters (percentage)", x = "Cluster number", y = NULL) +
  geom_hline(yintercept = 2.5) +
  geom_hline(yintercept = 4.5) +
  geom_hline(yintercept = 6.5) +
  geom_hline(yintercept = 9.5) +
  geom_hline(yintercept = 12.5) +
  geom_hline(yintercept = 15.5) +
  geom_hline(yintercept = 17.5) +
  geom_hline(yintercept = 19.5) +
  geom_hline(yintercept = 22.5) +
  geom_hline(yintercept = 26.5) +
  geom_hline(yintercept = 29.5) +
  geom_hline(yintercept = 31.5) +
  scale_fill_gradient2(
    low = "darkslategray1",
    mid = "yellow",
    high =  "turquoise4"
  )

heatmap.p
ggsave("figures/hclust_characteristics_qual.png",
       width = 10,
       height = 15)

## --------- Values of continuous variables per cluster ---------
library(viridis)

clust_df <- cbind(df, clust.num)

colnames(clust_df)[length(colnames(clust_df))] <- "cluster"

#change to character for plotting
clust_df$cluster <- as.character(clust_df$cluster)

# fertile stamens
c1 <- ggplot(clust_df,
             aes(x = cluster,
                 y = Numberoffertilestamens,
                 fill = cluster)) + geom_boxplot(outlier.shape = NA, width = 0.2) +
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
  scale_y_continuous(name = "Number of fertile stamens") +
  scale_x_discrete(name = "cluster") +
  scale_fill_viridis(discrete = TRUE, alpha = 0.3) +
  scale_color_viridis(discrete = TRUE, alpha = 0.8)

# ovules per functional carpel
c2 <- ggplot(clust_df,
             aes(x = cluster,
                 y = Numberofovulesperfunctionalcarpel,
                 fill = cluster)) + geom_boxplot(outlier.shape = NA, width = 0.2) +
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
  scale_y_continuous(name = "Number of ovules per functional carpel") +
  scale_x_discrete(name = "cluster") +
  scale_fill_viridis(discrete = TRUE, alpha = 0.3) +
  scale_color_viridis(discrete = TRUE, alpha = 0.8)

# Flower size
c3 <- ggplot(clust_df, aes(x = cluster,
                           y = flowerSize,
                           fill = cluster)) + geom_boxplot(outlier.shape = NA, width = 0.2) +
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
  scale_y_continuous(name = "Flower diameter") +
  scale_x_discrete(name = "cluster") +
  scale_fill_viridis(discrete = TRUE, alpha = 0.3) +
  scale_color_viridis(discrete = TRUE, alpha = 0.8)

# Structural carpels
c4 <- ggplot(clust_df,
             aes(x = cluster,
                 y = Numberofstructuralcarpels,
                 fill = cluster)) + geom_boxplot(outlier.shape = NA, width = 0.2) +
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
  scale_y_continuous(name = "Number of structural carpels") +
  scale_x_discrete(name = "cluster") +
  scale_fill_viridis(discrete = TRUE, alpha = 0.3) +
  scale_color_viridis(discrete = TRUE, alpha = 0.8)

# Fusion of ovaries
c5 <- ggplot(clust_df, aes(x = cluster,
                           y = Fusionofovaries,
                           fill = cluster)) + geom_boxplot(outlier.shape = NA, width = 0.2) +
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
  scale_y_continuous(name = "Fusion of ovaries") +
  scale_x_discrete(name = "cluster") +
  scale_fill_viridis(discrete = TRUE, alpha = 0.3) +
  scale_color_viridis(discrete = TRUE, alpha = 0.8)

# Maximum vertical height
c6 <- ggplot(clust_df,
             aes(x = cluster,
                 y = Maximumverticalheight,
                 fill = cluster)) + geom_boxplot(outlier.shape = NA, width = 0.2) +
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
  scale_y_continuous(name = "Maximum vertical height") +
  scale_x_discrete(name = "cluster") +
  scale_fill_viridis(discrete = TRUE, alpha = 0.3) +
  scale_color_viridis(discrete = TRUE, alpha = 0.8)


(c1 + c2) / (c3 + c4) / (c5 + c6)
ggsave("figures/hclust_characteristics_quant.png",
       width = 15,
       height = 20)


## --------- PCOA scatterplot with cluster annotation ---------

#convert gower df to distance matrix
dataset_dist <- stats::as.dist(gower_df)

#run PCOA
dataset_pcoa <- ape::pcoa(dataset_dist)

#select number of clusters
# divisive.clust = 3
# aggl.clust.c (complete) = 4
# aggl.clust.a (average) = 2
# aggl.clust.w (ward) = 5
clust.num <- cutree(aggl.clust.w, k = 5)

#make matrix for downstream use
clust.num.stack <- table(stack(clust.num))
clust.num.stack <- as.data.frame.matrix(clust.num.stack)
rownames(clust.num.stack) <- c("cluster1",
                               "cluster2",
                               "cluster3",
                               "cluster4",
                               "cluster5")
saveRDS(clust.num.stack, "outputs/clust.num.stack5.Rds")

#get subset of species names to highlight
sp_names <- rownames(dataset_pcoa$vectors)
#prune down species name to make readable
inds <- round (runif(320, 1, length(sp_names)))
sp_names[inds] <- NA

library(ggrepel)

#plot points on first two axes, coloured by cluster
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
  ) +  geom_text_repel(aes(label = sp_names, colour = as.factor(clust.num)),
                       size = 3.5) +  stat_ellipse(geom = "polygon",
                                                   aes(fill = as.factor(clust.num)),
                                                   alpha = 0.25) +
  xlab(paste(
    "Axis 1: relative eigenvalue =",
    round(dataset_pcoa$values$Relative_eig[1], 2)
  )) +
  ylab(paste(
    "Axis 2: relative eigenvalue =",
    round(dataset_pcoa$values$Relative_eig[2], 2)
  ))

ggsave("figures/pcoa_hclust_k5.png",
       width = 12,
       height = 10)

#same but for k = 3
clust.num <- cutree(aggl.clust.w, k = 3)
#make matrix for downstream use
clust.num.stack <- table(stack(clust.num))
clust.num.stack <- as.data.frame.matrix(clust.num.stack)
rownames(clust.num.stack) <- c("cluster1",
                               "cluster2",
                               "cluster3")
saveRDS(clust.num.stack, "outputs/clust.num.stack3.Rds")


#plot points on first two axes, coloured by cluster
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
  ) +
  stat_ellipse(geom = "polygon",
               aes(fill = as.factor(clust.num)),
               alpha = 0.25) +
  xlab(paste(
    "Axis 1: relative eigenvalue =",
    round(dataset_pcoa$values$Relative_eig[1], 2)
  )) +
  ylab(paste(
    "Axis 2: relative eigenvalue =",
    round(dataset_pcoa$values$Relative_eig[2], 2)
  ))

ggsave("figures/pcoa_hclust_k3.png",
       width = 12,
       height = 10)

####
# Sankey plot
####

#from: https://r-graph-gallery.com/321-introduction-to-interactive-sankey-diagram-2.html
#table of different k values (2-7)

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

saveRDS(clust.num.k.2.7.df, file = here::here("outputs/clust_num_k_2_7_ward.rds"))


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

saveNetwork(p, "figures/sankey_ward.html")
