rm(list=ls())

# load libraries
library(dplyr)
library(ggmosaic)
library(tibble)
library(factoextra)
library(dplyr)
library(cluster)
library(GGally)
library(corrr)
library(tidyverse)
library(rcompanion)
library(harrypotter)
library(vegan)
library(patchwork)

# load data
df <- readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

## Trait combinations ----

# One-liner to look at frequency of trait combinations
combo_df <-
  df %>% group_by(SexualSystem, FlowerSex, .drop = FALSE) %>%
  summarize(count = n())
combo_df

# mosaic plot to reperesent combinations visually
ggplot(data = df) +
  geom_mosaic(aes(x = product(Pollination,Mating), fill=Mating), na.rm=TRUE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#### 
## Hierarchical clustering for variables ----
#### 

# https://rpubs.com/pjmurphy/269609

# split data set into quant and qualitative columns
df_num <- df[sapply(df,is.numeric)]
df_cat <- df[sapply(df,is.factor)]

# do clustering and plot tree
tree <- ClustOfVar::hclustvar(df_num, df_cat)
plot(tree)

# get clustering for k = 4
order_df <- data.frame(tree$label,tree$clusmat[,4])
colnames(order_df)<-c("trait","cluster")
order_df <- order_df[order(order_df$cluster),]
order_df

# estimate number of clusters using stability of partitions
stab <- ClustOfVar::stability(tree,
                              B=100, # number of bootstraps
                              graph = T
)

# plot stability results
plot(stab)
boxplot(stab$matCR[,1:7])

#### 
## Heatmap ----
#### 

# load original data set
df <- readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

# Calculate a pairwise association between all variables in a data-frame. 
# In particular nominal vs nominal with Chi-square,
# numeric vs numeric with Pearson correlation,
# and nominal vs numeric with ANOVA.
source("R/mixed_assoc.R")

# correlation matrix
cor_mat_ori <- df %>%
  mixed_assoc() %>%
  select(x, y, assoc) %>%
  spread(y, assoc) %>%
  column_to_rownames("x") %>%
  as.matrix %>%
  as_cordf

# melt correlated matrix into long table
melted_cor_mat_ori <- reshape2::melt(cor_mat_ori)
head(melted_cor_mat_ori)

# get absolute values
# melted_cor_mat_ori$value <- abs(melted_cor_mat_ori$value)

# cluster variables and get order of clustering
# NOTE: could replace with ClustOfVar dendrogram
ord <- hclust(dist(cor_mat_ori, method = "euclidean"), method = "centroid" )$order

# reorder melted table based on clustering
melted_cor_mat_ori$term <- factor(melted_cor_mat_ori$term, levels=tree$labels[tree$order])
melted_cor_mat_ori$variable <- factor(melted_cor_mat_ori$variable, levels=tree$labels[tree$order])

# replace NA with 1
melted_cor_mat_ori$value[is.na(melted_cor_mat_ori$value)] <- 1

####
## Figure SX: joined correlation heatmap and dendrogram ----
####

# get colours
cols<-c(rev(harrypotter::hp(2,option="Ravenclaw")),harrypotter::hp(2,option="LunaLovegood"))

# plot correlation matrix
c1 <- ggplot(data = melted_cor_mat_ori, aes(x=term, y=variable, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2(low = cols[1], mid = "white", high = cols[2], midpoint = 0, name = "Correlation",
                       limits = c(-0.5, 1), 
                       breaks = c(-0.5, 0, 0.5, 1),
                       labels = c(-0.5, 0, 0.5, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position="left",) +
  xlab("") +
  ylab('')

c1

# add dendrogram
d1 <- ggdendro::ggdendrogram(as.dendrogram(tree), rotate=TRUE)
d1 + theme_classic()
d2 <- d1 + ggdendro::theme_dendro() + theme(plot.margin = unit(c(0,0,0,0), "cm"))
# + scale_y_reverse(expand = c(0.2, 0))
d2

# combined plot
# NOTE: dendrogram not exactly aligned.
patch <- c1 + d2 + plot_layout(widths=c(4, 1))
patch + plot_annotation(tag_levels = 'a',tag_prefix="(",tag_suffix=")") & theme(plot.tag = element_text(size = 14))

# save pdf to be modified
ggsave("figures/figure_SX_correlation_heatmap_tree.pdf",width=14,height=10)

# NOTE: could try an alternative method with
# install.packages("plotly")
# install.packages("heatmaply")
# https://talgalili.github.io/heatmaply/

# save
# ggsave("figures/one_hot_7_heatmap_abs.png",
#        width = 20,
#        height = 15,
#        units = 'cm')

#### 
## Correlation matrix ----
#### 

# load one-hot data set
df <- readRDS(file = here::here("outputs/one_hot_6_df_filt_trans.rds"))

# load original data set
# df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

df %>%
  mixed_assoc()

# correlation matrix
cor_mat<- df %>%
  mixed_assoc() %>%
  select(x, y, assoc) %>%
  spread(y, assoc) %>%
  column_to_rownames("x") %>%
  as.matrix %>%
  as_cordf

# look at correlations
cor_mat

#### 
## Figure S4: Network plot adapted from code of corrr::network_plot ----
#### 

# minimum correlation allowed
min_cor <- 0.3
curved <- FALSE

# reformat matrix
rdf <- as_matrix(cor_mat, diagonal = 1)
distance <- 1 - abs(rdf)

### Choose MDS approach ----

## PCOA
# points<-stats::cmdscale(distance, k = 2)
# points <- data.frame(points)
# colnames(points) <- c("x", "y")
# points$id <- rownames(points)

## NMDS
nmds <-
  metaMDS(distance,
          # distance = "gower",
          k = 2,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)
points <- data.frame(nmds$points)
colnames(points) <- c("x", "y")
points$id <- rownames(points)

# Create a proximity matrix of the paths to be plotted.
proximity <- abs(rdf)
proximity[upper.tri(proximity)] <- NA
diag(proximity) <- NA
proximity[proximity < min_cor] <- NA

# Produce a data frame of data needed for plotting the paths.
n_paths <- sum(!is.na(proximity))
paths <- data.frame(matrix(nrow = n_paths, ncol = 6))
colnames(paths) <- c("x", "y", "xend", "yend", "proximity", "sign")
path <- 1
for (row in 1:nrow(proximity)) {
  for (col in 1:ncol(proximity)) {
    path_proximity <- proximity[row, col]
    if (!is.na(path_proximity)) {
      path_sign <- sign(rdf[row, col])
      x <- points$x[row]
      y <- points$y[row]
      xend <- points$x[col]
      yend <- points$y[col]
      paths[path, ] <- c(x, y, xend, yend, path_proximity, path_sign)
      path <- path + 1
    }
  }
}

# look at output
paths

# add point size columns
points$freq<-rep(NA,length(points$id))

for(i in 1:length(rownames(points))){
  
  points$freq[i] <- sum(na.omit(df[,points$id[i]]) != 0)
  
}


# trait type column
points$type<-rep("reproductive",length(points$id))
points$type[grep("Aqua",points$id)] <- "vegetative"
points$type[grep("Climb",points$id)] <- "vegetative"
points$type[grep("Dispersal",points$id)] <- "vegetative"
points$type[grep("Lifespan",points$id)] <- "vegetative"
points$type[grep("height",points$id)] <- "vegetative"
points$type[grep("Wood",points$id)] <- "vegetative"
head(points)

# get colours
cols<-c(rev(harrypotter::hp(2,option="Ravenclaw")),harrypotter::hp(2,option="LunaLovegood"))

# make plot
ggplot() +
  geom_point(
    data = points,
    aes(x, y, size = freq),
    shape = 21,
    fill = "#EDEDED",
    colour = "#BDBDBD",
    alpha = 0.6,
  ) + 
  geom_curve(
    data = paths,
    alpha = 0.5,
    size = 1,
    curvature = 0.2,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      colour = proximity * sign
    )
  ) +
  # geom_segment(
  #   data = paths,
  #   alpha = 0.5,
  #   size = 1,
  #   aes(
  #     x = x,
  #     y = y,
  #     xend = xend,
  #     yend = yend,
  #     colour = proximity * sign
  #   )
# ) + 
geom_point(
  data = points,
  aes(x, y, fill = type),
  size=3.5,
  shape = 21,
  alpha=0.7,
  colour = "#CDCDCD"
) +
  geom_point(
    data = points,
    aes(x, y),
    size=1.5,
    shape = 16,
    alpha=0.7,
    colour = "black"
  ) +
  scale_fill_manual(values = c(cols[4],cols[3]), name="Trait type") +
  scale_colour_gradientn(colours = c(cols[1],"white",cols[2]), name = "Correlation", values = c(0,
                                                                                                abs(min(paths$proximity * paths$sign))/(max(paths$proximity * paths$sign)+abs(min(paths$proximity * paths$sign))), # to calculate where 0 point is in scale
                                                                                                1)) +
  scale_size(range=c(1,15), name = "Frequency",breaks = c(25,50,100,200,300)) + 
  scale_alpha(range=c(0,1)) + 
  ggrepel::geom_text_repel(
    data = points,
    aes(x, y, label = id),
    size = 5,
    segment.size = 0.0,
    segment.color = "white"
  ) + 
  theme_void() +
  guides(alpha = "none",
         size = guide_legend(override.aes = list(shape = 21, alpha=0.5, fill = "#EDEDED", colour = "#CDCDCD")),
         # shape = guide_legend(override.aes = list(shape = 21, fill = "#EDEDED", colour = "#CDCDCD")),
         # fill = guide_legend(override.aes = list(shape = 21, fill = "#EDEDED", colour = "#CDCDCD"))
  ) +
  theme(legend.key = element_blank(),
        legend.background = element_blank(),
        legend.box = "horizontal",
        legend.box.background = element_rect(color = "black"),
        legend.box.margin = margin(5,5,5,5),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12),
        legend.position = c(0.85,0.1),
        panel.background = element_rect(fill='white'))


ggsave("figures/figure_S4_correlation_network.png",width=14,height=14)

#### 
## Original corrr function for comparison ----
#### 

network_plot(cor_mat,colours = c("indianred2", "white", "skyblue1"),min_cor=0.3,repel=TRUE)

#### 
## ---- igraph network plot -----
#### 

library(GGally)
library(network)
library(sna)
library(ggplot2)
library(tidyverse)
library(corrr)
library(igraph)
library(ggraph)

net = rgraph(10, mode = "graph", tprob = 0.5)

# as matrix
cor_mat <- as.matrix(cor_mat,rownames.force=TRUE)
rownames(cor_mat)<-cor_mat[,1]
cor_mat <- cor_mat[,-1]
as.numeric(cor_mat)

# Keep only high correlations
cor_mat[cor_mat<0.25] <- 0

# Make an Igraph object from this matrix:
network <- graph_from_adjacency_matrix( cor_mat, weighted=T, mode="undirected", diag=F)

ggraph(network) +
 geom_edge_link(aes(edge_alpha = abs(weight), edge_width = 1, color = weight)) +
 guides(edge_alpha = "none", edge_width = "none") +
 scale_edge_colour_gradientn(limits = c(min(cor_mat), max(cor_mat)), colors = c("white", "dodgerblue2")) +
 geom_node_point(color = "grey", size = 2) +
 geom_node_text(aes(label = name), repel = TRUE) +
 theme_graph() +
 labs(title = "Correlations between DiveRS traits")

