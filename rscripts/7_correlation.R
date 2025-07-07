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

## Make names nice for plots ----
new_names <- c("Maximum height",
 "No. fertile stamens",
 "No. ovules per carpel",
 "No. structural carpels",
 "Fusion of ovaries",
 "Flower size",
 "Seed mass",
 "Woodiness",
 "Climbing",
 "Aquatic",
 "Sexual system",
 "Lifespan",
 "Mating system",
 "Pollination",
 "Dispersal mode",
 "Dispersal distance",
 "Flower sex",
 "Ovary position",
 "Floral reward",
 "Flower symmetry",
 "Showiness")

# check correspondence
data.frame(colnames(df), new_names)

# assign new names
colnames(df) <- new_names


## NOT RUN: Trait combinations ----

# # One-liner to look at frequency of trait combinations
# combo_df <-
#   df %>% group_by(SexualSystem, FlowerSex, .drop = FALSE) %>%
#   summarize(count = n())
# combo_df
# 
# # mosaic plot to reperesent combinations visually
# ggplot(data = df) +
#   geom_mosaic(aes(x = product(Pollination,Mating), fill=Mating), na.rm=TRUE) + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

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

# make names nice 

# melt correlated matrix into long table
melted_cor_mat_ori <- reshape2::melt(cor_mat_ori)
head(melted_cor_mat_ori)

# get absolute values
# melted_cor_mat_ori$value <- abs(melted_cor_mat_ori$value)

# reorder melted table based on clustering
melted_cor_mat_ori$term <- factor(melted_cor_mat_ori$term, levels=tree$labels[tree$order])
melted_cor_mat_ori$variable <- factor(melted_cor_mat_ori$variable, levels=tree$labels[tree$order])

# replace NA with 1
melted_cor_mat_ori$value[is.na(melted_cor_mat_ori$value)] <- 1

####
## Figure 2b: Joined correlation heatmap and dendrogram ----
####

# get colours
cols<-c(rev(harrypotter::hp(2,option="Ravenclaw")),harrypotter::hp(2,option="LunaLovegood"))

# plot correlation matrix
c1 <- ggplot(data = melted_cor_mat_ori, aes(x=term, y=variable, fill=abs(value))) + 
  geom_tile() +
  scale_fill_gradient2(low = "white", high = cols[2], midpoint = 0, name = "Correlation",
                       limits = c(0, 1), 
                       breaks = c(0, 0.5, 1),
                       labels = c(0, 0.5, 1)) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        legend.position="none",) +
  xlab("") +
  ylab('')

c1

# add dendrogram
d1 <- ggdendro::ggdendrogram(as.dendrogram(tree), rotate=F)
d1
d2 <- d1 + ggdendro::theme_dendro() + theme(plot.margin = unit(c(0,0,0,0), "cm"))
# + scale_y_reverse(expand = c(0.2, 0))
d2

# combined plot
# NOTE: dendrogram not exactly aligned.
patch <- ( d2 / c1 ) + plot_layout(heights=c(1, 4))
patch + plot_annotation(tag_levels = 'a',tag_prefix="(",tag_suffix=")") & theme(plot.tag = element_text(size = 14))

# save plots for combined figure with loadings plot
saveRDS(c1, file = "outputs/7_correlation_heatmap.rds")
saveRDS(d2, file = "outputs/7_correlation_dendrogram.rds")

# save pdf to be modified
# ggsave("figures/figure_SX_correlation_heatmap_tree.pdf",width=14,height=10)

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

## Not used: Network plot adapted from code of corrr::network_plot ----

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


# ggsave("figures/figure_S4_correlation_network.png",width=14,height=14)

#### 
### Original corrr function for comparison ----
#### 

network_plot(cor_mat,colours = c("indianred2", "white", "skyblue1"),min_cor=0.3,repel=TRUE)
