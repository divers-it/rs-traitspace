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
library(data.table)

#load formatted data
df <- readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

#dissimilarity matrix
gower_df <- daisy(df,
                  metric = "gower")
summary(gower_df)

## ------------ DIVISIVE CLUSTERING ------------
divisive.clust <- diana(as.matrix(gower_df),
                        diss = TRUE, keep.diss = TRUE)

pltree(divisive.clust, main = "Divisive", cex = 0.25)

## ------------ AGGLOMERATIVE CLUSTERING ------------

#running different methods to compare results
aggl.clust.c <- hclust(gower_df, method = "complete")
plot(aggl.clust.c, main = "Agglomerative, complete linkages", cex = 0.25)

aggl.clust.a <- hclust(gower_df, method = "average")
plot(aggl.clust.a, main = "Agglomerative, average", cex = 0.25)

aggl.clust.w <- hclust(gower_df, method = "ward.D2")
plot(aggl.clust.w, main = "Agglomerative, Ward", cex = 0.25)

## ------------ CLUSTERING STATS ------------

#source function
source("R/cstats.table.R")

# produce stats table, capping the maximum number of clusters at 10

#divisive
stats.df.d <- cstats.table(gower_df, divisive.clust, 10)
stats.df.d

#complete
stats.df.aggl.c <- cstats.table(gower_df, aggl.clust.c, 10)
stats.df.aggl.c

#average
stats.df.aggl.a <- cstats.table(gower_df, aggl.clust.a, 10)
stats.df.aggl.a

#ward
stats.df.aggl.w <- cstats.table(gower_df, aggl.clust.w, 10)
stats.df.aggl.w

#stats of one clustering approach / K value
cluster.stats(d = gower_df, clustering = cutree(aggl.clust.w, 4))

## --------- Choosing the number of clusters - elbow withinness ---------

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
# k = 3

## --------- Choosing the number of clusters - silhouette ---------

#When it comes to silhouette assessment, the rule is you should choose the number that maximizes the silhouette coefficient because you want clusters that are distinctive (far) enough to be considered separate. The silhouette coefficient ranges between -1 and 1, with 1 indicating good consistency within clusters, -1 — not so good.
#Seems to be biased towards two clusters, I prefer elbow

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

## --------- Plotting clusters on dendrogram  ---------

#set palette
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

# dendrogram with clusters agglomerative ward k = 6
dendro <- as.dendrogram(aggl.clust.w)
dendro.col <- dendro %>%
  dendextend::set("branches_k_color",
      k = 6,
      value = brewer.pal(7, "Dark2")) %>%
  dendextend::set("branches_lwd", 0.6) %>%
  dendextend::set("labels_colors",
      value = c("darkslategray")) %>%
  dendextend::set("labels_cex", 0.5)
ggd1 <- as.ggdend(dendro.col)
wk6 <- ggplot(ggd1, theme = theme_minimal()) +
  labs(x = "Num. observations", y = "Height", title = "Dendrogram agglomerative ward, k = 6")

#plot coloured dendrogram
wk3 + wk6
ggsave("figures/11_dendrograms_ward2_coloured_by_cluster.png",
       width = 10,
       height = 10)

## --------- PCOA scatterplot with cluster annotation ---------

#convert gower df to distance matrix
dataset_dist <- stats::as.dist(gower_df)

#run PCOA
dataset_pcoa <- ape::pcoa(dataset_dist)

#select number of clusters to plot (6)
clust.num <- cutree(aggl.clust.w, k = 6)

#get subset of species names to highlight
sp_names <- rownames(dataset_pcoa$vectors)
inds <- round (runif(320, 1, length(sp_names)))
sp_names[inds] <- NA

#plot points on first two axes, coloured by cluster with species names
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

#save image
ggsave("figures/11_scatterplot_pcoa_wardD2_k6_coloured_by_cluster.png",
       width = 12,
       height = 10)

#same as above but for k = 3
clust.num <- cutree(aggl.clust.w, k = 3)

#plot points on first two axes, coloured by cluster with species names
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

ggsave("figures/11_scatterplot_pcoa_wardD2_k3_coloured_by_cluster.png",
       width = 12,
       height = 10)

###
# ---- Sankey plot ----
###
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

#save output for downstream use
saveRDS(clust.num.k.2.7.df, file = here::here("outputs/11_clust_num_k_2_7_ward.rds"))


# Build a connection data frame - a list of flows (links between clusters)
# with intensity for each flow (number of individuals that move)
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

# With networkD3, connection must be provided using id, not using real name like in the links dataframe
# so this must be added
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1
links

#remove rows where values are 0
links<-links[links$value>0,]

# Make and plot the network
# Not that the lines do not correspond to individuals but groups
# so a single line cannot be followed across the plot
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE)
p

#save as html
saveNetwork(p, "figures/11_sankey_wardD2.html")

###
# ---- Boxplots and stacked barplots for robust groups ----
###

#make label
#robust_group<-paste("kpro_robust_",robust_vect_kpro_full,sep="")

#check names
rownames(df)==rownames(clust.num.k.2.7.df)

#add label to group
df_labelled<-cbind(df,clust.num.k.2.7.df$`3clusters`)

#change colname for label
colnames(df_labelled)[length(colnames(df_labelled))]<-"cluster"

#empty list for plots
plot_list <- list()

#make plots
for(i in 1:(length(colnames(df_labelled))-1)){
  
  #for quantitative
  if(is.numeric(df_labelled[,i])){
    plot_list[[i]]<-ggplot(df_labelled, aes(x=cluster, y=!!as.name(colnames(df_labelled)[i]), fill=cluster)) + 
      geom_boxplot() + 
      geom_jitter(shape=16, position=position_jitter(0.1)) + 
      theme(legend.position = "none",axis.text.x = element_text(angle = 90),axis.title.x = element_blank(),plot.margin = unit(c(1,1,1,1), "cm"))  
  } else {   #for qualitative
    plot_list[[i]]<-ggplot(df_labelled, aes(x=cluster, fill = !!as.name(colnames(df_labelled)[i]))) +
      geom_bar(stat="count") + theme(axis.text.x = element_text(angle = 90),axis.title.x = element_blank(),plot.margin = unit(c(1,1,1,1), "cm"))
  }
  
}

#plot as multiple pages in PDF
pdf("figures/11_boxplots_stacked_barplots_wardD2_clusters_by_trait.pdf",width = 15,height = 15)

print(grid.arrange(grobs=plot_list[1:4],ncol=2,nrow=2))
print(grid.arrange(grobs=plot_list[5:8],ncol=2,nrow=2))
print(grid.arrange(grobs=plot_list[9:12],ncol=2,nrow=2))
print(grid.arrange(grobs=plot_list[13:16],ncol=2,nrow=2))
print(grid.arrange(grobs=plot_list[17:19],ncol=2,nrow=2))

dev.off()

###
# ---- Plot qualitative stats of robust groups ----
###

#add group size to robust group label
for (i in 1:length(unique(df_labelled$cluster))) {
  df_labelled$cluster[df_labelled$cluster %in% sort(unique(df_labelled$cluster))[i]] <-
    paste(
      sort(unique(df_labelled$cluster))[i],
      " (n = ",
      table(df_labelled$cluster)[i],
      ")",
      sep = ""
    )
  
}

#make as factor for grouping
df_labelled$cluster<-as.factor(df_labelled$cluster)

#qualitative only
facts <- unlist(lapply(df_labelled, is.factor))
df_temp<-df_labelled[ , facts]

#change table to long form and count combinations
df_temp_melt<-data.table::melt(df_temp,id.vars="cluster")
df_temp_melt_counts <- df_temp_melt %>% group_by(cluster,variable,value) %>% summarise(count=n())

#add new column to remove text labels if counts are <5
df_temp_melt_counts$label<-df_temp_melt_counts$value
df_temp_melt_counts$label[df_temp_melt_counts$count<3]<-NA

# NOT RUN: make new column for text size
#df_temp_melt_counts$text_size<-df_temp_melt_counts$count^(1/2)

#plot stacked barplots per robust group for each qualitative trait, with labels
ggplot(df_temp_melt_counts, aes(variable, count, fill = value)) +
  geom_col(position = 'stack') + facet_wrap(. ~ cluster, scales = "free")  + theme(
    legend.position = "none",
    axis.text.x = element_text(vjust = 0.5, hjust=1),
    axis.title.x = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) + geom_text(aes(size = count,label = label),
                position = position_stack(vjust = .5)) + coord_flip()

#save plot 
ggsave("figures/11_stacked_barplots_wardD2_traits_by_cluster.png",width=15,height=10)

