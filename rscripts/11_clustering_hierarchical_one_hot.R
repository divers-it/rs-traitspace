rm(list = ls())
library(dplyr)

#load formatted data
df2<-readRDS(file = here::here("outputs/df_filt_trans_one_hot.rds"))
#df2[sapply(df2, is.factor)] <- lapply(df2[sapply(df2, is.factor)],
#                                     as.integer)

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

#stats of one clustering approach / K value
cluster.stats(d = gower_df2, clustering = cutree(aggl.clust.w, 4))

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

# dendrogram with clusters agglomerative ward k = 4
dendro <- as.dendrogram(aggl.clust.w)
dendro.col <- dendro %>%
  set("branches_k_color",
      k = 4,
      value = brewer.pal(7, "Dark2")) %>%
  set("branches_lwd", 0.6) %>%
  set("labels_colors",
      value = c("darkslategray")) %>%
  set("labels_cex", 0.5)
ggd1 <- as.ggdend(dendro.col)
wk4 <- ggplot(ggd1, theme = theme_minimal()) +
  labs(x = "Num. observations", y = "Height", title = "Dendrogram agglomerative ward, k = 4")

# dendrogram with clusters agglomerative ward k = 5
dendro <- as.dendrogram(aggl.clust.w)
dendro.col <- dendro %>%
  set("branches_k_color",
      k = 5,
      value = brewer.pal(7, "Dark2")) %>%
  set("branches_lwd", 0.6) %>%
  set("labels_colors",
      value = c("darkslategray")) %>%
  set("labels_cex", 0.5)
ggd1 <- as.ggdend(dendro.col)
wk5 <- ggplot(ggd1, theme = theme_minimal()) +
  labs(x = "Num. observations", y = "Height", title = "Dendrogram agglomerative ward, k = 5")

library(patchwork)
wk4 + wk5
ggsave("figures/dendro_ward_k4_k5_one_hot.png",
       width = 20,
       height = 10)


## --------- PCOA scatterplot with cluster annotation ---------

#convert gower df to distance matrix
dataset_dist <- stats::as.dist(gower_df2)

#run PCOA
dataset_pcoa <- ape::pcoa(dataset_dist)

#select number of clusters
# divisive.clust = 3
# aggl.clust.c (complete) = 4
# aggl.clust.a (average) = 2
# aggl.clust.w (ward) = 5
clust.num <- cutree(aggl.clust.w, k = 4)

#make matrix for downstream use
clust.num.stack <- table(stack(clust.num))
clust.num.stack <- as.data.frame.matrix(clust.num.stack)
rownames(clust.num.stack) <- c("cluster1",
                               "cluster2",
                               "cluster3",
                               "cluster4")
saveRDS(clust.num.stack, "outputs/clust.num.one_hot.stack4.Rds")

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

ggsave("figures/pcoa_hclust_k4_one_hot.png",
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

saveRDS(clust.num.k.2.7.df, file = here::here("outputs/clust_num_k_2_7_ward_one_hot.rds"))

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

library(networkD3)
library(dplyr)
# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE)

p

saveNetwork(p, "figures/sankey_ward_one_hot.html")


#####
#Boxplots and stacked barplots for robust groups
#####

# library
library(ggplot2)
library(gridExtra)

#make label
#robust_group<-paste("kpro_robust_",robust_vect_kpro_full,sep="")

#check names
rownames(df2)==rownames(clust.num.k.2.7.df)

#add label to group
df_labelled<-cbind(df2,clust.num.k.2.7.df$`3clusters`)

#change colname for label
colnames(df_labelled)[length(colnames(df_labelled))]<-"cluster"

#empty list for plots
plot_list <- list()

#make plots
for(i in 1:(length(colnames(df_labelled))-1)){
  
  #for quantitative
  if(is.numeric(df_labelled[,i])){
    plot_list[[i]]<-ggplot(df_labelled, aes(x=cluster, y=!!as.name(colnames(df_labelled)[i]), fill=cluster)) + 
      geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.1)) + theme(legend.position = "none",axis.text.x = element_text(angle = 90),axis.title.x = element_blank(),plot.margin = unit(c(1,1,1,1), "cm"))  
  } else {   #for qualitative
    plot_list[[i]]<-ggplot(df_labelled, aes(x=cluster, fill = !!as.name(colnames(df_labelled)[i]))) +
      geom_bar(stat="count") + theme(axis.text.x = element_text(angle = 90),axis.title.x = element_blank(),plot.margin = unit(c(1,1,1,1), "cm"))
  }
  
}

pdf("figures/k3_ward_plots.pdf",width = 15,height = 15)

print(grid.arrange(grobs=plot_list[1:4],ncol=2,nrow=2))
print(grid.arrange(grobs=plot_list[5:8],ncol=2,nrow=2))
print(grid.arrange(grobs=plot_list[9:12],ncol=2,nrow=2))
print(grid.arrange(grobs=plot_list[13:16],ncol=2,nrow=2))
print(grid.arrange(grobs=plot_list[17:19],ncol=2,nrow=2))

dev.off()


###
# Plot qualitative stats of robust groups
###
library(data.table)

df_labelled$cluster<-as.character(df_labelled$cluster)

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

#make new column for text size
#df_temp_melt_counts$text_size<-df_temp_melt_counts$count^(1/2)

#plot stacked barplots per robust group for each qualitative trait, with labels
ggplot(df_temp_melt_counts, aes(variable, count, fill = value)) +
  geom_col(position = 'stack') + facet_wrap(. ~ cluster, scales = "free")  + theme(
    legend.position = "none",
    axis.text.x = element_text(vjust = 0.5, hjust=1),
    axis.title.x = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) + geom_text(aes(label = label),
                size = 2,
                position = position_stack(vjust = .5)) + coord_flip()

ggsave("figures/stacked_barplots_k3_ward_one_hot.png",width=15,height=15)
