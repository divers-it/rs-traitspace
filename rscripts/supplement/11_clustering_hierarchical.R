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
wk3# + wk6
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

#Recalculate relative eigenvalues by removing negative eigenvalues as in Mouillot et al.  
ev_pcoa <- dataset_pcoa$values$Eigenvalues
ev_pcoa_g0 <- ev_pcoa[ev_pcoa>0]
rel_ev_pcoa_g0 <- ev_pcoa_g0/sum(ev_pcoa_g0)

#plot points on first two axes, coloured by cluster with species names
ggplot(data.frame(dataset_pcoa$vectors),
       aes(
         x = Axis.3,
         y = Axis.4,
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
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2))) + 
  theme(legend.position = "none")

#save image
ggsave("figures/11_scatterplot_pcoa_wardD2_k6_coloured_by_cluster.png",
       width = 10,
       height = 10)

#same as above but for k = 3
clust.num <- cutree(aggl.clust.w, k = 3)

cols<-harrypotter::hp(3,option="ronweasley2")

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


ggsave("figures/11_scatterplot_pcoa_wardD2_k3_coloured_by_cluster.png",
       width = 10,
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

####
# ---- Figure 3b: Quantitative trait boxplots for robust clusters ----
####

#palette from scatterplot
cols<-harrypotter::hp(3,option="ronweasley2")

b1 <- ggplot(df_labelled, aes(x=cluster, y=Maximumverticalheight, fill=cluster)) + 
  geom_boxplot(alpha=0.7) + 
  geom_jitter(shape=21, position=position_jitter(0.1),alpha=0.7) + 
  scale_fill_manual(values = c(cols)) +
  scale_y_continuous(limits = quantile(df_labelled$Maximumverticalheight, c(0.05, 0.95),na.rm = TRUE)) +
  ylab("Maximum vertical height") +
  theme(legend.position = "none",
        # add border 1)
        # panel.border = element_blank(),
        # color background 2)
        panel.background = element_rect(fill = "white"),
        # modify grid 3)
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color="grey70"),
        panel.grid.minor.y = element_blank(),
        # modify text, axis and colour 4) and 5)
        axis.line.y = element_line(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        axis.ticks.x = element_blank(),
        #axis.ticks.y = element_blank()
  )#   + 
#annotate("text", x = -.5, y = 1.5, label = "(b)", size = 8) +
#coord_cartesian(ylim = c(-1.5, 1.75), xlim = c(1, 8), clip = "off") +
#theme(plot.margin = unit(c(1,1,1,3), "lines"))

b1

b2 <- ggplot(df_labelled, aes(x=cluster, y=flowerSize, fill=cluster)) + 
  geom_boxplot(alpha=0.7) + 
  geom_jitter(shape=21, position=position_jitter(0.1),alpha=0.7) + 
  scale_fill_manual(values = c(cols,"grey")) +
  scale_y_continuous(limits = quantile(df_labelled$flowerSize, c(0.05, 0.95),na.rm = TRUE)) +
  ylab("Flower size") +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1), "cm"),
        # add border 1)
        #panel.border = element_blank(),
        # color background 2)
        panel.background = element_rect(fill = "white"),
        # modify grid 3)
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color="grey70"),
        panel.grid.minor.y = element_blank(),
        # modify text, axis and colour 4) and 5)
        axis.line.y = element_line(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=20),
        axis.ticks.x = element_blank(),
        #axis.ticks.y = element_blank()
  )  

b2

b1 / b2

#save plot 
ggsave("figures/11_boxplots_wardD2_traits_by_cluster.png",width=10,height=10)

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


#theme
my_theme <- function() {
  theme(
    # add border 1)
    panel.border = element_blank(),
    # color background 2)
    panel.background = element_rect(fill = "white"),
    # modify grid 3)
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    # modify text, axis and colour 4) and 5)
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=16),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_line(),
    axis.ticks.length=unit(.25, "cm"),
    axis.ticks.y = element_blank(),
    # legend
    legend.position = "none",
    # margin
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
  )
}

## Robust group 1

#subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$cluster=="k_3_cluster_1 (n = 172)",]

#get palette based on max counts
pal1<-colorRampPalette(c("white",cols[1]))(max(rob_df$count))

#make new column for colours based on count
rob_df$pal <- pal1[rob_df$count]
rob_df$pal[is.na(rob_df$value)]<-"grey20"

p1 <- ggplot(rob_df, aes(x = variable, y = count, fill=pal, alpha=0.98)) +
  geom_bar(position="stack", stat="identity",col="black") +
  scale_fill_identity() + 
  geom_text(aes(size = count, label = label), position = position_stack(vjust = 0.5)) + 
  coord_flip() +
  my_theme()

p1


## Robust group 2

#subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$cluster=="k_3_cluster_2 (n = 87)",]

#get palette based on max counts
pal2<-colorRampPalette(c("white",cols[2]))(max(rob_df$count))

#make new column for colours based on count
rob_df$pal <- pal2[rob_df$count]
rob_df$pal[is.na(rob_df$value)]<-"grey20"

p2 <- ggplot(rob_df, aes(x = variable, y = count, fill=pal, alpha=0.98)) +
  geom_bar(position="stack", stat="identity",col="black") +
  scale_fill_identity() + 
  geom_text(aes(size = count, label = label), position = position_stack(vjust = 0.5)) + 
  coord_flip() +
  my_theme()

p2


## Robust group 3

#subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$cluster=="k_3_cluster_3 (n = 101)",]

#get palette based on max counts
pal3<-colorRampPalette(c("white",cols[3]))(max(rob_df$count))

#make new column for colours based on count
rob_df$pal <- pal3[rob_df$count]
rob_df$pal[is.na(rob_df$value)]<-"grey20"

p3 <- ggplot(rob_df, aes(x = variable, y = count, fill=pal, alpha=0.98)) +
  geom_bar(position="stack", stat="identity",col="black") +
  scale_fill_identity() + 
  geom_text(aes(size = count, label = label), position = position_stack(vjust = 0.5)) + 
  coord_flip() +
  my_theme()

p3

p1 + p2 + p3

#save plot 
ggsave("figures/11_stacked_barplots_wardD2_traits_by_cluster.png",width=20,height=10)

#combined plot
patch <- ( p1 + p2 + p3 ) / ( b1 + b2 ) + plot_layout(heights=c(2, 1))

patch + plot_annotation(tag_levels = 'a',tag_prefix="(",tag_suffix=")") & theme(plot.tag = element_text(size = 14))
#patch

ggsave("figures/11_scatterplot_boxplots_and_stacked_barplots.png",width=20,height=15)
