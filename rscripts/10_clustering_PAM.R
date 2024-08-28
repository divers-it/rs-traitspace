rm(list = ls())

#load packages
library(dplyr)
library(gridExtra)
library(Rtsne)
library(ggplot2)
library(cluster)
library(RColorBrewer)
library(ggrepel)
library(networkD3)
library(patchwork)

# load formatted data
df <- readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

# build distance matrix
gower_df <- daisy(df,
                  metric = "gower" )
summary(gower_df)

## Select number of clusters (k) ----

# Silhouette Width to select the optimal number of clusters
# The silhouette width is one of the very popular choices when it comes to selecting the optimal number of clusters. 
# It measures the similarity of each point to its cluster, and compares that to the similarity of the point with the closest neighboring cluster. 
# This metric ranges between -1 to 1, where a higher value implies better similarity of the points to their clusters. 
# Therefore, a higher value of the Silhouette Width is desirable. 
# We calculate this metric for a range of cluster numbers and find where it is maximized. 
# The following code shows the implementation in R:

#empty vector
silhouette <- c()

#run PAM with different values of K from 2-10 and calculate silhouette width
for(i in 2:10){
  pam_clusters <- pam(as.matrix(gower_df),
                      diss = TRUE,
                      k = i)
  silhouette <- c(silhouette ,pam_clusters$silinfo$avg.width)
}

### Plot silhouette width ----
plot(2:10, silhouette,
     xlab = "Clusters",
     ylab = "Silhouette Width")
lines(2:10, silhouette)

# make df with cluster membership for each value of k from 2-7
for(i in 2:7){
  if(i == 2){
    pam.gower = pam(gower_df, diss = TRUE, k = i)
    pam_df<-pam.gower$clustering
  } else {
    pam.gower = pam(gower_df, diss = TRUE, k = i)
    pam_df<-cbind(pam_df,pam.gower$clustering)
  }
}

colnames(pam_df)<-c("2clusters",
                    "3clusters",
                    "4clusters",
                    "5clusters",
                    "6clusters",
                    "7clusters")
clust.num.k.2.7.df <-as.data.frame(pam_df)
rownames(clust.num.k.2.7.df)<-names(pam.gower$clustering)

####
## Sankey plot ----
####

#from: https://r-graph-gallery.com/321-introduction-to-interactive-sankey-diagram-2.html

#table of different k values (2-7)
for (i in 1:6) {
  if (i == 1) {
    clust.num.k.2.7 <- paste("k",i+1,"cluster",as.character(clust.num.k.2.7.df[,i]),sep="_")
  } else {
    clust.num.k.2.7 <-
      cbind(clust.num.k.2.7, paste("k",i+1,"cluster",as.character(clust.num.k.2.7.df[,i]),sep="_"))
  }
}

colnames(clust.num.k.2.7)<-c("2clusters",
                             "3clusters",
                             "4clusters",
                             "5clusters",
                             "6clusters",
                             "7clusters")
clust.num.k.2.7.df <-as.data.frame(clust.num.k.2.7)

# fix rownames again
rownames(clust.num.k.2.7.df)<-names(pam.gower$clustering)

# save RDS
saveRDS(clust.num.k.2.7.df, file = here::here("outputs/10_clust_num_k_2_7_pam.rds"))

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

# remove rows where values are 0
links<-links[links$value>0,]

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", fontSize = 14,
                   sinksRight=FALSE)

p

# save as HTML
saveNetwork(p, "figures/10_sankey_pam.html")

####
## Identify robust groups ----
####
# Robust groups are those that consistently stay together as values of k change

# make data frame of combo frequencies
combos <- as.data.frame(table(clust.num.k.2.7.df))

# remove no existant combos
combos <- combos[combos$Freq > 0, ]

# order
combos <- combos[order(combos$Freq, decreasing = T), ]

# change to strings
combos <-data.frame(lapply(combos, as.character), stringsAsFactors = FALSE)
plot(combos$Freq)

# empty list
robust<-list()

# empty vector
robust_vect_pam<-rep(NA,length(rownames(df)))
names(robust_vect_pam)<-rownames(df)

# loop through ordered table to extract robust groups
for(i in 1:length(combos$Freq[as.numeric(combos$Freq)>0])){
  foo<-rownames(clust.num.k.2.7.df[clust.num.k.2.7.df[, 1] == combos[i, 1] & 
                                     clust.num.k.2.7.df[, 2] == combos[i, 2] &
                                     clust.num.k.2.7.df[, 3] == combos[i, 3] &
                                     clust.num.k.2.7.df[, 4] == combos[i, 4] &
                                     clust.num.k.2.7.df[, 5] == combos[i, 5] &
                                     clust.num.k.2.7.df[, 6] == combos[i, 6],])
  
  robust[[i]]<-foo
  
  robust_vect_pam[foo]<-i
  
}

# robust groups
robust

# keep 80% of the species in robust clusters ; others are NA
sum = 0

for (i in 1:max(robust_vect_pam)) {
  
  if (sum > 0.8 * length(robust_vect_pam)) {
    robust_vect_pam[robust_vect_pam == i] <- NA
    
  }
  
  sum = sum + length(robust[[i]])
  
}

table(robust_vect_pam)

# complete vector of robust groups and non-robust 
robust_vect_pam_full<-robust_vect_pam

saveRDS(robust_vect_pam_full, file = here::here("outputs/10_robust_vect_pam_full.rds"))

# remove species not in robust groups
robust_vect_pam<-na.omit(robust_vect_pam)

## PCOA ----

# run PCOA
dataset_dist <- stats::as.dist(gower_df)
dataset_pcoa <- ape::pcoa(dataset_dist)

# Recalculate relative eigenvalues by removing negative eigenvalues as in Mouillot et al.  
ev_pcoa <- dataset_pcoa$values$Eigenvalues
ev_pcoa_g0 <- ev_pcoa[ev_pcoa>0]
rel_ev_pcoa_g0 <- ev_pcoa_g0/sum(ev_pcoa_g0)

# reformat clustering df, assign NA to non-robust species
clusters_pcoa <- clust.num.k.2.7.df
colnames(clusters_pcoa) <- paste("n",colnames(clusters_pcoa),sep="")

# check alignment
rownames(clusters_pcoa) == names(robust_vect_pam_full)

# make non-robust NA
clusters_pcoa[is.na(robust_vect_pam_full),]<-NA

# empty list to store plots
pcoa_plots <- list()

# loop through different values of k
for(i in 1:length(clusters_pcoa[1,])){
  
  # Plot points on first two PCoA axes, coloured by cluster
  pcoa_plots[[i]] <- local({
    i<-i
    ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2, fill = as.factor(clusters_pcoa[,i]))) +
      geom_point(
        color="black",
        shape=21,
        alpha=0.5,
        size=3,
        stroke = 0.5
      ) + 
      ggtitle(paste(colnames(clusters_pcoa)[i])) +
      # stat_ellipse(geom = "polygon",
      #              aes(fill =  as.factor(clusters_pcoa[,i])), 
      #              alpha = 0.2) +
      xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
      ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2))) + 
      theme(legend.position = "none")})
  
  
}


( pcoa_plots[[1]] + pcoa_plots[[2]] ) / ( pcoa_plots[[3]] + pcoa_plots[[4]] ) / ( pcoa_plots[[5]] + pcoa_plots[[6]] )

# ggsave("figures/10_scatterplots_pcoa_pam_coloured_by_cluster.png",width=10,height=15)

####
## Figure S7: Scatterplot PAM with clusters ----
####

#convert gower df to distance matrix
dataset_dist <- stats::as.dist(gower_df)

#run PCOA
dataset_pcoa <- ape::pcoa(dataset_dist)

#Recalculate relative eigenvalues by removing negative eigenvalues as in Mouillot et al.  
ev_pcoa <- dataset_pcoa$values$Eigenvalues
ev_pcoa_g0 <- ev_pcoa[ev_pcoa>0]
rel_ev_pcoa_g0 <- ev_pcoa_g0/sum(ev_pcoa_g0)

# colors
cols<-harrypotter::hp(3,option="ronweasley2")

# plot points on first two axes, coloured by cluster with species names
ggplot(data.frame(dataset_pcoa$vectors),
       aes(
         x = Axis.1,
         y = Axis.2,
         fill = as.factor(clust.num.k.2.7.df$`3clusters`)
       )) +
  geom_point(
    color = "black",
    shape = 21,
    alpha = 0.5,
    size = 3,
    stroke = 0.5
  ) +  stat_ellipse(geom = "polygon", aes(fill = as.factor(clust.num.k.2.7.df$`3clusters`)), alpha = 0.25) +
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = cols) +
  xlab(paste("Axis 1: relative eigenvalue =", round(rel_ev_pcoa_g0[1], 2))) +
  ylab(paste("Axis 2: relative eigenvalue =", round(rel_ev_pcoa_g0[2], 2))) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16)
  )

ggsave("figures/figure_S7_pcoa_pam.png",
       width = 10,
       height = 10)


####
### Scatterplot PAM with robust groups showing unisexual, abiotic outliers ----
####

#empty matrix
rob_mat<-matrix(nrow = length(unique(robust_vect_pam)), ncol=length(df[1,]))

#empty matrix
rob_mat_names<-matrix(nrow = length(unique(robust_vect_pam)), ncol=length(df[1,]))

#loop through different robust groups
for(i in 1:length(unique(robust_vect_pam))){
  
  #names of species in robust group
  grp<-names(robust_vect_pam)[robust_vect_pam==i]
  
  #data from group
  grp_df<-df[rownames(df)%in%grp,]
  
  #loop through table
  for(j in 1:length(colnames(grp_df))){
    
    #for quantitative traits
    if(is.factor(grp_df[,j])){
      
      #frequency of most frequent state
      rob_mat[i,j]<-sort(table(grp_df[,j]),decreasing = T)[1] / length(grp_df[,j])
      
      #name of most frequent state
      names(sort(table(grp_df[,j]),decreasing = T)[1])
      rob_mat_names[i,j]<-names(sort(table(grp_df[,j]),decreasing = T)[1])
      
    } else {
      
      #mean of values
      rob_mat[i,j]<-mean(na.omit(grp_df[,j]))
      
    }
    
  }
  
}

#add row and column names
rownames(rob_mat)<-paste("robust",c(1:length(unique(robust_vect_pam))),sep="")
colnames(rob_mat)<-colnames(df)
rob_mat

rownames(rob_mat_names)<-paste("robust",c(1:length(unique(robust_vect_pam))),sep="")
colnames(rob_mat_names)<-colnames(df)
rob_mat_names

#check order
rownames(dataset_pcoa$vectors)==names(robust_vect_pam_full)

#get dataframe of robust groups
pcoa_robust<-cbind(dataset_pcoa$vectors,robust_vect_pam_full)

#check order
rownames(df)==rownames(pcoa_robust)

#plot points on first two PCoA axes, coloured by robust group and shaped by cluster
ggplot(
  data.frame(pcoa_robust),
  aes(
    x = Axis.1,
    y = Axis.2,
    fill = as.factor(robust_vect_pam_full),
    col = as.factor(robust_vect_pam_full)
  )
) +
  geom_point(
    aes(shape = as.factor(df$Pollination),
        size = as.factor(df$FlowerSex)),
    alpha = 0.5,
    stroke = 0.5
  ) +
  scale_shape_manual(values=c(23,22,8,24,4,21)) +
  labs(
    colour = "Robust group",
    size = "Flower sex",
    shape = "Pollination"
  ) + 
  theme_bw() +
  guides(col = guide_legend(override.aes = list(size = 5)),
         shape = guide_legend(override.aes = list(size = 5)),
         fill="none") +
  xlab("PCoA Axis 1") +
  ylab("PCoA Axis 2")

# ggsave("figures/10_pcoa_pam_coloured_by_robust.png",width=12,height=10)

####
# Figure 5a: UMAP PAM clustering / robust ----- 
####

devtools::install_github("ropenscilabs/ochRe")
library(ochRe)
library(umap)

# NOTE: Can comment out if already run
custom_config <- umap.defaults
custom_config$n_components <-  2 # number of dimensions targeted
custom_config$n_neighbors <- 10 # number of dimensions targeted
custom_config$input <- "dist" # The input matrix is a distance matrix

# run umap
umap_final <- umap(d = as.matrix(gower_df),config = custom_config)

# get output
df_umap_final <- data.frame(umap_final$layout)
rownames(df_umap_final) <- rownames(df) #CHECK
df_umap_final

# save output
saveRDS(df_umap_final,file="outputs/10_umap_nc2_nn10.rds")

# load previously run output
df_umap_final <- readRDS("outputs/10_umap_nc2_nn10.rds")

#PCoA scatterplot with density polygons
s1 <- ggplot(
  data.frame(df_umap_final),
  aes(
    x = X1,
    y = X2,
    fill = as.factor(robust_vect_pam_full)
  )
) +
  geom_point(
    aes(shape = as.factor(clust.num.k.2.7.df$`3clusters`)),
    alpha = 0.7,
    size = 7,
    stroke = 0.5) + 
  scale_fill_manual(values=c(ochre_pal("healthy_reef")(7)[1:4],"thistle",ochre_pal("healthy_reef")(7)[7])) +
  #scale_color_ochre(palette = "healthy_reef") +
  scale_shape_manual(values=c(21,22,23), labels=c('1', '2', '3')) +
  xlab("UMAP Axis 1") +
  ylab("UMAP Axis 2") +
  theme_bw() + theme(
    panel.border = element_blank(),
    #panel.grid.major = element_line(colour = "darkgrey"),
    #panel.grid.minor = element_line(colour = "grey"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.125, 0.75),
    legend.text = element_text(size=16),
    legend.title = element_text(size=18),
    axis.text.y = element_text(size=14),
    axis.title.y = element_text(size=20),
    axis.text.x = element_text(size=14),
    axis.title.x = element_text(size=20),
    axis.line = element_line(colour = "black")
  ) + 
  labs(
    colour = "Robust group",
    shape = "PAM cluster"
  ) + 
  guides(fill = guide_legend(override.aes = list(size = 8,
                                                 shape=21,
                                                 fill=c(ochre_pal("healthy_reef")(7)[1:4],
                                                        "thistle",
                                                        ochre_pal("healthy_reef")(7)[7],
                                                        "darkgrey")), title = "Robust group"),
         shape = guide_legend(override.aes = list(size = 8))
  ) + 
  #annotate("text", x = -4.5, y = 5.75, label = "(a)", size = 8) +
  coord_cartesian(xlim = c(-5, 5), ylim = c(-4, 6), clip = "off") +
  theme(plot.margin = unit(c(3,1,1,3), "lines"))


s1

# save plot
# ggsave("figures/10_scatterplot_umap_pam_coloured_by_robust.png",width=10,height=10)

#Proportion of missing data 
#species that dont belong to robust group
df_not_robust<-df[is.na(robust_vect_pam_full),]
mean(is.na(df_not_robust))

#species that do
df_robust<-df[!is.na(robust_vect_pam_full),]
mean(is.na(df_robust))

####
## Figure 5b: Quantitative trait boxplots for robust clusters ----
####

# make label
robust_group<-paste("pam_robust_",robust_vect_pam_full,sep="")

# add label to group
df_labelled<-cbind(df,robust_group)

# empty list for plots
plot_list <- list()

# palette from scatterplot
cols <- c(ochre_pal("healthy_reef")(7)[1:4],"thistle",ochre_pal("healthy_reef")(7)[7])

# make plots
b1 <- ggplot(df_labelled, aes(x=robust_group, y=Maximumverticalheight, fill=robust_group)) + 
  geom_boxplot(alpha=0.7) + 
  geom_jitter(shape=21, position=position_jitter(0.1),alpha=0.7) + 
  scale_fill_manual(values = c(cols,"grey")) +
  scale_y_continuous(limits = quantile(df_labelled$Maximumverticalheight, c(0.05, 0.95),na.rm = TRUE)) +
  ylab("Maximum height") +
  theme(legend.position = "none",
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
        axis.title.y = element_text(size=16),
        axis.ticks.x = element_blank(),
        #axis.ticks.y = element_blank()
  )#   + 
  #annotate("text", x = -.5, y = 1.5, label = "(b)", size = 8) +
  #coord_cartesian(ylim = c(-1.5, 1.75), xlim = c(1, 8), clip = "off") +
  #theme(plot.margin = unit(c(1,1,1,3), "lines"))

b1

b2 <- ggplot(df_labelled, aes(x=robust_group, y=flowerSize, fill=robust_group)) + 
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
        axis.text.y = element_text(size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        axis.ticks.x = element_blank(),
        #axis.ticks.y = element_blank()
  )  

b2

b1 / b2

# ggsave("figures/10_robust_boxplots.png",width=15,height=10)

####
## Figure 5c: Stacked barplots ----
####

# reset margins
par(mar=c(3,3,3,3))

# add group size to robust group label
for (i in 1:length(unique(df_labelled$robust_group))) {
  df_labelled$robust_group[df_labelled$robust_group %in% sort(unique(df_labelled$robust_group))[i]] <-
    paste(
      sort(unique(df_labelled$robust_group))[i],
      " (n = ",
      table(df_labelled$robust_group)[i],
      ")",
      sep = ""
    )
  
}

# make as factor for grouping
df_labelled$robust_group<-as.factor(df_labelled$robust_group)

# qualitative only
facts <- unlist(lapply(df_labelled, is.factor))
df_temp<-df_labelled[ , facts]

# change table to long form and count combinations
df_temp_melt<-reshape2::melt(df_temp,id.vars="robust_group")
df_temp_melt_counts <- df_temp_melt %>% group_by(robust_group,variable,value) %>% summarise(count=n())

# add new column to remove text labels if counts are <5
df_temp_melt_counts$label<-df_temp_melt_counts$value
df_temp_melt_counts$label[df_temp_melt_counts$count<10]<-NA

#NOT RUN: make new column for text size
#df_temp_melt_counts$text_size<-df_temp_melt_counts$count^(1/2)

# palette from scatterplot
cols <- c(ochre_pal("healthy_reef")(7)[1:4],"thistle",ochre_pal("healthy_reef")(7)[7])

# theme
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

# set NA labels to blank
# df_temp_melt_counts$label[is.na(df_temp_melt_counts$label)]<-""

# remove unwanted traits for plot
df_temp_melt_counts <- df_temp_melt_counts[!df_temp_melt_counts$variable%in%"Climbing",]
df_temp_melt_counts <- df_temp_melt_counts[!df_temp_melt_counts$variable%in%"Aquatic",]

df_temp_melt_counts$variable <- factor(df_temp_melt_counts$variable)

#reorder factors in for plotting
df_temp_melt_counts$variable <- factor(df_temp_melt_counts$variable,
                                       levels = rev(c("SexualSystem",
                                                      "FlowerSex",
                                                      "Mating",
                                                      "FlowerSymmetry",
                                                      "Showiness",
                                                      "FloralReward",
                                                      "OvaryPosition",
                                                      "Pollination",
                                                      "DispersalMode",
                                                      "DispersalDist",
                                                      "Lifespan",
                                                      "Woodiness")))

## Robust group 1

#subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$robust_group=="pam_robust_1 (n = 71)",]

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
  my_theme() +
  theme(
    # axis.text.y = element_blank()
  )

p1

## Robust group 2

#subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$robust_group=="pam_robust_2 (n = 56)",]

#get palette based on max counts
pal<-colorRampPalette(c("white",cols[2]))(max(rob_df$count))

#make new column for colours based on count
rob_df$pal <- pal[rob_df$count]
rob_df$pal[is.na(rob_df$value)]<-"grey20"

p2 <- ggplot(rob_df, aes(x = variable, y = count, fill=pal, alpha=0.98)) +
  geom_bar(position="stack", stat="identity",col="black") +
  scale_fill_identity() + 
  geom_text(aes(size = count, label = label), position = position_stack(vjust = 0.5)) + 
  coord_flip() +
  my_theme() +
  theme(
    axis.text.y = element_blank()
  )

p2

## Robust group 3

#subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$robust_group=="pam_robust_3 (n = 53)",]

#get palette based on max counts
pal<-colorRampPalette(c("white",cols[3]))(max(rob_df$count))

#make new column for colours based on count
rob_df$pal <- pal[rob_df$count]
rob_df$pal[is.na(rob_df$value)]<-"grey20"

p3 <- ggplot(rob_df, aes(x = variable, y = count, fill=pal, alpha=0.98)) +
  geom_bar(position="stack", stat="identity",col="black") +
  scale_fill_identity() + 
  geom_text(aes(size = count, label = label), position = position_stack(vjust = 0.5)) + 
  coord_flip() +
  my_theme() + 
  theme(
    # axis.text.y = element_blank()
  )

p3

## Robust group 4

#subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$robust_group=="pam_robust_4 (n = 48)",]

#get palette based on max counts
pal<-colorRampPalette(c("white",cols[4]))(max(rob_df$count))

#make new column for colours based on count
rob_df$pal <- pal[rob_df$count]
rob_df$pal[is.na(rob_df$value)]<-"grey20"

p4 <- ggplot(rob_df, aes(x = variable, y = count, fill=pal, alpha=0.98)) +
  geom_bar(position="stack", stat="identity",col="black") +
  scale_fill_identity() + 
  scale_y_continuous(breaks=c(0, 10, 20, 30, 40)) +
  geom_text(aes(size = count, label = label), position = position_stack(vjust = 0.5)) + 
  coord_flip() +
  my_theme() +
  theme(
    axis.text.y = element_blank()
  )

p4

## Robust group 5

#subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$robust_group=="pam_robust_5 (n = 37)",]

#get palette based on max counts
pal<-colorRampPalette(c("white",cols[5]))(max(rob_df$count))

#make new column for colours based on count
rob_df$pal <- pal[rob_df$count]
rob_df$pal[is.na(rob_df$value)]<-"grey20"

p5 <- ggplot(rob_df, aes(x = variable, y = count, fill=pal, alpha=0.98)) +
  geom_bar(position="stack", stat="identity",col="black") +
  scale_fill_identity() + 
  geom_text(aes(size = count, label = label), position = position_stack(vjust = 0.5)) + 
  coord_flip()+
  my_theme() +
  theme(
   # axis.text.y = element_blank(),
    axis.title.x = element_text(size=20)
  )

p5

## Robust group 6

#subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$robust_group=="pam_robust_6 (n = 28)",]

#get palette based on max counts
pal<-colorRampPalette(c("white",cols[6]))(max(rob_df$count))

#make new column for colours based on count
rob_df$pal <- pal[rob_df$count]
rob_df$pal[is.na(rob_df$value)]<-"grey20"

p6 <- ggplot(rob_df, aes(x = variable, y = count, fill=pal, alpha=0.98)) +
  geom_bar(position="stack", stat="identity",col="black") +
  scale_fill_identity() + 
  geom_text(aes(size = count, label = label), position = position_stack(vjust = 0.5)) + 
  coord_flip() +
  my_theme() +
  theme(
    axis.text.y = element_blank(),
    axis.title.x = element_text(size=20)
  )

p6

# NA group

# Get name of NA group (depends on number of robust groups kept)
na_name <- unique(grep("NA",df_temp_melt_counts$robust_group,value = TRUE))

#subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$robust_group==na_name,]

#get palette based on max counts
pal<-colorRampPalette(c("white","grey"))(max(rob_df$count))

#make new column for colours based on count
rob_df$pal <- pal[rob_df$count]
rob_df$pal[is.na(rob_df$value)]<-"grey20"

pNA <- ggplot(rob_df, aes(x = variable, y = count, fill=pal, alpha=0.98)) +
  geom_bar(position="stack", stat="identity",col="black") +
  scale_fill_identity() + 
  geom_text(aes(size = count, label = label), position = position_stack(vjust = 0.5)) + 
  #ggtitle("No Robust Group") + 
  coord_flip() +
  ylab("Count") +
  my_theme() +
  theme(
    # axis.text.y = element_blank(),
    axis.title.x = element_text(size=20)
  )

pNA


(p1 + p2) / (p3 + p4) | (p5 + p6) / (pNA + plot_spacer())

# ggsave("figures/10_robust_stacked_barplots.png",width=20,height=15)

### Assemble Figure 5 ----
patch <- (s1 / b1 / b2 ) + plot_layout(heights=c(4, 1, 1)) | (p1 + p2) / (p3 + p4) / (p5 + p6) / (pNA + plot_spacer())

patch + plot_annotation(tag_levels = 'a',tag_prefix="(",tag_suffix=")") & theme(plot.tag = element_text(size = 14))

ggsave("figures/figure_5_robust_groups.png",width=25,height=20)

####
## Figure S8: Qualitative and quantitative values for PAM clusters ----
####

#### 
### Quantitative trait boxplots per cluster ----
#### 

#check names
rownames(df)==rownames(clust.num.k.2.7.df)

#add label to group
df_labelled<-cbind(df,clust.num.k.2.7.df$`3clusters`)

#change colname for label
colnames(df_labelled)[length(colnames(df_labelled))]<-"cluster"

# palette from scatterplot
cols <- harrypotter::hp(3,option="ronweasley2")

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
        # axis.ticks.y = element_blank()
  )#   + 
# annotate("text", x = -.5, y = 1.5, label = "(b)", size = 8) +
# coord_cartesian(ylim = c(-1.5, 1.75), xlim = c(1, 8), clip = "off") +
# theme(plot.margin = unit(c(1,1,1,3), "lines"))

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
        # axis.ticks.y = element_blank()
  )  

b2

b1 / b2

####
### Qualitative trait stacked barplots per cluster ----
####

# add group size to  group label
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

# make as factor for grouping
df_labelled$cluster<-as.factor(df_labelled$cluster)

# qualitative only
facts <- unlist(lapply(df_labelled, is.factor))
df_temp<-df_labelled[ , facts]

# change table to long form and count combinations
df_temp_melt<-reshape2::melt(df_temp,id.vars="cluster")
df_temp_melt_counts <- df_temp_melt %>% group_by(cluster,variable,value) %>% summarise(count=n())

# add new column to remove text labels if counts are <5
df_temp_melt_counts$label<-df_temp_melt_counts$value
df_temp_melt_counts$label[df_temp_melt_counts$count<3]<-NA

# theme
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


# Cluster 1

# subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$cluster==levels(df_temp_melt_counts$cluster)[1],]

# get palette based on max counts
pal1<-colorRampPalette(c("white",cols[1]))(max(rob_df$count))

# make new column for colours based on count
rob_df$pal <- pal1[rob_df$count]
rob_df$pal[is.na(rob_df$value)]<-"grey20"

p1 <- ggplot(rob_df, aes(x = variable, y = count, fill=pal, alpha=0.98)) +
  geom_bar(position="stack", stat="identity",col="black") +
  scale_fill_identity() + 
  geom_text(aes(size = count, label = label), position = position_stack(vjust = 0.5)) + 
  coord_flip() +
  my_theme()

p1


# Cluster 2

# subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$cluster==levels(df_temp_melt_counts$cluster)[2],]

# get palette based on max counts
pal2<-colorRampPalette(c("white",cols[2]))(max(rob_df$count))

# make new column for colours based on count
rob_df$pal <- pal2[rob_df$count]
rob_df$pal[is.na(rob_df$value)]<-"grey20"

p2 <- ggplot(rob_df, aes(x = variable, y = count, fill=pal, alpha=0.98)) +
  geom_bar(position="stack", stat="identity",col="black") +
  scale_fill_identity() + 
  geom_text(aes(size = count, label = label), position = position_stack(vjust = 0.5)) + 
  coord_flip() +
  my_theme()

p2


# Cluster 3

# subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$cluster==levels(df_temp_melt_counts$cluster)[3],]

# get palette based on max counts
pal3<-colorRampPalette(c("white",cols[3]))(max(rob_df$count))

# make new column for colours based on count
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

# assemble figure

patch <- ( p1 + p2 + p3 ) / ( b1 + b2 ) + plot_layout(heights=c(2, 1))

patch + plot_annotation(tag_levels = 'a',tag_prefix="(",tag_suffix=")") & theme(plot.tag = element_text(size = 14))

ggsave("figures/figure_S8_pam_clusters.png",width=20,height=15)
