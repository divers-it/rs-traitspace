rm(list = ls())

#load packages
library(densityClust)
library(cluster)
library(ggplot2)
library(networkD3)
library(dplyr)
library(gridExtra)
library(data.table)
library(ggrepel)

#load formatted data
df <- readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

#calc gower distance
gower_df <- daisy(df,
                  metric = "gower" )
summary(gower_df)

#convert to distance matrix
dataset_dist <- stats::as.dist(gower_df)

#run PCOA
dataset_pcoa <- ape::pcoa(dataset_dist)

#Recalculate relative eigenvalues by removing negative eigenvalues as in Mouillot et al.  
ev_pcoa <- dataset_pcoa$values$Eigenvalues
ev_pcoa_g0 <- ev_pcoa[ev_pcoa>0]
rel_ev_pcoa_g0 <- ev_pcoa_g0/sum(ev_pcoa_g0)

#run initial density clustering
protClust <- densityClust(gower_df, gaussian=TRUE)

# Inspect clustering attributes to define thresholds
plot(protClust) 

#set values of rho and delta to choose numbers of clusters
#k=4
protClust <- findClusters(protClust, rho = 30,delta = 0.07,verbose = FALSE,plot = TRUE)

#species that represent cluster peaks
protClust$peaks

#check order
rownames(dataset_pcoa$vectors)==names(protClust$halo)

#make df for plotting
pcoa_df<-data.frame(dataset_pcoa$vectors)

#add names to data to be plotted
pcoa_df$names<-NA

#put in name of medoids
for(i in 1:length(protClust$peaks)){
  pcoa_df$names[protClust$peaks[i]]<-names(protClust$peaks[i])
}

#plot points on first two axes, coloured by cluster
ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, fill = as.factor(protClust$clusters))) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.5,
    size=3,
    stroke = 0.5
  ) + geom_text_repel(aes(label = pcoa_df$names, colour = as.factor(protClust$clusters)),
                        size = 3.5)  + 
  stat_ellipse(geom = "polygon",
               aes(fill = as.factor(protClust$clusters)), 
               alpha = 0.25) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2))) +
  theme(legend.position="none")

#save image
ggsave("figures/11.4_scatterplot_pcoa_dens_k4_coloured_by_cluster.png",width = 10,height=10)

#view cluster membership of species
split(rownames(df), protClust$clusters)

#vectors of values of stats for decision graphs that yield 2 to 7 clusters
rho_v<-c(30,30,30,15,8,8)
delta_v<-c(0.6,0.15,0.1,0.1,0.08,0.07)

#make df with different numbers of clusters
for(i in 1:length(rho_v)){
  if(i == 1){
    protClust <- findClusters(protClust, rho = rho_v[i],delta = delta_v[i],verbose = FALSE,plot = FALSE)
    dens_df<-protClust$clusters
  } else {
    protClust <- findClusters(protClust, rho = rho_v[i],delta = delta_v[i],verbose = FALSE,plot = FALSE)
    dens_df<-cbind(dens_df,protClust$clusters)
  }
}

colnames(dens_df)<-c("2clusters",
                    "3clusters",
                    "4clusters",
                    "5clusters",
                    "6clusters",
                    "7clusters")
clust.num.k.2.7.df <-as.data.frame(dens_df)

#check names
names(protClust$halo)==rownames(df)

#add rownames
rownames(clust.num.k.2.7.df)<-rownames(df)
head(clust.num.k.2.7.df)

###
# ---- Sankey plot ----
###

#from: https://r-graph-gallery.com/321-introduction-to-interactive-sankey-diagram-2.html

#table of different k values (2-7) adding label
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
rownames(clust.num.k.2.7.df)<-rownames(df)
head(clust.num.k.2.7.df)

#save RDS for downstream use
saveRDS(clust.num.k.2.7.df, file = here::here("outputs/11.4_clust_num_k_2_7_density.rds"))

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

# Make the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE)
p

#save as html
saveNetwork(p, "figures/11.4_sankey_dens.html")

###
# ---- Robust combinations ----
###

#make data frame of combo frequencies
combos <- as.data.frame(table(clust.num.k.2.7.df))

#remove no existant combos
combos <- combos[combos$Freq > 0, ]

#order
combos <- combos[order(combos$Freq, decreasing = T), ]

#change to strings
combos <-data.frame(lapply(combos, as.character), stringsAsFactors = FALSE)
plot(combos$Freq)

#empty list
robust<-list()

#empty vector CHECK ROWNAMES CORRECT
robust_vect_dens<-rep(NA,length(rownames(df)))
names(robust_vect_dens)<-rownames(df)

#loop through ordered table to extract robust groups
for(i in 1:length(combos$Freq[as.numeric(combos$Freq)>20])){
  foo<-rownames(clust.num.k.2.7.df[clust.num.k.2.7.df[, 1] == combos[i, 1] & 
                                                clust.num.k.2.7.df[, 2] == combos[i, 2] &
                                                clust.num.k.2.7.df[, 3] == combos[i, 3] &
                                                clust.num.k.2.7.df[, 4] == combos[i, 4] &
                                                clust.num.k.2.7.df[, 5] == combos[i, 5] &
                                                clust.num.k.2.7.df[, 6] == combos[i, 6],])
  
  robust[[i]]<-rownames(df)[foo]
  
  robust_vect_dens[foo]<-i
  
}

#robust groups
head(robust)

#complete vector of robust groups and non-robust 
robust_vect_dens_full<-robust_vect_dens
saveRDS(robust_vect_dens_full, file = here::here("outputs/11.4_robust_vect_dens_full.rds"))

#remove species not in robust groups
robust_vect_dens<-na.omit(robust_vect_dens)

#empty matrix
rob_mat<-matrix(nrow = length(unique(robust_vect_dens)), ncol=length(df[1,]))

#empty matrix
rob_mat_names<-matrix(nrow = length(unique(robust_vect_dens)), ncol=length(df[1,]))

#loop through different robust groups
for(i in 1:length(unique(robust_vect_dens))){
  
  #names of species in robust group
  grp<-names(robust_vect_dens)[robust_vect_dens==i]
  
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
rownames(rob_mat)<-paste("robust",c(1:length(unique(robust_vect_dens))),sep="")
colnames(rob_mat)<-colnames(df)
rob_mat

rownames(rob_mat_names)<-paste("robust",c(1:length(unique(robust_vect_dens))),sep="")
colnames(rob_mat_names)<-colnames(df)
rob_mat_names

###
# ---- Boxplots and stacked barplots for robust groups ----
###

#make label
robust_group<-paste("dens_robust_",robust_vect_dens_full,sep="")

#check names
rownames(df)==names(robust_vect_dens_full)

#add label to group
df_labelled<-cbind(df,as.factor(robust_group))

#change colname for label
colnames(df_labelled)[length(colnames(df_labelled))]<-"cluster"

#empty list for plots
plot_list <- list()

#make plots
for(i in 1:(length(colnames(df_labelled))-1)){
  
  #for quantitative
  if(is.numeric(df_labelled[,i])){
    plot_list[[i]]<-ggplot(df_labelled, aes(x=cluster, y=!!as.name(colnames(df_labelled)[i]), fill=cluster)) + geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.1)) + theme(legend.position = "none",axis.text.x = element_text(angle = 90),axis.title.x = element_blank(),plot.margin = unit(c(1,1,1,1), "cm"))  
  } else {   #for qualitative
    plot_list[[i]]<-ggplot(df_labelled, aes(x=cluster, fill = !!as.name(colnames(df_labelled)[i]))) +
      geom_bar(stat="count") + theme(axis.text.x = element_text(angle = 90),axis.title.x = element_blank(),plot.margin = unit(c(1,1,1,1), "cm"))
  }
  
}

pdf("figures/11.4_boxplots_stacked_barplots_dens_clusters_by_trait.pdf",width = 15,height = 15)

print(grid.arrange(grobs=plot_list[1:4],ncol=2,nrow=2))
print(grid.arrange(grobs=plot_list[5:8],ncol=2,nrow=2))
print(grid.arrange(grobs=plot_list[9:12],ncol=2,nrow=2))
print(grid.arrange(grobs=plot_list[13:16],ncol=2,nrow=2))
print(grid.arrange(grobs=plot_list[17:19],ncol=2,nrow=2))

dev.off()

###
# Plot qualitative stats of robust groups
###

#reset margins
par(mar=c(3,3,3,3))

#make cluster as character
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

#NOT RUN: make new column for text size
#df_temp_melt_counts$text_size<-df_temp_melt_counts$count^(1/2)

#plot stacked barplots per robust group for each qualitative trait, with labels
ggplot(df_temp_melt_counts, aes(variable, count, fill = value)) +
  geom_col(position = 'stack') + facet_wrap(. ~ cluster, scales = "free")  + theme(
    legend.position = "none",
    axis.text.x = element_text( vjust = 0.5, hjust=1),
    axis.title.x = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) + geom_text(aes(size = count,label = label),
                                position = position_stack(vjust = .5)) + coord_flip()

#save image
ggsave("figures/11.4_stacked_barplots_dens_traits_by_cluster.png",width=15,height=10)

