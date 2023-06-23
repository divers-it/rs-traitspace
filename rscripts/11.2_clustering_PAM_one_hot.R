rm(list = ls())
library(dplyr)
library(gridExtra)
library(Rtsne)
library(ggplot2)

#load formatted data
df<-readRDS(file = here::here("outputs/df_filt_trans_one_hot.rds"))

#numeric columns only
nums <- unlist(lapply(df, is.numeric))
facts <- unlist(lapply(df, is.factor))
df<-cbind(df[ , nums],df[ , facts])

library(cluster)
gower_df <- daisy(df,
                  metric = "gower" )

summary(gower_df)

#introduced NAs - need to solve better
#gower_df[is.na(gower_df)]<-0

#Silhouette Width to select the optimal number of clusters
#The silhouette width is one of the very popular choices when it comes to selecting the optimal number of clusters. It measures the similarity of each point to its cluster, and compares that to the similarity of the point with the closest neighboring cluster. This metric ranges between -1 to 1, where a higher value implies better similarity of the points to their clusters. Therefore, a higher value of the Silhouette Width is desirable. We calculate this metric for a range of cluster numbers and find where it is maximized. The following code shows the implementation in R:
silhouette <- c()
silhouette = c(silhouette, NA)
for(i in 2:10){
  pam_clusters = pam(as.matrix(gower_df),
                     diss = TRUE,
                     k = i)
  silhouette = c(silhouette ,pam_clusters$silinfo$avg.width)
}

plot(1:10, silhouette,
     xlab = "Clusters",
     ylab = "Silhouette Width")
lines(1:10, silhouette)


#make df with each value of k
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

saveRDS(clust.num.k.2.7.df, file = here::here("outputs/clust_num_k_2_7_pam_one_hot.rds"))

#construct a PAM model with X clusters, and try to interpret the behavior of these clusters with the help of the medoids.
pam.gower = pam(gower_df, diss = TRUE, k =3)
df[pam.gower$medoids, ]
write.csv(df[pam.gower$medoids, ], "outputs/pam_medoids_k3_one_hot.csv")

#To dig deeper into the characteristics of each cluster, we find the summary stats.
pam_summary <- df %>%
  mutate(cluster = pam.gower$clustering) %>%
  group_by(cluster) %>%
  do(cluster_summary = summary(.))

#summary stats of cluster 1
pam_summary$cluster_summary[[1]]

#set palette
library(RColorBrewer)
palette(brewer.pal(6,"Dark2"))

#the t-SNE or the t-Distributed Stochastic Neighbor Embedding technique
#alternative to PCOA
tsne_object <- Rtsne(gower_df, is_distance = TRUE)
tsne_df <- tsne_object$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(pam.gower$clustering))

#add rownames (NOT SURE IF MATCH)
tsne_df$names<-rownames(df)

#prune down species name to make readable
inds <- round ( runif(320, 1, length(tsne_df$names)) )
tsne_df$names[inds]<-NA
#tsne_df$names[sample(seq_along(tsne_df$names), 280, replace = FALSE)] <- NA

#put in medoids
for(i in 1:length(pam.gower$medoids)){
  ind<-grep(pam.gower$medoids[i],rownames(df))
  tsne_df$names[ind]<-rownames(df)[ind]
}

#plot points on first two axes, coloured by cluster
library(ggrepel)

ggplot(tsne_df, aes(x = X, y = Y, fill = as.factor(cluster))) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.5,
    size=3,
    stroke = 0.5
  ) + 
  geom_text_repel(aes(label = tsne_df$names, colour = as.factor(cluster)),
                  size = 3.5)  + 
  stat_ellipse(geom = "polygon",
               aes(fill = cluster), 
               alpha = 0.25) +
  xlab("t-SNE Axis 1") +
  ylab("t-SNE Axis 2")

ggsave("figures/scatter_tsne_pam_clusters_one_hot.png",width=12,height=10)

#CHECK ROWNAMES ARE IN CORRECT ORDER
clust_df<-data.frame(rownames(df),tsne_df$cluster,row.names=NULL)
colnames(clust_df)<-c("Species","Cluster")
clust_df<-clust_df[order(clust_df$Species),]
head(clust_df)

clust_df$Species<-gsub(" ","_",clust_df$Species)

write.csv(clust_df,"outputs/pam_clustering_one_hot.csv")

#save.image("outputs/pam_clustering.Rdata")

#pcoa
dataset_dist <- stats::as.dist(gower_df)
dataset_pcoa <- ape::pcoa(dataset_dist)

#plot points on first two axes, coloured by cluster
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2, fill = as.factor(pam.gower$clustering))) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.5,
    size=3,
    stroke = 0.5
  ) + 
  stat_ellipse(geom = "polygon",
               aes(fill =  as.factor(pam.gower$clustering)), 
               alpha = 0.25) +
  xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))

ggsave("figures/scatter_pcoa_pam_clusters_one_hot.png",width=12,height=10)


####
# Sankey plot
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

#fix rownames again
rownames(clust.num.k.2.7.df)<-names(pam.gower$clustering)

saveRDS(clust.num.k.2.7.df, file = here::here("outputs/clust_num_k_2_7_pam_one_hot.rds"))

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

saveNetwork(p, "figures/sankey_pam_one_hot.html")

###
# Robust combinations
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

#empty vector
robust_vect_pam<-rep(NA,length(rownames(dataset_pcoa$vectors)))
names(robust_vect_pam)<-rownames(dataset_pcoa$vectors)


#loop through ordered table to extract robust groups
for(i in 1:length(combos$Freq[as.numeric(combos$Freq)>20])){
  foo<-as.numeric(rownames(clust.num.k.2.7.df[clust.num.k.2.7.df[, 1] == combos[i, 1] & 
                                                clust.num.k.2.7.df[, 2] == combos[i, 2] &
                                                clust.num.k.2.7.df[, 3] == combos[i, 3] &
                                                clust.num.k.2.7.df[, 4] == combos[i, 4] &
                                                clust.num.k.2.7.df[, 5] == combos[i, 5] &
                                                clust.num.k.2.7.df[, 6] == combos[i, 6],]))
  
  robust[[i]]<-rownames(dataset_pcoa$vectors)[foo]
  
  robust_vect_pam[foo]<-i
  
}

#robust groups
robust

#complete vector of robust groups and non-robust 
robust_vect_pam_full<-robust_vect_pam
saveRDS(robust_vect_pam_full, file = here::here("outputs/robust_vect_pam_full_one_hot.rds"))


#remove species not in robust groups
robust_vect_pam<-na.omit(robust_vect_pam)

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

#Plot robust groups
tsne_df_robust<-cbind(tsne_df,robust_vect_pam_full)


#check order
rownames(tsne_df_robust)==rownames(clust.num.k.2.7.df)

#Plot robust groups
#plot points on first two axes, coloured by cluster
ggplot(
  data.frame(tsne_df_robust),
  aes(
    x = X,
    y = Y,
    col = as.factor(robust_vect_pam_full)
  )
) +
  geom_point(
    aes(shape = as.factor(clust.num.k.2.7.df$`3clusters`)),
    alpha = 0.5,
    size = 3,
    stroke = 0.5
  )

ggsave("figures/scatter_pcoa_pam_robust_one_hot.png",width=12,height=10)

#species that dont belong to robust group
df_not_robust<-df[is.na(robust_vect_pam_full),]
mean(is.na(df_not_robust))

df_robust<-df[!is.na(robust_vect_pam_full),]
mean(is.na(df_robust))

#####
#Boxplots and stacked barplots for robust groups
#####

# library
library(ggplot2)

#make label
robust_group<-paste("pam_robust_",robust_vect_pam_full,sep="")

#add label to group
df_labelled<-cbind(df,robust_group)

#empty list for plots
plot_list <- list()

#make plots
for(i in 1:(length(colnames(df_labelled))-1)){
  
  #for quantitative
  if(is.numeric(df_labelled[,i])){
    plot_list[[i]]<-ggplot(df_labelled, aes(x=robust_group, y=!!as.name(colnames(df_labelled)[i]), fill=robust_group)) + 
      geom_boxplot() + geom_jitter(shape=16, position=position_jitter(0.1)) + theme(legend.position = "none",axis.text.x = element_text(angle = 90),axis.title.x = element_blank(),plot.margin = unit(c(1,1,1,1), "cm"))  
  } else {   #for qualitative
    plot_list[[i]]<-ggplot(df_labelled, aes(x=robust_group, fill = !!as.name(colnames(df_labelled)[i]))) +
      geom_bar(stat="count") + theme(axis.text.x = element_text(angle = 90),axis.title.x = element_blank(),plot.margin = unit(c(1,1,1,1), "cm"))
  }
  
}

pdf("figures/robust_pam_plots_one_hot.pdf",width = 15,height = 15)

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

#add group size to robust group label
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

#make as factor for grouping
df_labelled$robust_group<-as.factor(df_labelled$robust_group)

#qualitative only
facts <- unlist(lapply(df_labelled, is.factor))
df_temp<-df_labelled[ , facts]

#change table to long form and count combinations
df_temp_melt<-data.table::melt(df_temp,id.vars="robust_group")
df_temp_melt_counts <- df_temp_melt %>% group_by(robust_group,variable,value) %>% summarise(count=n())

#add new column to remove text labels if counts are <5
df_temp_melt_counts$label<-df_temp_melt_counts$value
df_temp_melt_counts$label[df_temp_melt_counts$count<3]<-NA

#make new column for text size
#df_temp_melt_counts$text_size<-df_temp_melt_counts$count^(1/2)

#plot stacked barplots per robust group for each qualitative trait, with labels
ggplot(df_temp_melt_counts, aes(variable, count, fill = value)) +
  geom_col(position = 'stack') + facet_wrap(. ~ robust_group, scales = "free")  + theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.title.x = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) + geom_text(aes(label = label),
                size = 2,
                angle = 90,
                position = position_stack(vjust = .5))

ggsave("figures/stacked_barplots_robust_groups_pam_one_hot.pdf",width=15,height=15)

