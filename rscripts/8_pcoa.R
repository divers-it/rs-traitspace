rm(list=ls())

#load libraries
library(dplyr)
library(ggplot2)
library(cluster)
library(patchwork)

#read in imputed data set
df<-read.csv("outputs/imputed_with_phylo.csv",row.names=1,stringsAsFactors = TRUE)

#dissimilarity matrix calculation
gower_df_no_miss <- daisy(df,
                          metric = "gower" )
summary(gower_df_no_miss)

#load data set
df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

#dissimilarity matrix calculation
gower_df <- daisy(df,
                  metric = "gower" )
summary(gower_df)

#check names
labels(gower_df)==gsub("_"," ",labels(gower_df_no_miss))

#compare pairwaise distances of matrices with missing data and with imputed
png("figures/8_scatterplot_dist_missing_vs_imputed.png",width = 500,height = 500)
plot(gower_df,gower_df_no_miss,xlim=c(0,1),ylim=c(0,1)) + abline(0,1,lty=2,col="red")
dev.off()

#make into distance object
dataset_dist <- stats::as.dist(gower_df)

#run PCoA on distance matrix
dataset_pcoa <- ape::pcoa(dataset_dist)

#Recalculate relative eigenvalues by removing negative eigenvalues as in Mouillot et al.  
ev_pcoa <- dataset_pcoa$values$Eigenvalues
ev_pcoa_g0 <- ev_pcoa[ev_pcoa>0]
rel_ev_pcoa_g0 <- ev_pcoa_g0/sum(ev_pcoa_g0)

#check names
rownames(dataset_pcoa$vectors)==rownames(df)

#missing data per row (species) for plot
missDat<-rowSums(apply(is.na(df),2,as.numeric))/ncol(df)

#plot PCOA points on first two axes coloured by missing data
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    aes(color=missDat),
    shape=16,
    alpha=0.75,
    size=2,
    stroke = 0.5
  ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

#save plot
ggsave("figures/8_scatterplot_pcoa_missing.png",
       width = 20,
       height = 15,
       units = 'cm')

#Make data frame of first 9 relative eigenvalues
eig_df<-data.frame(c(1:9),rel_ev_pcoa_g0[1:9])
colnames(eig_df)<-c("pcoa_axis","relative_eigenvalue")
eig_df$pcoa_axis<-as.character(eig_df$pcoa_axis)

#plot barplot of first 9 relative eigenvalues
ggplot(eig_df, aes(x=pcoa_axis, y=relative_eigenvalue)) + 
  geom_bar(stat = "identity")
ggsave("figures/8_barplot_relative_eigenvalues_pcoa.png")

#plot points on first two axes with point style changed by two variables
#reproductive systems (color) and woodiness (shape)
p1 <- ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2, fill = as.factor(df$SexualSystem), shape=as.factor(df$Woodiness))) +
  geom_point(
    color="black",
#    shape=21,
    alpha=0.5,
    size=4,
    stroke = 0.5
  ) + 
  scale_shape_manual(values=c(21,22,23)) +
  guides(fill = guide_legend(override.aes = list(shape = 24) ) ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2))) +
  theme(legend.position="bottom",
        legend.title = element_blank(),
        legend.box="vertical", 
        legend.margin=margin())

# mating system (colour) and woodiness (shape)
p2 <- ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2, fill = as.factor(df$Mating), shape=as.factor(df$Woodiness))) +
  geom_point(
    color="black",
#    shape=21,
    alpha=0.5,
    size=4,
    stroke = 0.5
  ) + 
  scale_shape_manual(values=c(21,22,23)) +
  guides(fill = guide_legend(override.aes = list(shape = 24) ) ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2))) +
  theme(legend.position="bottom",
        legend.title = element_blank(),
        legend.box="vertical", 
        legend.margin=margin())

#combine plots
p1 | p2

#save
ggsave("figures/8_scatterplot_pcoa_coloured_by_traits.png",
       width = 30,
       height = 15,
       units = 'cm')

#PCoA scatterplot with density polygons
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2)) + 
  stat_density_2d(aes(fill = ..level..), geom = "polygon", colour="white") +
  geom_point(
    aes(color = as.factor(df$Pollination), shape=as.factor(df$Woodiness)),
    #    shape=21,
    alpha=0.5,
    size=4,
    stroke = 0.5
  )

ggsave("figures/8_scatterplot_pcoa_density_polygon.png",
       width = 20,
       height = 15,
       units = 'cm')

#Plot density raster with PCoA scatterplot points overlain
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2)) + 
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=3, direction=1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    legend.position='none'
  ) +
  geom_point(
    aes(color = as.factor(df$Mating), shape=as.factor(df$Woodiness)),
    #    shape=21,
    alpha=0.5,
    size=4,
    stroke = 0.5
  )

# ---- Taxonomy ----

#load taxonomy
tax<-readRDS(file = here::here("outputs/taxonomy.rds"))

#check row order
rownames(df)==rownames(tax)

#empty column for common orders
tax$order_common<-NA

#top 10 most common orders in data
order_common<-names(sort(table(tax$order),decreasing = T)[1:8])

for(i in 1:length(order_common)){
  tax$order_common[grep(order_common[i],tax$order)]<-order_common[i]
}

#plot PCOA points on first two axes coloured by order
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    aes(color=tax$order_common),
    shape=16,
    alpha=0.75,
    size=3,
    stroke = 0.5
  ) + 
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

ggsave("figures/8_scatterplot_pcoa_taxonomy.png",
       width = 20,
       height = 15,
       units = 'cm')
