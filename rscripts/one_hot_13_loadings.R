rm(list = ls())
library(dplyr)
library(cluster)
library(vegan)
library(ggrepel)
library(wesanderson)
library(ggimage)
library(png)
library(patchwork)

# load formatted data
df2<-readRDS(file = here::here("outputs/one_hot_6_df_filt_trans.rds"))

# make factors integers
df2[sapply(df2, is.factor)] <- lapply(df2[sapply(df2, is.factor)],as.integer)

# calculate dissimilarity matrix

gower_df2 <- daisy(df2, metric = "gower" )
summary(gower_df2)

# convert gower df to distance matrix
dataset_dist2 <- stats::as.dist(gower_df2)

####
## Run PCOA ----
####

dataset_pcoa2 <- wcmdscale(d = dataset_dist2, eig = TRUE)

# get scores for calculating eigenvectors
traitvectors12=vegan::scores(envfit(dataset_pcoa2, df2, na.rm = T, add = T),display="vectors")
traitvectors34=vegan::scores(envfit(dataset_pcoa2, df2, na.rm = T, add = T,choices=c(3,4)),display="vectors")
traitd=as.data.frame(traitvectors12)
traitd$Dim3=as.data.frame(traitvectors34)$Dim3
traitd$Dim4=as.data.frame(traitvectors34)$Dim4
traitd$trait=rownames(traitd)

speciesv=vegan::scores(dataset_pcoa2,display="species")
speciesd=as.data.frame(speciesv)

# Recalculate relative eigenvalues by removing negative eigenvalues
ev_pcoa <- dataset_pcoa2$eig
ev_pcoa_g0 <- ev_pcoa[ev_pcoa>0]
rel_ev_pcoa_g0 <- ev_pcoa_g0/sum(ev_pcoa_g0)

# add original traits (multistate) for visualization
if(sum(rownames(df2)==rownames(speciesd))==nrow(df2)){
	df_ord=as.data.frame(c(df2,speciesd))
	rownames(df_ord)=rownames(df2)
} else {
	stop("Row (species) names are not identical")
}

# read in clustering results
clusters_pam_one_hot=readRDS(file = here::here("outputs/one_hot_10_clust_num_k_2_7_pam.rds"))
clusters_pam=readRDS(file = here::here("outputs/10_clust_num_k_2_7_pam.rds"))
clusters_kpro_one_hot=readRDS(file = here::here("outputs/one_hot_10.1_clust_num_k_2_7_kpro.rds"))
clusters_kpro=readRDS(file = here::here("outputs/10.1_clust_num_k_2_7_kpro.rds"))
clusters_ward_one_hot=readRDS(file = here::here("outputs/one_hot_10.2_clust_num_k_2_7_ward.rds"))
clusters_ward=readRDS(file = here::here("outputs/10.2_clust_num_k_2_7_ward.rds"))

# select number of clusters
# divisive.clust = 3
# aggl.clust.c (complete) = 4
# aggl.clust.a (average) = 2
# aggl.clust.w (ward) = 5

# clust.num[clust.num==1]="showy perennials"
# clust.num[clust.num==2]="dioecious trees"
# clust.num[clust.num==3]="small/short-lived herbs"
# clust.num[clust.num==4]="showy trees"
# clust.num[clust.num==5]="showy trees"
# clust.num[clust.num==6]="wind-pollinated sexually monomorphic"

# minimal value for arrow length
minarrow=0.2

# number of clusters to use for plotting
clust.num=data.frame(clust.num=clusters_pam[,"3clusters"])
rownames(clust.num)=rownames(clusters_ward)
df_ord_clust=as.data.frame(c(df_ord,clust.num))

#### 
# Figure 3: Loadings scatterplot without cluster Axes 1 and 2 ----
#### 

# Images found in AJH_DiveRS/trait_circles.svg
# trait images
woody_img <- readPNG(paste(here::here(),"/figures/images/woody1.png",sep=""))
herb_img <- readPNG(paste(here::here(),"/figures/images/herb1.png",sep=""))
dull_img <- readPNG(paste(here::here(),"/figures/images/dull_flower1.png",sep=""))
showy_img <- readPNG(paste(here::here(),"/figures/images/showy_flower1.png",sep=""))

# rename traits
traitd$trait<-gsub("_"," ",traitd$trait)
traitd$trait<-gsub("\\."," ",traitd$trait)

# remove unimportant traits

# hypotenuse 
traitd$hypo <- abs(traitd$Dim1)^2+abs(traitd$Dim2)^2
hist(traitd$hypo)
traitd_labels <- traitd[traitd$hypo > 0.15,]

# plot
ggplot() + 
  geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2),shape=21,fill="grey",alpha=0.4,size=6) + 
  geom_segment(data=traitd[sqrt(traitd$Dim1^2+traitd$Dim2^2)>minarrow,],
               aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2), col="grey30",
               # aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2,col=signal),
               # aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2,col=log(transition_rates)),
               arrow=arrow(length = unit(0.4, "cm")),
               alpha=0.7,
               linewidth=0.85,
               lineend='round',
               linejoin='round') +
  geom_text_repel(data=traitd_labels[sqrt(traitd_labels$Dim1^2+traitd_labels$Dim2^2)>minarrow,],
                  aes(x=Dim1/2,y=Dim2/2,label=trait),size=4) +
  scale_color_gradientn(colours = rainbow(3)) +
  xlim(-0.5,0.5) +
  ylim(-0.5,0.5) +
  theme_bw() + theme(
      panel.border = element_blank(),
      # panel.grid.major = element_line(colour = "darkgrey"),
      # panel.grid.minor = element_line(colour = "grey"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = c(0.9, 0.8),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size=14),
      axis.title = element_text(size=18)) +
  xlab(paste("PCoA Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("PCoA Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2))) +
  annotation_raster(dull_img, ymin = 0.34,ymax= 0.49,xmin = -0.4,xmax = -0.25) +
  annotation_raster(herb_img, ymin = 0.375,ymax= 0.525,xmin = 0.225,xmax = 0.375) +
  annotation_raster(woody_img, ymin = -0.5,ymax= -0.35,xmin = -0.475,xmax = -0.325) +
  annotation_raster(showy_img, ymin = -0.4,ymax= -0.25,xmin = 0.26,xmax = 0.41)

ggsave("figures/figure_3_loadings.png",width=12.5,height=12.5)
ggsave("figures/figure_3_loadings.pdf",width=12.5,height=12.5)

#### 
## Figure S3: Loadings scatterplot without cluster Axes 3 and 4 ----
#### 

# remove unimportant traits
traitd$absDim3Dim4<-abs(traitd$Dim3)+abs(traitd$Dim4)
hist(traitd$absDim3Dim4)
traitd_labels <- traitd[traitd$absDim3Dim4>0.55,]

# plot
ggplot() + 
  geom_point(data=df_ord_clust,aes(x=Dim3,y=Dim4),shape=21,fill=wes_palette("Darjeeling1")[3],alpha=0.4,size=6) + 
  geom_segment(data=traitd[sqrt(traitd$Dim3^2+traitd$Dim4^2)>minarrow,],
               aes(x=0,y=0,xend=Dim3/2,yend=Dim4/2),
               arrow=arrow(length = unit(0.4, "cm")),
               col="grey40",
               alpha=0.7,
               linewidth=0.85,
               lineend='round',
               linejoin='round') +
  geom_text_repel(data=traitd_labels[sqrt(traitd_labels$Dim3^2+traitd_labels$Dim4^2)>minarrow,],
                  aes(x=Dim3/2,y=Dim4/2,label=trait),size=6) +
  # xlim(-0.4,0.4) +
  # ylim(-0.4,0.4) +
  theme_bw() + theme(
    panel.border = element_blank(),
    # panel.grid.major = element_line(colour = "darkgrey"),
    # panel.grid.minor = element_line(colour = "grey"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.9, 0.8),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size=14),
    axis.title = element_text(size=16)) +
  xlab(paste("PCoA Axis 3: relative eigenvalue =",round(rel_ev_pcoa_g0[3],2))) +
  ylab(paste("PCoA Axis 4: relative eigenvalue =",round(rel_ev_pcoa_g0[4],2)))

ggsave("figures/figure_S3_loadings_axes_3_4.png",width=12.5,height=12.5)

####
## Plot loadings with clustering on first 4 PCOA axes ----
####

# plot points on first two axes, coloured by cluster
ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(clust.num)))+geom_segment(data=traitd[sqrt(traitd$Dim1^2+traitd$Dim2^2)>minarrow,],aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd[sqrt(traitd$Dim1^2+traitd$Dim2^2)>minarrow,],aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust, geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(clust.num)), alpha = 0.25)+theme_bw() + theme(legend.position="top",legend.title=element_blank()) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

# plot points on second and third axis, coloured by cluster
ggplot() + geom_point(data=df_ord_clust,aes(x=Dim2,y=Dim3,color = as.factor(clust.num)))+geom_segment(data=traitd[sqrt(traitd$Dim2^2+traitd$Dim3^2)>minarrow,],aes(x=0,y=0,xend=Dim2/2,yend=Dim3/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd[sqrt(traitd$Dim2^2+traitd$Dim3^2)>minarrow,],aes(x=Dim2/2,y=Dim3/2,label=trait))+  stat_ellipse(data=df_ord_clust, geom = "polygon", aes(x=Dim2,y=Dim3,fill = as.factor(clust.num)), alpha = 0.25)+theme_bw() + theme(legend.position="none") +
  xlab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2))) +
  ylab(paste("Axis 3: relative eigenvalue =",round(rel_ev_pcoa_g0[3],2)))

# plot points on third and fourth axis, coloured by cluster
ggplot() + geom_point(data=df_ord_clust,aes(x=Dim3,y=Dim4,color = as.factor(clust.num)))+geom_segment(data=traitd[sqrt(traitd$Dim3^2+traitd$Dim4^2)>minarrow,],aes(x=0,y=0,xend=Dim3/2,yend=Dim4/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd[sqrt(traitd$Dim3^2+traitd$Dim4^2)>minarrow,],aes(x=Dim3/2,y=Dim4/2,label=trait))+  stat_ellipse(data=df_ord_clust, geom = "polygon", aes(x=Dim3,y=Dim4,fill = as.factor(clust.num)), alpha = 0.25)+theme_bw() + theme(legend.position="none") +
  xlab(paste("Axis 3: relative eigenvalue =",round(rel_ev_pcoa_g0[3],2))) +
  ylab(paste("Axis 4: relative eigenvalue =",round(rel_ev_pcoa_g0[4],2)))

####
## Compare clustering methods ----
####

### Original encoding ----

numclust="6clusters"

clust.num=data.frame(clust.num=clusters_ward[,numclust])
rownames(clust.num)=rownames(clusters_ward)
df_ord_clust=as.data.frame(c(df_ord,clust.num))

p1=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(clust.num)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust, geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(clust.num)), alpha = 0.25) + theme(legend.position="none") + ggtitle("ward") +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

clust.num=data.frame(clust.num=clusters_kpro[,numclust])
rownames(clust.num)=rownames(clusters_kpro)
df_ord_clust=as.data.frame(c(df_ord,clust.num))

p2=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(clust.num)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust, geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(clust.num)), alpha = 0.25) + theme(legend.position="none") + ggtitle("kproto") +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

clust.num=data.frame(clust.num=clusters_pam[,numclust])
rownames(clust.num)=rownames(clusters_pam)
df_ord_clust=as.data.frame(c(df_ord,clust.num))

p3=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(clust.num)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust, geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(clust.num)), alpha = 0.25) + theme(legend.position="none") + ggtitle("pam") +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

p1 + p2 + p3

# ggsave("figures/one_hot_13_scatterplot_pcoa_k6_loadings_by_clustering_method.png",
#        width = 40,
#        height = 20,
#        units = 'cm')

### One-hot based clustering ----

numclust="6clusters"

clust.num=data.frame(clust.num=clusters_ward_one_hot[,numclust])
rownames(clust.num)=rownames(clusters_ward_one_hot)
df_ord_clust=as.data.frame(c(df_ord,clust.num))

p1=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(clust.num)))+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust, geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(clust.num)), alpha = 0.25) + theme(legend.position="none") + ggtitle("ward") +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

clust.num=data.frame(clust.num=clusters_kpro_one_hot[,numclust])
rownames(clust.num)=rownames(clusters_kpro_one_hot)
df_ord_clust=as.data.frame(c(df_ord,clust.num))

p2=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(clust.num)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust, geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(clust.num)), alpha = 0.25) + theme(legend.position="none") + ggtitle("kproto") +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

clust.num=data.frame(clust.num=clusters_pam_one_hot[,numclust])
rownames(clust.num)=rownames(clusters_pam_one_hot)
df_ord_clust=as.data.frame(c(df_ord,clust.num))

p3=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(clust.num)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust, geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(clust.num)), alpha = 0.25) + theme(legend.position="none") + ggtitle("pam") +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

p1 + p2 + p3

# ggsave("figures/one_hot_13_scatterplot_pcoa_k6_loadings_by_one_hot_clustering_method.png",
#        width = 40,
#        height = 20,
#        units = 'cm')

####
## Compare robust groups with loadings ----
####

### Original encoding ----

# read robust groups 
groups_kpro_one_hot=readRDS(file = here::here("outputs/one_hot_10.1_robust_vect_kpro_full.rds"))
groups_kpro=readRDS(file = here::here("outputs/10.1_robust_vect_kpro_full.rds"))
groups_pam_one_hot=readRDS(file = here::here("outputs/one_hot_10_robust_vect_pam_full.rds"))
groups_pam=readRDS(file = here::here("outputs/10_robust_vect_pam_full.rds"))

rownames(df_ord)==rownames(groups_kpro)

df_ord_clust=as.data.frame(c(df_ord,as.data.frame(groups_kpro)))
names(df_ord_clust)[ncol(df_ord_clust)]="group"

p1=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(group)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust[!is.na(df_ord_clust$group),], geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(group)), alpha = 0.25) + theme(legend.position="none") + ggtitle("kproto") +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

df_ord_clust=as.data.frame(c(df_ord,as.data.frame(groups_pam)))
names(df_ord_clust)[ncol(df_ord_clust)]="group"

p2=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(group)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust[!is.na(df_ord_clust$group),], geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(group)), alpha = 0.25) + theme(legend.position="none") + ggtitle("pam") +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

p1 + p2

# ggsave("figures/one_hot_13_scatterplot_pcoa_robust_loadings_by_clustering_method.png",
#        width = 40,
#        height = 20,
#        units = 'cm')

### One-hot ----

df_ord_clust=as.data.frame(c(df_ord,as.data.frame(groups_kpro_one_hot)))
names(df_ord_clust)[ncol(df_ord_clust)]="group"

p1=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(group)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust[!is.na(df_ord_clust$group),], geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(group)), alpha = 0.25) + theme(legend.position="none") + ggtitle("kproto") +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

df_ord_clust=as.data.frame(c(df_ord,as.data.frame(groups_pam_one_hot)))
names(df_ord_clust)[ncol(df_ord_clust)]="group"

p2=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(group)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust[!is.na(df_ord_clust$group),], geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(group)), alpha = 0.25) + theme(legend.position="none") + ggtitle("pam") +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

p1 + p2

# ggsave("figures/one_hot_13_scatterplot_pcoa_robust_loadings_by_one_hot_clustering_method.png",
#        width = 40,
#        height = 20,
#        units = 'cm')

####
## Big table for checking characters ----
####

df <- readRDS(file = here::here("outputs/5_df_filt.rds"))
tax <- readRDS(file = here::here("outputs/taxonomy.rds"))

if(sum(rownames(df)==rownames(tax))==nrow(df)){
	df_tax=as.data.frame(c(tax,df))
	rownames(df_tax)=rownames(df)
} else {
	stop("Row (species) names are not identical")
}


if(sum(rownames(df_tax)==rownames(as.data.frame(groups_pam)))==nrow(df_tax)){
	df_tax_clust=as.data.frame(c(df_tax,as.data.frame(groups_pam)))
	rownames(df_tax_clust)=rownames(df_tax)
} else {
	stop("Row (species) names are not identical")
}
if(sum(rownames(df_tax_clust)==rownames(as.data.frame(groups_pam_one_hot)))==nrow(df_tax_clust)){
	df_tax_clust=as.data.frame(c(df_tax_clust,as.data.frame(groups_pam_one_hot)))
	rownames(df_tax_clust)=rownames(df_tax)
} else {
	stop("Row (species) names are not identical")
}
if(sum(rownames(df_tax_clust)==rownames(as.data.frame(groups_kpro)))==nrow(df_tax_clust)){
	df_tax_clust=as.data.frame(c(df_tax_clust,as.data.frame(groups_kpro)))
	rownames(df_tax_clust)=rownames(df_tax)
} else {
	stop("Row (species) names are not identical")
}
if(sum(rownames(df_tax_clust)==rownames(as.data.frame(groups_kpro_one_hot)))==nrow(df_tax_clust)){
	df_tax_clust=as.data.frame(c(df_tax_clust,as.data.frame(groups_kpro_one_hot)))
	rownames(df_tax_clust)=rownames(df_tax)
} else {
	stop("Row (species) names are not identical")
}

write.csv(df_tax_clust,file="outputs/13_robust_groups_with_trait_data.csv")
