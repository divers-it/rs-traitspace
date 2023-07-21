rm(list = ls())
library(dplyr)

#load formatted data
df2<-readRDS(file = here::here("outputs/df_filt_trans_one_hot.rds"))
df2[sapply(df2, is.factor)] <- lapply(df2[sapply(df2, is.factor)],
                                     as.integer)

#dissimilarity matrix
library(cluster)
gower_df2 <- daisy(df2, metric = "gower" )

summary(gower_df2)

#convert gower df to distance matrix
dataset_dist2 <- stats::as.dist(gower_df2)

library(vegan)
#run PCOA
#dataset_pcoa <- ape::pcoa(dataset_dist)
dataset_pcoa2 <- wcmdscale(d = dataset_dist2, eig = TRUE, add = "lingoes")
traitvectors12=scores(envfit(dataset_pcoa2, df2, na.rm = T, add = T),display="vectors")
traitvectors34=scores(envfit(dataset_pcoa2, df2, na.rm = T, add = T,choices=c(3,4)),display="vectors")
traitd=as.data.frame(traitvectors12)
traitd$Dim3=as.data.frame(traitvectors34)$Dim3
traitd$Dim4=as.data.frame(traitvectors34)$Dim4
traitd$trait=rownames(traitd)

speciesv=scores(dataset_pcoa2,display="species")
speciesd=as.data.frame(speciesv)

# use original traits (multistate) for visualization
dft<-readRDS(file = here::here("outputs/df_filt_trans_one_hot.rds"))
df_ord=as.data.frame(c(dft,speciesd))

# read clusters
clusters_ward_one_hot=readRDS(file = here::here("outputs/clust_num_k_2_7_ward_one_hot.rds"))
clusters_ward=readRDS(file = here::here("outputs/clust_num_k_2_7_ward.rds"))
clusters_kpro_one_hot=readRDS(file = here::here("outputs/clust_num_k_2_7_kpro_one_hot.rds"))
clusters_kpro=readRDS(file = here::here("outputs/clust_num_k_2_7_kpro.rds"))
clusters_pam_one_hot=readRDS(file = here::here("outputs/clust_num_k_2_7_pam_one_hot.rds"))
clusters_pam=readRDS(file = here::here("outputs/clust_num_k_2_7_pam.rds"))
clusters_density_one_hot=readRDS(file = here::here("outputs/clust_num_k_2_7_density_one_hot.rds"))
clusters_density=readRDS(file = here::here("outputs/clust_num_k_2_7_density.rds"))


#select number of clusters
# divisive.clust = 3
# aggl.clust.c (complete) = 4
# aggl.clust.a (average) = 2
# aggl.clust.w (ward) = 5

#clust.num[clust.num==1]="showy perennials"
#clust.num[clust.num==2]="dioecious trees"
#clust.num[clust.num==3]="small/short-lived herbs"
#clust.num[clust.num==4]="showy trees"
#clust.num[clust.num==5]="showy trees"
#clust.num[clust.num==6]="wind-pollinated sexually monomorphic"

library(ggrepel)

#minimal value for arrow length
minarrow=0.2

#number of clusters to use for plotting
clust.num=data.frame(clust.num=clusters_ward[,"6clusters"])
rownames(clust.num)=rownames(clusters_ward)
df_ord_clust=as.data.frame(c(df_ord,clust.num))

#plot points on first two axes, coloured by cluster
p1=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(clust.num)))+geom_segment(data=traitd[sqrt(traitd$Dim1^2+traitd$Dim2^2)>minarrow,],aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd[sqrt(traitd$Dim1^2+traitd$Dim2^2)>minarrow,],aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust, geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(clust.num)), alpha = 0.25)+theme_bw() + theme(legend.position="top",legend.title=element_blank()) 

#plot points on second and third axis, coloured by cluster
p2=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim2,y=Dim3,color = as.factor(clust.num)))+geom_segment(data=traitd[sqrt(traitd$Dim2^2+traitd$Dim3^2)>minarrow,],aes(x=0,y=0,xend=Dim2/2,yend=Dim3/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd[sqrt(traitd$Dim2^2+traitd$Dim3^2)>minarrow,],aes(x=Dim2/2,y=Dim3/2,label=trait))+  stat_ellipse(data=df_ord_clust, geom = "polygon", aes(x=Dim2,y=Dim3,fill = as.factor(clust.num)), alpha = 0.25)+theme_bw() + theme(legend.position="none")

#plot points on third and fourth axis, coloured by cluster
p3=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim3,y=Dim4,color = as.factor(clust.num)))+geom_segment(data=traitd[sqrt(traitd$Dim3^2+traitd$Dim4^2)>minarrow,],aes(x=0,y=0,xend=Dim3/2,yend=Dim4/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd[sqrt(traitd$Dim3^2+traitd$Dim4^2)>minarrow,],aes(x=Dim3/2,y=Dim4/2,label=trait))+  stat_ellipse(data=df_ord_clust, geom = "polygon", aes(x=Dim3,y=Dim4,fill = as.factor(clust.num)), alpha = 0.25)+theme_bw() + theme(legend.position="none")

library(patchwork)
p1 / p2 / p3

ggsave("figures/pcoa_hclust_k6_loadings_one_hot.png",
       width = 20,
       height = 40,
       units = 'cm')

clust.num=data.frame(clust.num=clusters_ward[,"3clusters"])
rownames(clust.num)=rownames(clusters_ward)
df_ord_clust=as.data.frame(c(df_ord,clust.num))

library(ggrepel)

#plot points on first two axes, coloured by cluster
p1=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(clust.num)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust, geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(clust.num)), alpha = 0.25) + theme(legend.position="none")

#plot points on second and third axis, coloured by cluster
p2=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim2,y=Dim3,color = as.factor(clust.num)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim2/2,yend=Dim3/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim2/2,y=Dim3/2,label=trait))+  stat_ellipse(data=df_ord_clust, geom = "polygon", aes(x=Dim2,y=Dim3,fill = as.factor(clust.num)), alpha = 0.25) + theme(legend.position="none")

#plot points on third and fourth axis, coloured by cluster
p3=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim3,y=Dim4,color = as.factor(clust.num)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim3/2,yend=Dim4/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim3/2,y=Dim4/2,label=trait))+  stat_ellipse(data=df_ord_clust, geom = "polygon", aes(x=Dim3,y=Dim4,fill = as.factor(clust.num)), alpha = 0.25) + theme(legend.position="none")

library(patchwork)
p1 / p2 / p3

ggsave("figures/pcoa_hclust_k3_loadings_one_hot.png",
       width = 20,
       height = 40,
       units = 'cm')


#compare clustering methods
# "normal encoding" first

numclust="6clusters"

clust.num=data.frame(clust.num=clusters_ward[,numclust])
rownames(clust.num)=rownames(clusters_ward)
df_ord_clust=as.data.frame(c(df_ord,clust.num))

p1=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(clust.num)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust, geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(clust.num)), alpha = 0.25) + theme(legend.position="none") + ggtitle("ward")

clust.num=data.frame(clust.num=clusters_kpro[,numclust])
rownames(clust.num)=rownames(clusters_kpro)
df_ord_clust=as.data.frame(c(df_ord,clust.num))

p2=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(clust.num)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust, geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(clust.num)), alpha = 0.25) + theme(legend.position="none") + ggtitle("kproto")

clust.num=data.frame(clust.num=clusters_pam[,numclust])
rownames(clust.num)=rownames(clusters_pam)
df_ord_clust=as.data.frame(c(df_ord,clust.num))

p3=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(clust.num)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust, geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(clust.num)), alpha = 0.25) + theme(legend.position="none") + ggtitle("pam")

clust.num=data.frame(clust.num=clusters_density[,numclust])
rownames(clust.num)=rownames(clusters_density)
df_ord_clust=as.data.frame(c(df_ord,clust.num))

p4=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(clust.num)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust, geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(clust.num)), alpha = 0.25) + theme(legend.position="none") + ggtitle("density")

p1 + p2 + p3 + p4

ggsave("figures/pcoa_clustering_k6_loadings.png",
       width = 40,
       height = 20,
       units = 'cm')

#one-hot based clustering

numclust="6clusters"

clust.num=data.frame(clust.num=clusters_ward_one_hot[,numclust])
rownames(clust.num)=rownames(clusters_ward_one_hot)
df_ord_clust=as.data.frame(c(df_ord,clust.num))

p1=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(clust.num)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust, geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(clust.num)), alpha = 0.25) + theme(legend.position="none") + ggtitle("ward")

clust.num=data.frame(clust.num=clusters_kpro_one_hot[,numclust])
rownames(clust.num)=rownames(clusters_kpro_one_hot)
df_ord_clust=as.data.frame(c(df_ord,clust.num))

p2=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(clust.num)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust, geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(clust.num)), alpha = 0.25) + theme(legend.position="none") + ggtitle("kproto")

clust.num=data.frame(clust.num=clusters_pam_one_hot[,numclust])
rownames(clust.num)=rownames(clusters_pam_one_hot)
df_ord_clust=as.data.frame(c(df_ord,clust.num))

p3=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(clust.num)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust, geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(clust.num)), alpha = 0.25) + theme(legend.position="none") + ggtitle("pam")

clust.num=data.frame(clust.num=clusters_density_one_hot[,numclust])
rownames(clust.num)=rownames(clusters_density_one_hot)
df_ord_clust=as.data.frame(c(df_ord,clust.num))

p4=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(clust.num)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust, geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(clust.num)), alpha = 0.25) + theme(legend.position="none") + ggtitle("density")

p1 + p2 + p3 + p4

ggsave("figures/pcoa_clustering_one_hot_k6_loadings.png",
       width = 40,
       height = 20,
       units = 'cm')

# read robust groups 
groups_kpro_one_hot=readRDS(file = here::here("outputs/robust_vect_kpro_full_one_hot.rds"))
groups_kpro=readRDS(file = here::here("outputs/robust_vect_kpro_full.rds"))
groups_pam_one_hot=readRDS(file = here::here("outputs/robust_vect_pam_full_one_hot.rds"))
groups_pam=readRDS(file = here::here("outputs/robust_vect_pam_full.rds"))
groups_density_one_hot=readRDS(file = here::here("outputs/robust_vect_dens_full_one_hot.rds"))
groups_density=readRDS(file = here::here("outputs/robust_vect_dens_full.rds"))


rownames(df_ord)=rownames(groups_kpro)

df_ord_clust=as.data.frame(c(df_ord,as.data.frame(groups_kpro)))
names(df_ord_clust)[401]="group"

p1=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(group)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust[!is.na(df_ord_clust$group),], geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(group)), alpha = 0.25) + theme(legend.position="none") + ggtitle("kproto")

df_ord_clust=as.data.frame(c(df_ord,as.data.frame(groups_pam)))
names(df_ord_clust)[401]="group"

p2=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(group)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust[!is.na(df_ord_clust$group),], geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(group)), alpha = 0.25) + theme(legend.position="none") + ggtitle("pam")

df_ord_clust=as.data.frame(c(df_ord,as.data.frame(groups_density)))
names(df_ord_clust)[401]="group"

p3=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(group)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust[!is.na(df_ord_clust$group),], geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(group)), alpha = 0.25) + theme(legend.position="none") + ggtitle("density")

p1 + p2 + p3

ggsave("figures/pcoa_clustering_robust_loadings.png",
       width = 40,
       height = 20,
       units = 'cm')

#one-hot

df_ord_clust=as.data.frame(c(df_ord,as.data.frame(groups_kpro_one_hot)))
names(df_ord_clust)[401]="group"

p1=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(group)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust[!is.na(df_ord_clust$group),], geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(group)), alpha = 0.25) + theme(legend.position="none") + ggtitle("kproto")

df_ord_clust=as.data.frame(c(df_ord,as.data.frame(groups_pam_one_hot)))
names(df_ord_clust)[401]="group"

p2=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(group)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust[!is.na(df_ord_clust$group),], geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(group)), alpha = 0.25) + theme(legend.position="none") + ggtitle("pam")

df_ord_clust=as.data.frame(c(df_ord,as.data.frame(groups_density_one_hot)))
names(df_ord_clust)[401]="group"

p3=ggplot() + geom_point(data=df_ord_clust,aes(x=Dim1,y=Dim2,color = as.factor(group)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))+  stat_ellipse(data=df_ord_clust[!is.na(df_ord_clust$group),], geom = "polygon", aes(x=Dim1,y=Dim2,fill = as.factor(group)), alpha = 0.25) + theme(legend.position="none") + ggtitle("density")

p1 + p2 + p3

ggsave("figures/pcoa_clustering_robust_one_hot_loadings.png",
       width = 40,
       height = 20,
       units = 'cm')

#3d

#clust.num <- cutree(aggl.clust.w, k = 6)

#library(rgl)
#open3d()
#get_colors <- function(groups, group.col = palette()){
#  groups <- as.factor(groups)
#  ngrps <- length(levels(groups))
#  if(ngrps > length(group.col)) 
#    group.col <- rep(group.col, ngrps)
#  color <- group.col[as.numeric(groups)]
#  names(color) <- as.vector(groups)
#  return(color)
#}

#plot_loadings3d <- function(loadings,scale=0.5){
#for(i in 1:nrow(loadings)){
#arrow3d(p0=c(0,0,0),p1=c(loadings[i,"Dim1"]*scale,loadings[i,"Dim2"]*scale,loadings[i,"Dim3"]*scale),type="lines",barblen=0.01)
#text3d(x=loadings[i,"Dim1"]*scale,y=loadings[i,"Dim2"]*scale,z=loadings[i,"Dim3"]*scale,loadings$trait[i],cex=0.5)
#}
#}

#rgl.spheres(df_ord[,"Dim1"], df_ord[,"Dim2"], df_ord[,"Dim3"], r = 0.005,color = get_colors(clust.num)) 
#plot_loadings3d(traitd)
#clear3d()

#end 3d


