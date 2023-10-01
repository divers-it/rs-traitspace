rm(list=ls())

#load libraries
library(dplyr)
library(ggplot2)
library(cluster)
library(patchwork)
library(vegan)
library(ggrepel)

#load original data set
df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

#dissimilarity matrix calculation
gower_df <- daisy(df,
                  metric = "gower" )
summary(gower_df)

#load one-hot data set
df2<-readRDS(file = here::here("outputs/one_hot_6_df_filt_trans.rds"))

#NOTE: not sure this is necessary
#change factors to integers
df2[sapply(df2, is.factor)] <- lapply(df2[sapply(df2, is.factor)],as.integer)

#dissimilarity matrix calculation
gower_df2 <- daisy(df2,
                  metric = "gower" )
summary(gower_df2)

#check names
setdiff(labels(gower_df),labels(gower_df2))
labels(gower_df)==labels(gower_df2)

#linear model
mod<-summary(lm(gower_df~gower_df2))
mod

#compare pairwaise distances of matrices with missing data and with imputed
png("figures/one_hot_8_scatterplot_dist_og_vs_one_hot.png",width = 500,height = 500)
plot(gower_df,gower_df2,xlim=c(0,1),ylim=c(0,1),xlab="Original",ylab="One-hot")
abline(0,1,lty=2,col="red",lwd=2)
text(x=0.15, y=0.9, labels=paste("R-squared =",round(mod$r.squared,3)))
dev.off()


#make into distance matrix
dataset_dist2 <- stats::as.dist(gower_df2)

#NOT RUN:
#original ape pcoa
dataset_pcoa2b <- ape::pcoa(dataset_dist2)
dataset_pcoa2b$values$Eigenvalues

#Recalculate relative eigenvalues by removing negative eigenvalues as in Mouillot et al.  
ev_pcoa <- dataset_pcoa2b$values$Eigenvalues
ev_pcoa_g0 <- ev_pcoa[ev_pcoa>0]
rel_ev_pcoa_g0 <- ev_pcoa_g0/sum(ev_pcoa_g0)

#Make data frame of first 9 relative eigenvalues
eig_df<-data.frame(c(1:9),rel_ev_pcoa_g0[1:9])
colnames(eig_df)<-c("pcoa_axis","relative_eigenvalue")
eig_df$pcoa_axis<-as.character(eig_df$pcoa_axis)

#prop variation first two axes explain
sum(rel_ev_pcoa_g0[1:2])

#plot barplot of first 9 relative eigenvalues
ggplot(eig_df, aes(x=pcoa_axis, y=relative_eigenvalue)) + 
  geom_bar(stat = "identity")
ggsave("figures/barplot_relative_eigenvalues_pcoa_one_hot.png")

#Change pcoa method (using the method from the "vegan" package)
#Weight version
dataset_pcoa2 <- wcmdscale(d = dataset_dist2, eig = TRUE)
dataset_pcoa2

#difference between ape and vegan pcoa eigenvalues
dataset_pcoa2$eig[1:20]-dataset_pcoa2b$values$Eigenvalues[1:20]

#extract locations of species on PCoA axes 
species_scores <- as.data.frame(scores(dataset_pcoa2,display="species"))

#check names
rownames(species_scores)==rownames(df2)

#missing data per row (species) for plot
missDat<-rowSums(apply(is.na(df2),2,as.numeric))/ncol(df2)

#plot PCOA points on first two axes coloured by missing data
ggplot(species_scores, aes(x = Dim1, y = Dim2)) +
  geom_point(
    aes(color=missDat),
    shape=16,
    alpha=0.75,
    size=2,
    stroke = 0.5
  )  +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

#save plot
ggsave("figures/scatter_pcoa_missing_one_hot.png",
       width = 20,
       height = 15,
       units = 'cm')

#check matching names of df an df2
rownames(df)==rownames(df2)

#NOTE: unsure why one row is being removed for missing values
#plot points on first two axes with point style changed by two variables
#reproductive systems (color) and woodiness (shape)
p1 <- ggplot(species_scores, aes(x = Dim1, y = Dim2, fill = as.factor(df$SexualSystem), shape=as.factor(df$Woodiness))) +
  geom_point(
    color="black",
    #    shape=21,
    alpha=0.5,
    size=4,
    stroke = 0.5
  ) + 
  scale_shape_manual(values=c(21,22,23)) +
  guides(fill = guide_legend(override.aes = list(shape = 24))) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2))) + 
  theme(legend.position="bottom",
        legend.title = element_blank(),
        legend.box="vertical", 
        legend.margin=margin())

# mating system (colour) and woodiness (shape)
p2 <- ggplot(species_scores, aes(x = Dim1, y = Dim2, fill = as.factor(df$Mating), shape=as.factor(df$Woodiness))) +
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
ggsave("figures/scatter_pcoa_coloured_by_traits_one_hot.png",
       width = 30,
       height = 15,
       units = 'cm')


#make data frame with scores of environemntal vectors
traitd <- as.data.frame(scores(envfit(dataset_pcoa2, df2, na.rm = T, add = T,choices=c(1,2,3,4)),display="vectors"))
traitd$trait <- rownames(traitd)

#plot PCoA scatterplot with "loadings" arrows
p1 <- ggplot() + geom_point(data=species_scores,aes(x=Dim1,y=Dim2))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))

p2 <- ggplot() + geom_point(data=species_scores,aes(x=Dim2,y=Dim3))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim2/2,yend=Dim3/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim2/2,y=Dim3/2,label=trait))

p3 <- ggplot() + geom_point(data=species_scores,aes(x=Dim3,y=Dim4))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim3/2,yend=Dim4/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim3/2,y=Dim4/2,label=trait))

#combine
p1 / p2 / p3

ggsave("figures/pcoa_loadings_one_hot.png",
       width = 20,
       height = 40,
       units = 'cm')

#NOT RUN: combine two plots above 

#combine trait data and scores
#df_ord <- cbind(df,species_scores)

#ggplot() + geom_point(data=df_ord,aes(x=Dim1,y=Dim2,shape=as.factor(Woodiness),color=as.factor(Mating)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))

#ggsave("figures/pcoa_woody_mating_one_hot.pdf")

#ggplot() + geom_point(data=df_ord,aes(x=Dim1,y=Dim2,shape=as.factor(Woodiness),color=as.factor(FlowerSex)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))

#ggsave("figures/pcoa_woody_flowersex_one_hot.pdf")

#ggplot() + geom_point(data=df_ord,aes(x=Dim1,y=Dim2,shape=as.factor(Woodiness),color=as.factor(SexualSystem)))+geom_segment(data=traitd,aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2),arrow=arrow(),col="blue")+geom_text_repel(data=traitd,aes(x=Dim1/2,y=Dim2/2,label=trait))

#ggsave("figures/pcoa_woody_sexualsystem_one_hot.pdf")



