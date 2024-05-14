rm(list=ls())

#load libraries
library(dplyr)
library(ggplot2)
library(cluster)
library(patchwork)
library(wesanderson)
library(gawdis)

#load data set
df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

#dissimilarity matrix calculation
gower_df <- daisy(df,
                  metric = "gower" )
summary(gower_df)

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

#plot PCOA points on first two axes coloured by missing data
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    #aes(color=missDat),
    shape=16,
    alpha=0.75,
    size=2,
    stroke = 0.5
  ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))


#Make data frame of first 9 relative eigenvalues
eig_df<-data.frame(c(1:9),rel_ev_pcoa_g0[1:9])
colnames(eig_df)<-c("pcoa_axis","relative_eigenvalue")
eig_df$pcoa_axis<-as.character(eig_df$pcoa_axis)

#prop variation first two axes explain
sum(rel_ev_pcoa_g0[1:2])

#plot barplot of first 9 relative eigenvalues
ggplot(eig_df, aes(x=pcoa_axis, y=relative_eigenvalue)) + 
  geom_bar(stat = "identity")

####
# ---- Gawdis ----
####

# https://cran.r-project.org/web/packages/gawdis/vignettes/gawdis.html

ex1 <- gowdis(df)
View(as.matrix(round(ex1, 3))) ##just to see only 3 decimals, enough

straightgowdis.2<-gawdis(df, w.type = "equal", silent = T)#we compute 'normal' gower with the new function because it provides more results
cors.gow<-attr(straightgowdis.2,"correls")
cors.gow

#must use because of NA (or equal weights)
iterations <- gawdis(df, w.type ="optimized")
attr(iterations, "correls")
attr(iterations, "weights")

#make into distance object
dataset_dist <- stats::as.dist(iterations)

#run PCoA on distance matrix
dataset_pcoa <- ape::pcoa(dataset_dist)

#Recalculate relative eigenvalues by removing negative eigenvalues as in Mouillot et al.  
ev_pcoa <- dataset_pcoa$values$Eigenvalues
ev_pcoa_g0 <- ev_pcoa[ev_pcoa>0]
rel_ev_pcoa_g0 <- ev_pcoa_g0/sum(ev_pcoa_g0)

#check names
rownames(dataset_pcoa$vectors)==rownames(df)

#plot PCOA points on first two axes coloured by missing data
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    aes(color=df$FlowerSex,shape=df$Woodiness),
    alpha=0.75,
    size=2,
    stroke = 0.5
  ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))


#Make data frame of first 9 relative eigenvalues
eig_df<-data.frame(c(1:9),rel_ev_pcoa_g0[1:9])
colnames(eig_df)<-c("pcoa_axis","relative_eigenvalue")
eig_df$pcoa_axis<-as.character(eig_df$pcoa_axis)

#prop variation first two axes explain
sum(rel_ev_pcoa_g0[1:2])

#plot barplot of first 9 relative eigenvalues
ggplot(eig_df, aes(x=pcoa_axis, y=relative_eigenvalue)) + 
  geom_bar(stat = "identity")

