rm(list=ls())

# load libraries
library(dplyr)
library(ggplot2)
library(cluster)
library(patchwork)
library(wesanderson)
library(ape)
library(vegan)
library(rphylopic)

# load data set
df <- readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

# remove some dioecious species
dio_sp <- rownames(df[df$SexualSystem == "dimorphic" & df$FlowerSex == "unisexual",])
df <- df[!rownames(df)%in%sample(dio_sp, 30),]


# dissimilarity matrix calculation
gower_df <- daisy(df,
                  metric = "gower" )
summary(gower_df)

## Run and plot PCOA ----

# make into distance object
dataset_dist <- stats::as.dist(gower_df)

# run PCoA on distance matrix
dataset_pcoa <- ape::pcoa(dataset_dist)

# Recalculate relative eigenvalues by removing negative eigenvalues as in Mouillot et al.  
ev_pcoa <- dataset_pcoa$values$Eigenvalues
ev_pcoa_g0 <- ev_pcoa[ev_pcoa>0]
rel_ev_pcoa_g0 <- ev_pcoa_g0/sum(ev_pcoa_g0)

# check names
table(rownames(dataset_pcoa$vectors)==rownames(df))


# plot PCOA points on first two axes coloured by missing data
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    aes(color=df$FlowerSex, shape=df$Woodiness),
    alpha=0.75,
    size=3,
    stroke = 0.5
  ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))
