#Same as 8_pcoa.R but removes reproductive traits
rm(list=ls())
library(dplyr)
library(ggplot2)

#load formatted data
df_orig<-readRDS(file = here::here("outputs/df_filt_trans.rds"))

#remove reproductive traits
df<-subset(df_orig,select=-c(SexualSystem,Mating,FlowerSex))

#dissimilarity matrix calculation
library(cluster)
gower_df <- daisy(df,
                  metric = "gower" )

summary(gower_df)

dataset_dist <- stats::as.dist(gower_df)
dataset_pcoa <- ape::pcoa(dataset_dist)

#plot PCOA points on first two axes
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.75,
    size=2,
    stroke = 0.5
  ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))

ggsave("figures/scatter_pcoa.png",
       width = 15,
       height = 15,
       units = 'cm')

#relative eigenvalues
eig_df<-data.frame(c(1:9),dataset_pcoa$values$Relative_eig[1:9])
colnames(eig_df)<-c("axis","relative_eigenvalue")
eig_df$axis<-as.character(eig_df$axis)

ggplot(eig_df, aes(x=axis, y=relative_eigenvalue)) + 
  geom_bar(stat = "identity")
ggsave("figures/rel_eig_pcoa_full.png")

#plot points on first two axes coloured by different traits
p1 <- ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2, fill = as.factor(df_orig$SexualSystem))) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.5,
    size=4,
    stroke = 0.5
  ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))

p2 <- ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2, fill = as.factor(df_orig$Mating))) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.5,
    size=4,
    stroke = 0.5
  ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))

p3 <- ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2, fill = as.factor(df_orig$FlowerSex))) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.5,
    size=4,
    stroke = 0.5
  ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))


###
# Combined
###
library(patchwork)
p1 / p2 / p3
ggsave("figures/scatter_pcoa_no_repro.png",
       width = 20,
       height = 40,
       units = 'cm')

#plot points on first two axes with two variables: reproductive system  (color) and woodiness (shape)
p1 <- ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2, fill = as.factor(df_orig$SexualSystem), shape=as.factor(df_orig$Woodiness))) +
  geom_point(
    color="black",
    #    shape=21,
    alpha=0.5,
    size=4,
    stroke = 0.5
  ) + 
  scale_shape_manual(values=c(21,22,23)) +
  guides(fill = guide_legend(override.aes = list(shape = 24) ) ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))

p2 <- ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2, fill = as.factor(df_orig$Mating), shape=as.factor(df_orig$Woodiness))) +
  geom_point(
    color="black",
    #    shape=21,
    alpha=0.5,
    size=4,
    stroke = 0.5
  ) + 
  scale_shape_manual(values=c(21,22,23)) +
  guides(fill = guide_legend(override.aes = list(shape = 24) ) ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))

p3 <- ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2, fill = as.factor(df_orig$FlowerSex), shape=as.factor(df_orig$Woodiness))) +
  geom_point(
    color="black",
    #    shape=21,
    alpha=0.5,
    size=4,
    stroke = 0.5
  ) + 
  scale_shape_manual(values=c(21,22,23)) +
  guides(fill = guide_legend(override.aes = list(shape = 24) ) ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))

#combine
p1 / p2 / p3

ggsave("figures/scatter_pcoa_no_repro_wood.png",
       width = 20,
       height = 40,
       units = 'cm')

