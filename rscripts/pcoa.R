rm(ls=list())
library(dplyr)
library(ggplot2)

#load formatted data
df<-readRDS(file = here::here("outputs/df_filt.rds"))

#dissimilarity matrix calculation
library(cluster)
gower_df <- daisy(df,
                  metric = "gower" )

summary(gower_df)

dataset_dist <- stats::as.dist(gower_df)
dataset_pcoa <- ape::pcoa(dataset_dist)

#relative eigenvalues
eig_df<-data.frame(c(1:9),dataset_pcoa$values$Relative_eig[1:9])
colnames(eig_df)<-c("axis","relative_eigenvalue")
eig_df$axis<-as.character(eig_df$axis)

ggplot(eig_df, aes(x=axis, y=relative_eigenvalue)) + 
  geom_bar(stat = "identity")
ggsave("figures/rel_eig_pcoa_full.pdf")

#plot points on first two axes
p1 <- ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2, fill = as.factor(df$SexualSystem))) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.5,
    size=4,
    stroke = 0.5
  ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))

p2 <- ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2, fill = as.factor(df$Mating))) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.5,
    size=4,
    stroke = 0.5
  ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))

p3 <- ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2, fill = as.factor(df$FlowerSex))) +
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
ggsave("figures/scatter_pcoa_full.pdf",
       width = 20,
       height = 40,
       units = 'cm')

#plot points on first two axes with two variables: reproductive system and woodiness

p1 <- ggplot(data.frame(dataset_pcoa$vecto
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
  xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))

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
  xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))

p3 <- ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2, fill = as.factor(df$FlowerSex), shape=as.factor(df$Woodiness))) +
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

p1 / p2 / p3
ggsave("figures/scatter_pcoa_full_wood.pdf",
       width = 20,
       height = 40,
       units = 'cm')


