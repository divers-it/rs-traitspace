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

## Make traits ordinal ----

# increasing outcrossing rate
df$Mating <- factor(df$Mating, order = TRUE, levels = c("selfing", "mixed", "outcrossing"))

# increasing woodiness
# df$Woodiness <- factor(df$Woodiness, order = TRUE, levels = c("herbaceous","herbaceous_woody", "woody"))

# increasing sexual segregation
# df$SexualSystem <- factor(df$SexualSystem, order = TRUE, levels = c("monomorphic","dimorphic_monomorphic", "dimorphic"))

# increasing lifespan
df$Lifespan <- factor(df$Lifespan, order = TRUE, levels = c("short","long_short","long"))

# increasingly biotic pollination
# df$Pollination <- factor(df$Pollination, order = TRUE, levels = c("abiotic", "abiotic_autonomous","abiotic_biotic","autonomous", "autonomous_biotic", "biotic"))

# increasingly biotic dispersal
# df$DispersalMode <- factor(df$DispersalMode, order = TRUE, levels = c("abiotic", "abiotic_autonomous", "abiotic_autonomous_biotic", "abiotic_biotic","autonomous", "autonomous_biotic", "biotic"))

# increasing dispersal distance
df$DispersalDist <- factor(df$DispersalDist, order = TRUE, levels = c("short","long_short","long"))

# increasing sexual segregation
# df$FlowerSex <- factor(df$DispersalDist, order = TRUE, levels = c("bisexual", "bisexual_unisexual", "unisexual"))

# increasing ovary height
df$OvaryPosition <- factor(df$OvaryPosition, order = TRUE, levels = c("inferior", "inferior_superior", "intermediate", "intermediate_superior", "superior"))

# increasingly bright flowers
df$Showiness <- factor(df$Showiness, order = TRUE, levels = c("dull", "bright_dull", "bright"))

# NOT MADE ORDINAL
# df$Climbing 
# df$Aquatic
# df$FloralReward

# check structure
str(df)
df_ord <- df

# dissimilarity matrix calculation
gower_ord <- daisy(df_ord,
                  metric = "gower" )
summary(gower_ord)

## Compare original and ordinal data sets ----

# load data set
df <- readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

# dissimilarity matrix calculation
gower_df <- daisy(df, metric = "gower" )
summary(gower_df)

# check names
table(labels(gower_ord)==labels(gower_df))

# check correlation
mant <- vegan::mantel(as.matrix(gower_df),as.matrix(gower_ord))
mant

# compare pairwise distances of matrices with missing data and with imputed
# png("figures/8.0_scatterplot_dist_missing_vs_imputed.png",width = 500,height = 500)
plot(gower_df,gower_ord,xlim=c(0,1),ylim=c(0,1),xlab="Original",ylab="Ordinal")
abline(0,1,lty=2,col="red",lwd=2)
text(x=0.15, y=0.9, labels=paste("Mantel statistic r =",round(mant$statistic,3)))
# dev.off()

## Run and plot PCOA ----

# make into distance object
dataset_dist <- stats::as.dist(gower_ord)

# run PCoA on distance matrix
dataset_pcoa <- ape::pcoa(dataset_dist)

# Recalculate relative eigenvalues by removing negative eigenvalues as in Mouillot et al.  
ev_pcoa <- dataset_pcoa$values$Eigenvalues
ev_pcoa_g0 <- ev_pcoa[ev_pcoa>0]
rel_ev_pcoa_g0 <- ev_pcoa_g0/sum(ev_pcoa_g0)

# check names
table(rownames(dataset_pcoa$vectors)==rownames(df))

# read original reproductive systems
rs <- read.csv("outputs/original_reproductive_systems.csv")

table(rs$species == rownames(df))

# plot PCOA points on first two axes coloured by missing data
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    aes(color=rs$RS,
        shape=df$Woodiness),
    alpha=0.75,
    size=2,
    stroke = 0.5
  ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

# Make data frame of first 9 relative eigenvalues
eig_df<-data.frame(c(1:9),rel_ev_pcoa_g0[1:9])
colnames(eig_df)<-c("pcoa_axis","relative_eigenvalue")
eig_df$pcoa_axis<-as.character(eig_df$pcoa_axis)

# Proportion of variation first two axes explain
sum(rel_ev_pcoa_g0[1:2])

# plot barplot of first 9 relative eigenvalues
ggplot(eig_df, aes(x=pcoa_axis, y=relative_eigenvalue)) + 
  geom_bar(stat = "identity")
# ggsave("figures/8.0_barplot_relative_eigenvalues_pcoa.png")

# reorder is close to order, but is made to change the order of the factor levels.
rs$RS <- factor(rs$RS, levels = c("Dioecy", "Monoecy", "Monocliny: selfing", "Monocliny: mixed", "Monocliny: outcrossing", "unknown"))

# palette
mypal <- c("#009E73","#CC79A7", colorRampPalette(c( "#E69F00","#D55E00"))(3), "black")

# scatterplot with reprdocutive systems
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2))+
  ggforce::geom_mark_hull(concavity = 5,expand=0,radius=0, aes(colour=as.factor(rs$RS), fill=as.factor(rs$RS),filter = rs$RS != 'unknown')) +
  ggforce::geom_mark_hull(concavity = 5,expand=0,radius=0, colour = "darkgrey") +
  scale_fill_manual(values=mypal) +
  scale_colour_manual(values=mypal, name = "Reproductive system") +
  geom_point(aes(colour=as.factor(rs$RS), shape=as.factor(rs$RS))) +
  scale_shape_manual(values=c(16,16,16,16,16,21)) +
  guides(fill ="none", shape = "none", color = guide_legend(override.aes = list(shape = c(16,16,16,16,16,21) ) ) ) +
  xlab(paste("PCoA Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("PCoA Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2))) +
  theme_bw() + theme(
    panel.border = element_blank(),
    #panel.grid.major = element_line(colour = "darkgrey"),
    #panel.grid.minor = element_line(colour = "grey"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.85, 0.85),
    legend.text = element_text(size=12),
    legend.title = element_text(size=14),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size=12),
    axis.title = element_text(size=14))



