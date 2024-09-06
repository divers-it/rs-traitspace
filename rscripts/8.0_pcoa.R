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

# dissimilarity matrix calculation
gower_df <- daisy(df,
                  metric = "gower" )
summary(gower_df)

## Compare original and imputed data sets ----

# read in imputed data set
df_no_miss <-read.csv("outputs/imputed_with_phylo.csv",row.names=1,stringsAsFactors = TRUE)

# dissimilarity matrix calculation
gower_df_no_miss <- daisy(df_no_miss,
                          metric = "gower" )
summary(gower_df_no_miss)

# check names
table(labels(gower_df)==gsub("_"," ",labels(gower_df_no_miss)))

# check correlation
mant <- vegan::mantel(as.matrix(gower_df),as.matrix(gower_df_no_miss))
mant

# compare pairwise distances of matrices with missing data and with imputed
# png("figures/8.0_scatterplot_dist_missing_vs_imputed.png",width = 500,height = 500)
plot(gower_df,gower_df_no_miss,xlim=c(0,1),ylim=c(0,1),xlab="Original",ylab="Imputed")
abline(0,1,lty=2,col="red",lwd=2)
text(x=0.15, y=0.9, labels=paste("Mantel statistic r =",round(mant$statistic,3)))
# dev.off()

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

# missing data per row (species) for plot
missDat <- rowSums(apply(is.na(df),2,as.numeric))/ncol(df)

# plot PCOA points on first two axes coloured by missing data
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

# save plot
# ggsave("figures/8.0_scatterplot_pcoa_missing.png",
#        width = 20,
#        height = 15,
#        units = 'cm')

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

####
##  Phylogenetic distance vs. trait distance ----
####

# read in phylogenetic tree
phy<-read.tree("outputs/pruned_tree.tre")

# add underscores to names
df2 <- df
rownames(df2) <- gsub(" ", "_", x = rownames(df2))

# in dataset but not in phylo
setdiff(rownames(df2),phy$tip.label)

# in phylo but not in dataset
setdiff(phy$tip.label,rownames(df2))

# get pairwise phylogenetic distances
phy_dm <- cophenetic.phylo(phy)

#dissimilarity matrix calculation
gower_df2 <- daisy(df2, metric = "gower" )
gower_df2 <- as.matrix(gower_df2)

#order rows and columns and check
phy_dm <- phy_dm[order(rownames(phy_dm)), order(colnames(phy_dm))]
gower_df2 <- gower_df2[order(rownames(gower_df2)), order(colnames(gower_df2))]
table(rownames(gower_df2)==rownames(phy_dm))
table(colnames(gower_df2)==colnames(phy_dm))

#compare distance matrices
plot(gower_df2,phy_dm,xlab="Trait distance",ylab="Phylogenetic distance")
vegan::mantel(gower_df2,phy_dm,method = "spear")

####
## Get phylopics ----
####

#choose from the different ones for a grass species
img <- pick_phylopic(name = "Oryza sativa", n = 5)
# Get a single image uuid
uuid <- get_uuid(name = "Oryza sativa", n = 1)
# Get the image for that uuid
oryza_pp <- get_phylopic(uuid = uuid)

uuid <- get_uuid(name = "Zea mays", n = 5)
zea_pp <- get_phylopic(uuid = uuid[4])

uuid <- get_uuid(name = "Phoenix dactylifera", n = 1)
phoenix_pp <- get_phylopic(uuid = uuid)

uuid <- get_uuid(name = "Solanum dulcamara", n = 2)
solanum_pp <- get_phylopic(uuid = uuid[2])

uuid <- get_uuid(name = "Plantago lanceolata", n = 2)
plantago_pp <- get_phylopic(uuid = uuid)

uuid <- get_uuid(name = "Trithuria submersa", n = 2)
trithuria_pp <- get_phylopic(uuid = uuid)

uuid <- get_uuid(name = "Commelina communis", n = 2)
commelina_pp <- get_phylopic(uuid = uuid[1])

uuid <- get_uuid(name = "Coffea arabica", n = 3)
coffea_pp <- get_phylopic(uuid = uuid[2])

uuid <- get_uuid(name = "Cornus florida", n = 1)
cornus_pp <- get_phylopic(uuid = uuid)

####
## Figure 2: PCoA with traits and density ----
####

# Edit factors for plotting to control colors
# add new factor level "none"
df$Woodiness = factor(df$Woodiness, levels=c(levels(df$Woodiness), "None"))
#df$SexualSystem = factor(df$SexualSystem, levels=c(levels(df$SexualSystem), "None"))

# convert all NA's to None
df$Woodiness[is.na(df$Woodiness)] = "None"
#df$SexualSystem = factor(df$SexualSystem, levels=c(levels(df$SexualSystem), "None"))

# y-axis scaling
# yax <- 11/19

# PCoA scatterplot with density polygons
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", col=NA , n=100, bins=20) +
  scale_fill_distiller(palette = "Greys", direction = 1, guide = "none") +
  geom_point(
    aes(
      color = as.factor(df$FlowerSex),
      shape = as.factor(df$Woodiness)),
    # shape=21,
    alpha = 0.7,
    size = 2.5,
    stroke = 0.5) + 
  scale_color_manual(values=wes_palette("FantasticFox1", 3),
                     labels=c('Bisexual', 'Bisexual & Unisexual', 'Unisexual','No Data')) +
  scale_shape_discrete(labels=c('Herbaceous', 'Herbaceous & Woody', 'Woody','No Data')) +
  xlab(paste("PCoA Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("PCoA Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2))) +
  xlim(-0.6,0.6) + 
  ylim(-0.44,0.6) + 
  theme_bw() + theme(
    panel.border = element_blank(),
    #panel.grid.major = element_line(colour = "darkgrey"),
    #panel.grid.minor = element_line(colour = "grey"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.85, 0.8),
    legend.text = element_text(size=12),
    legend.title = element_text(size=14),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size=12),
    axis.title = element_text(size=14)) + 
  labs(
    colour = "Flower sex",
    shape = "Woodiness"
  ) +   
  add_phylopic(img=solanum_pp,x = 0.28, y=-0.225, ysize = 0.15,col = "grey30") +
  add_phylopic(img=commelina_pp,x = 0.44, y=0, ysize = 0.15,col = "grey30") +
  add_phylopic(img=phoenix_pp,x = -0.52, y=0, ysize = 0.15,col = "grey30") +
  add_phylopic(img=coffea_pp,x = -0.07, y=-0.35, ysize = 0.125,col = "grey30") +
  add_phylopic(img=trithuria_pp,x = 0.07, y=0.55, ysize = 0.125,col = "grey30") +
  add_phylopic(img=cornus_pp,x = -0.35, y=-0.26, ysize = 0.125,col = "grey30") +
  add_phylopic(img=zea_pp,x = -0.39, y=0.38, ysize = 0.125,col = "grey30") +
  annotate("text", x=-0.07, y=-0.43, label= "paste(italic(Coffea), ' ',italic(arabica))",
               col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=-0.52, y=-0.09, label= "paste(italic(Phoenix), ' ',italic(dactylifera))",
             col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=0.28, y=-0.32, label= "paste(italic(Solanum), ' ',italic(dulcamara))",
             col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=0.44, y=-0.1, label= "paste(italic(Commelina), ' ',italic(communis))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=0.07, y=0.47, label= "paste(italic(Trithuria), ' ',italic(submersa))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=-0.35, y=-0.34, label= "paste(italic(Cornus), ' ',italic(florida))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=-0.39, y=0.3, label= "paste(italic(Zea), ' ',italic(mays))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("segment", #phoenix
           #linetype=2,
           linewidth=0.75,
           x=-0.46,
           xend=-0.3305978702, 
           y=0.03,
           yend=4.934202e-02, 
           color = "black",
           alpha=0.5) + 
  annotate("segment", #solanum
           #linetype=2,
           linewidth=0.75,
           x=0.25,
           xend=0.0966219059, 
           y=-0.2,
           yend=-3.542270e-02, 
           color = "black",
           alpha=0.5) +
  annotate("segment", #coffea
           #linetype=2,
           linewidth=0.75,
           x=-0.045,
           xend=-0.0139421126, 
           y=-0.28,
           yend=-0.1305380241, 
           color = "black",
           alpha=0.5) + 
  annotate("segment", #trithuria
           #linetype=2,
           linewidth=0.75,
           x=0.07,
           xend=0.1230980263,
           y=0.45,
           yend=3.744929e-01, 
           color = "black",
           alpha=0.5) + 
  annotate("segment", #cornus
           #linetype=2,
           linewidth=0.75,
           x=-0.3,
           xend=-0.1152021084,
           y=-0.21,
           yend=-0.1413927222, 
           color = "black",
           alpha=0.5) +
  annotate("segment", #zea
           #linetype=2,
           linewidth=0.75,
           x=-0.35,
           xend=-0.1219254049,
           y=0.35,
           yend=2.600799e-01, 
           color = "black",
           alpha=0.5) +
  annotate("segment", #Commelina
           #linetype=2,
           linewidth=0.75,
           x=0.38,
           xend=0.2257547576, 
           y=0,
           yend=0.0947481926, 
           color = "black",
           alpha=0.5)
  
ggsave("figures/figure_2_pcoa.png",
       width = 35,
       #height = 35*yax,
       height = 35,
       units = 'cm')


####
## Alternate Figure 2: Plot reproductive systems (back-engineer them first) ----
####

# load package
library(ggnewscale)

# make copy of df
df_ors <- df

# recode flower sex, sexual system and mating system
df_ors$RS=rep("unknown",length(df_ors[,1]))
df_ors[df_ors$SexualSystem %in% "dimorphic" & df_ors$FlowerSex %in% "unisexual",]$RS="dioecy"
df_ors[df_ors$SexualSystem %in% "monomorphic" & df_ors$FlowerSex %in% "unisexual",]$RS="monoecy"
df_ors[df_ors$SexualSystem %in% "monomorphic" & df_ors$FlowerSex %in% "bisexual" & df_ors$Mating %in% "selfing",]$RS="monocliny: selfing"
df_ors[df_ors$SexualSystem %in% "monomorphic" & df_ors$FlowerSex %in% "bisexual" & df_ors$Mating %in% "mixed",]$RS="monocliny: mixed"
df_ors[df_ors$SexualSystem %in% "monomorphic" & df_ors$FlowerSex %in% "bisexual" & df_ors$Mating %in% "outcrossing",]$RS="monocliny: outcrossing"

# reorder is close to order, but is made to change the order of the factor levels.
df_ors$RS <- factor(df_ors$RS, levels = c("dioecy", "monoecy", "monocliny: selfing", "monocliny: mixed", "monocliny: outcrossing", "unknown"))

# color palette
# one colour for dioecy, one for monoecy, then three on a gradient for bisexual selfing->outcrossing, and one for NA

mypal <- c("darkred","darkorange",
           colorRampPalette(c( "lightblue","navyblue"))(3),
           "grey")

# PCoA scatterplot with density polygons
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", col=NA , n=100, bins=20) +
  scale_fill_distiller(palette = "Greys", direction = 1, guide = "none") +
  new_scale_fill() +
  geom_point(
    aes(
      fill = as.factor(df_ors$RS)),
    shape = 21,
    alpha = 0.6,
    size = 2,
    stroke = 0.5) + 
  scale_fill_manual(values=mypal) +
  scale_shape_manual(values = c(21), guide = "none") +
  xlab(paste("PCoA Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("PCoA Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2))) +
  xlim(-0.6,0.6) + 
  ylim(-0.44,0.6) + 
  theme_bw() + theme(
    panel.border = element_blank(),
    #panel.grid.major = element_line(colour = "darkgrey"),
    #panel.grid.minor = element_line(colour = "grey"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.85, 0.8),
    legend.text = element_text(size=12),
    legend.title = element_text(size=14),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size=12),
    axis.title = element_text(size=14)) + 
  labs(
    fill = "Reproductive system",
  ) + add_phylopic(img=solanum_pp,x = 0.28, y=-0.225, ysize = 0.15,col = "grey30") +
  add_phylopic(img=commelina_pp,x = 0.44, y=0, ysize = 0.15,col = "grey30") +
  add_phylopic(img=phoenix_pp,x = -0.52, y=0, ysize = 0.15,col = "grey30") +
  add_phylopic(img=coffea_pp,x = -0.07, y=-0.35, ysize = 0.125,col = "grey30") +
  add_phylopic(img=trithuria_pp,x = 0.07, y=0.55, ysize = 0.125,col = "grey30") +
  add_phylopic(img=cornus_pp,x = -0.35, y=-0.26, ysize = 0.125,col = "grey30") +
  add_phylopic(img=zea_pp,x = -0.39, y=0.38, ysize = 0.125,col = "grey30") +
  annotate("text", x=-0.07, y=-0.43, label= "paste(italic(Coffea), ' ',italic(arabica))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=-0.52, y=-0.09, label= "paste(italic(Phoenix), ' ',italic(dactylifera))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=0.28, y=-0.32, label= "paste(italic(Solanum), ' ',italic(dulcamara))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=0.44, y=-0.1, label= "paste(italic(Commelina), ' ',italic(communis))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=0.07, y=0.47, label= "paste(italic(Trithuria), ' ',italic(submersa))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=-0.35, y=-0.34, label= "paste(italic(Cornus), ' ',italic(florida))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=-0.39, y=0.3, label= "paste(italic(Zea), ' ',italic(mays))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("segment", #phoenix
           #linetype=2,
           linewidth=0.75,
           x=-0.46,
           xend=-0.3305978702, 
           y=0.03,
           yend=4.934202e-02, 
           color = "black",
           alpha=0.5) + 
  annotate("segment", #solanum
           #linetype=2,
           linewidth=0.75,
           x=0.25,
           xend=0.0966219059, 
           y=-0.2,
           yend=-3.542270e-02, 
           color = "black",
           alpha=0.5) +
  annotate("segment", #coffea
           #linetype=2,
           linewidth=0.75,
           x=-0.045,
           xend=-0.0139421126, 
           y=-0.28,
           yend=-0.1305380241, 
           color = "black",
           alpha=0.5) + 
  annotate("segment", #trithuria
           #linetype=2,
           linewidth=0.75,
           x=0.07,
           xend=0.1230980263,
           y=0.45,
           yend=3.744929e-01, 
           color = "black",
           alpha=0.5) + 
  annotate("segment", #cornus
           #linetype=2,
           linewidth=0.75,
           x=-0.3,
           xend=-0.1152021084,
           y=-0.21,
           yend=-0.1413927222, 
           color = "black",
           alpha=0.5) +
  annotate("segment", #zea
           #linetype=2,
           linewidth=0.75,
           x=-0.35,
           xend=-0.1219254049,
           y=0.35,
           yend=2.600799e-01, 
           color = "black",
           alpha=0.5) +
  annotate("segment", #Commelina
           #linetype=2,
           linewidth=0.75,
           x=0.38,
           xend=0.2257547576, 
           y=0,
           yend=0.0947481926, 
           color = "black",
           alpha=0.5)

# ggsave("figures/8_scatterplot_pcoa_coloured_original_RS.png",
#        width = 15,
#        height = 15,
#        units = 'cm')

####
## Plot taxonomy on PCOA ----
####

# load taxonomy
tax<-readRDS(file = here::here("outputs/taxonomy.rds"))

# check row order
table(rownames(df)==rownames(tax))

# empty column for common orders
tax$order_common<-NA

# top 10 most common orders in data
order_common<-names(sort(table(tax$order),decreasing = T)[1:8])
for(i in 1:length(order_common)){
  tax$order_common[grep(order_common[i],tax$order)]<-order_common[i]
}

# plot PCOA points on first two axes coloured by order
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

# ggsave("figures/8.0_scatterplot_pcoa_taxonomy.png",
#        width = 20,
#        height = 15,
#        units = 'cm')

####
## Correlations of traits with PCOA axes ----
####

# make empty matrix
corr_mat <- matrix(nrow=4,ncol=length(colnames(df)))

# loop through traits and calculate correlations with first four PCOA axes
for(i in 1:length(colnames(df))){
  
  if(is.numeric(df[,i])){
    
    corr_mat[1,i] <- cor(df[,i], dataset_pcoa$vectors[,1], method="pearson", use="complete.obs")
    corr_mat[2,i] <- cor(df[,i], dataset_pcoa$vectors[,2], method="pearson", use="complete.obs")
    corr_mat[3,i] <- cor(df[,i], dataset_pcoa$vectors[,3], method="pearson", use="complete.obs")
    corr_mat[4,i] <- cor(df[,i], dataset_pcoa$vectors[,4], method="pearson", use="complete.obs")
    
  } else {
    
    #https://stats.stackexchange.com/questions/119835/correlation-between-a-nominal-iv-and-a-continuous-dv-variable/124618#124618
    corr_mat[1,i] <- summary(lm(dataset_pcoa$vectors[,1] ~ df[,i]))$r.squared
    corr_mat[2,i] <- summary(lm(dataset_pcoa$vectors[,2] ~ df[,i]))$r.squared
    corr_mat[3,i] <- summary(lm(dataset_pcoa$vectors[,3] ~ df[,i]))$r.squared
    corr_mat[4,i] <- summary(lm(dataset_pcoa$vectors[,4] ~ df[,i]))$r.squared
    
  }
  
}

# set row/colnames
rownames(corr_mat)<-c("Axis.1","Axis.2","Axis.3","Axis.4")
colnames(corr_mat)<-colnames(df)

# transpose and order
corr_mat<-t(corr_mat)
corr_mat<-corr_mat[order(rownames(corr_mat)),]

# write output correlation table
write.csv(round(corr_mat,3),"outputs/8.0_correlation_pcoa_axes_traits.csv")
