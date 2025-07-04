rm(list=ls())

# load libraries
library(dplyr)
library(ggplot2)
library(cluster)
library(patchwork)
library(wesanderson)
library(vegan)

# load data set
df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

# dissimilarity matrix calculation
gower_df <- daisy(df,
                  metric = "gower" )
summary(gower_df)

# make into distance object
dataset_dist <- stats::as.dist(gower_df)

# run PCoA on distance matrix
dataset_pcoa <- ape::pcoa(dataset_dist)

# Recalculate relative eigenvalues by removing negative eigenvalues as in Mouillot et al.  
ev_pcoa <- dataset_pcoa$values$Eigenvalues
ev_pcoa_g0 <- ev_pcoa[ev_pcoa>0]
rel_ev_pcoa_g0 <- ev_pcoa_g0/sum(ev_pcoa_g0)

# Running NMDS in vegan (metaMDS)
nmds2 <-
  metaMDS(gower_df,
          distance = "gower",
          k = 2,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)

# Shepards test/goodness of fit
goodness(nmds2) # Produces a results of test statistics for goodness of fit for each point
stressplot(nmds2) # Produces a Shepards diagram

# PCoA scatterplot with density polygons
n2 <-ggplot(data.frame(nmds2$points), aes(x = MDS1, y = MDS2)) +
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
  xlab("nMDS Axis 1") +
  ylab("nMDS Axis 2") +
  xlim(-0.6,0.6) + 
  ylim(-0.44,0.6) + 
  theme_bw() + theme(
    panel.border = element_blank(),
    # panel.grid.major = element_line(colour = "darkgrey"),
    # panel.grid.minor = element_line(colour = "grey"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.text = element_text(size=13),
    legend.title = element_text(size=15),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size=15),
    axis.title = element_text(size=18)) + 
  labs(
    colour = "Flower sex",
    shape = "Woodiness"
  ) +  ggtitle("nMDS k = 2") 


# Running NMDS in vegan (metaMDS)
nmds3 <-
  metaMDS(gower_df,
          distance = "gower",
          k = 3,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)


# Shepards test/goodness of fit
goodness(nmds3) # Produces a results of test statistics for goodness of fit for each point
stressplot(nmds3) # Produces a Shepards diagram

# PCoA scatterplot with density polygons
n3 <-ggplot(data.frame(nmds3$points), aes(x = MDS1, y = MDS2)) +
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
  xlab("nMDS Axis 1") +
  ylab("nMDS Axis 2") +
  xlim(-0.6,0.6) + 
  ylim(-0.44,0.6) + 
  theme_bw() + theme(
    panel.border = element_blank(),
    # panel.grid.major = element_line(colour = "darkgrey"),
    # panel.grid.minor = element_line(colour = "grey"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.text = element_text(size=13),
    legend.title = element_text(size=15),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size=15),
    axis.title = element_text(size=18)) + 
  labs(
    colour = "Flower sex",
    shape = "Woodiness"
  ) +  ggtitle("nMDS k = 3") 


# Running NMDS in vegan (metaMDS)
nmds4 <-
  metaMDS(gower_df,
          distance = "gower",
          k = 4,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)


# Shepards test/goodness of fit
goodness(nmds4) # Produces a results of test statistics for goodness of fit for each point
stressplot(nmds4) # Produces a Shepards diagram

# PCoA scatterplot with density polygons
n4 <-ggplot(data.frame(nmds4$points), aes(x = MDS1, y = MDS2)) +
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
  xlab("nMDS Axis 1") +
  ylab("nMDS Axis 2") +
  xlim(-0.6,0.6) + 
  ylim(-0.44,0.6) + 
  theme_bw() + theme(
    panel.border = element_blank(),
    # panel.grid.major = element_line(colour = "darkgrey"),
    # panel.grid.minor = element_line(colour = "grey"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.text = element_text(size=13),
    legend.title = element_text(size=15),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size=15),
    axis.title = element_text(size=18)) + 
  labs(
    colour = "Flower sex",
    shape = "Woodiness"
  ) +  ggtitle("nMDS k = 4") 



# Running NMDS in vegan (metaMDS)
nmds5 <-
  metaMDS(gower_df,
          distance = "gower",
          k = 5,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)


# Shepards test/goodness of fit
goodness(nmds5) # Produces a results of test statistics for goodness of fit for each point
stressplot(nmds5) # Produces a Shepards diagram

# PCoA scatterplot with density polygons
n5 <-ggplot(data.frame(nmds5$points), aes(x = MDS1, y = MDS2)) +
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
  xlab("nMDS Axis 1") +
  ylab("nMDS Axis 2") +
  xlim(-0.6,0.6) + 
  ylim(-0.44,0.6) + 
  theme_bw() + theme(
    panel.border = element_blank(),
    # panel.grid.major = element_line(colour = "darkgrey"),
    # panel.grid.minor = element_line(colour = "grey"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.text = element_text(size=13),
    legend.title = element_text(size=15),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size=15),
    axis.title = element_text(size=18)) + 
  labs(
    colour = "Flower sex",
    shape = "Woodiness"
  ) +  ggtitle("nMDS k = 5") 


# Running NMDS in vegan (metaMDS)
nmds10 <-
  metaMDS(gower_df,
          distance = "gower",
          k = 10,
          maxit = 999, 
          trymax = 500,
          wascores = TRUE)


# Shepards test/goodness of fit
goodness(nmds10) # Produces a results of test statistics for goodness of fit for each point
stressplot(nmds10) # Produces a Shepards diagram

# PCoA scatterplot with density polygons
n10 <-ggplot(data.frame(nmds10$points), aes(x = MDS1, y = MDS2)) +
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
  xlab("nMDS Axis 1") +
  ylab("nMDS Axis 2") +
  xlim(-0.6,0.6) + 
  ylim(-0.44,0.6) + 
  theme_bw() + theme(
    panel.border = element_blank(),
    # panel.grid.major = element_line(colour = "darkgrey"),
    # panel.grid.minor = element_line(colour = "grey"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    legend.text = element_text(size=13),
    legend.title = element_text(size=15),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size=15),
    axis.title = element_text(size=18)) + 
  labs(
    colour = "Flower sex",
    shape = "Woodiness"
  ) +  ggtitle("nMDS k = 10") 


# PCoA scatterplot with density polygons
p1 <- ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2)) +
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
    # panel.grid.major = element_line(colour = "darkgrey"),
    # panel.grid.minor = element_line(colour = "grey"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.85, 0.8),
    legend.text = element_text(size=10),
    legend.title = element_text(size=12),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size=15),
    axis.title = element_text(size=18)) + 
  labs(
    colour = "Flower sex",
    shape = "Woodiness"
  ) + ggtitle("PCoA")  


n2 + p1
n3 + p1
n4 + p1
n5 + p1
n10 + p1

