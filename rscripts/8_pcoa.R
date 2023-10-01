rm(list=ls())

#load libraries
library(dplyr)
library(ggplot2)
library(cluster)
library(patchwork)
library(wesanderson)

#read in imputed data set
df<-read.csv("outputs/imputed_with_phylo.csv",row.names=1,stringsAsFactors = TRUE)

#dissimilarity matrix calculation
gower_df_no_miss <- daisy(df,
                          metric = "gower" )
summary(gower_df_no_miss)

#load data set
df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

#dissimilarity matrix calculation
gower_df <- daisy(df,
                  metric = "gower" )
summary(gower_df)

#check names
labels(gower_df)==gsub("_"," ",labels(gower_df_no_miss))

#linear model
mod<-summary(lm(gower_df~gower_df_no_miss))
mod$r.squared

#compare pairwaise distances of matrices with missing data and with imputed
png("figures/8_scatterplot_dist_missing_vs_imputed.png",width = 500,height = 500)
plot(gower_df,gower_df_no_miss,xlim=c(0,1),ylim=c(0,1),xlab="Original",ylab="Imputed")
abline(0,1,lty=2,col="red",lwd=2)
text(x=0.15, y=0.9, labels=paste("R-squared =",round(mod$r.squared,3)))
dev.off()

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

#missing data per row (species) for plot
missDat<-rowSums(apply(is.na(df),2,as.numeric))/ncol(df)

#plot PCOA points on first two axes coloured by missing data
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

#save plot
ggsave("figures/8_scatterplot_pcoa_missing.png",
       width = 20,
       height = 15,
       units = 'cm')

#Make data frame of first 9 relative eigenvalues
eig_df<-data.frame(c(1:9),rel_ev_pcoa_g0[1:9])
colnames(eig_df)<-c("pcoa_axis","relative_eigenvalue")
eig_df$pcoa_axis<-as.character(eig_df$pcoa_axis)

#prop variation first two axes explain
sum(rel_ev_pcoa_g0[1:2])

#plot barplot of first 9 relative eigenvalues
ggplot(eig_df, aes(x=pcoa_axis, y=relative_eigenvalue)) + 
  geom_bar(stat = "identity")
ggsave("figures/8_barplot_relative_eigenvalues_pcoa.png")

#plot points on first two axes with point style changed by two variables
#reproductive systems (color) and woodiness (shape)
p1 <- ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.3, y = Axis.4, fill = as.factor(df$DispersalMode), shape=as.factor(df$Woodiness))) +
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

p1

# mating system (colour) and woodiness (shape)
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
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2))) +
  theme(legend.position="bottom",
        legend.title = element_blank(),
        legend.box="vertical", 
        legend.margin=margin())

#combine plots
p1 | p2

#save
ggsave("figures/8_scatterplot_pcoa_coloured_by_traits.png",
       width = 30,
       height = 15,
       units = 'cm')


library(rphylopic)

#choose from the different ones for a grass species
img <- pick_phylopic(name = "Oryza sativa", n = 5)
# Get a single image uuid
uuid <- get_uuid(name = "Oryza sativa", n = 1)
# Get the image for that uuid
oryza_pp <- get_phylopic(uuid = uuid)

uuid <- get_uuid(name = "Zea mays", n = 4)
zea_pp <- get_phylopic(uuid = uuid[3])

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

#add new factor level e.g.  None 
df$Woodiness = factor(df$Woodiness, levels=c(levels(df$Woodiness), "None"))
#df$SexualSystem = factor(df$SexualSystem, levels=c(levels(df$SexualSystem), "None"))

#convert all NA's to None
df$Woodiness[is.na(df$Woodiness)] = "None"
#df$SexualSystem = factor(df$SexualSystem, levels=c(levels(df$SexualSystem), "None"))

####
# ---- Figure 1: PCoA with traits and density ----
####

#placement of species with phylopics
tmp<-readRDS("outputs/phylopic_uuids.rds")
dataset_pcoa$vectors[rownames(dataset_pcoa$vectors)%in%names(tmp),][,c(1:2)]


#PCoA scatterplot with density polygons
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
    legend.text = element_text(size=13),
    legend.title = element_text(size=15),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size=15),
    axis.title = element_text(size=18)) + 
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
  

ggsave("figures/8_scatterplot_pcoa_density_polygon.png",
       width = 25,
       height = 25,
       units = 'cm')

#Plot density raster with PCoA scatterplot points overlain
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2)) + 
  stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
  scale_fill_distiller(palette=3, direction=1) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    legend.position='none'
  ) +
  geom_point(
    aes(color = as.factor(df$Mating), shape=as.factor(df$Woodiness)),
    #    shape=21,
    alpha=0.5,
    size=4,
    stroke = 0.5
  )

# ---- Taxonomy ----

#load taxonomy
tax<-readRDS(file = here::here("outputs/taxonomy.rds"))

#check row order
rownames(df)==rownames(tax)

#empty column for common orders
tax$order_common<-NA

#top 10 most common orders in data
order_common<-names(sort(table(tax$order),decreasing = T)[1:8])

for(i in 1:length(order_common)){
  tax$order_common[grep(order_common[i],tax$order)]<-order_common[i]
}

#plot PCOA points on first two axes coloured by order
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

ggsave("figures/8_scatterplot_pcoa_taxonomy.png",
       width = 20,
       height = 15,
       units = 'cm')
