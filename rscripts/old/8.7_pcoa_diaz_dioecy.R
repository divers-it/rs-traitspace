rm(list=ls())

#load libraries
library(dplyr)
library(ggplot2)
library(ggalt)
library(cluster)

#NOTE: UNCOMMENT TO LOAD WHEN NEEDED
#load most recent PCoA image
load("outputs/8.3_diaz_pcoa.Rdata")

#read in dioecy info
dio_df<-read.csv("data/diosis_24juin2013_cleaned.csv")
head(dio_df)

#get proportion of genus dioecious
dio_df$prop<-as.numeric(dio_df$n_dio)/as.numeric(dio_df$n)

#when proportion is >0.5 classify genus as dioecious
dio_df$classed_dio<-dio_df$prop>0.5
dio_genera<-dio_df$genus[dio_df$classed_dio==1]

#make label for plot
dio_label<-as.factor(as.numeric(diaz_cf$Genus%in%dio_genera))

#check order
rownames(diaz_cf)==rownames(data.frame(dataset_pcoa$vectors))

#plot PCOA points on first two axes coloured by dioecy and sized by plant height
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    aes(col=dio_label,size=diaz_cf$Plant_height_m),
    shape=16,
    alpha=0.5,
  ) + 
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

#plot PCOA points on axes 3&4 coloured by dioecy and sized by plant height
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.3, y = Axis.4)) +
  geom_point(
    aes(col=dio_label,size=diaz_cf$Plant_height_m),
    shape=16,
    alpha=0.5,
  ) + 
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

#Plot on seperate panels
#combined dfs
df<-data.frame(dataset_pcoa$vectors,dio_label)

#plot PCOA points on first two axes coloured by SI
ggplot(data.frame(df), aes(x = Axis.1, y = Axis.2)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", col=NA , n=100, bins=20) +
  scale_fill_distiller(palette = "Greys", direction = 1, guide = "none") +
  geom_point(
    aes(color=dio_label),
    shape=16,
    alpha=0.25,
    size=2
  ) + 
  theme_bw()+
  theme(legend.position = 'none') +
  facet_wrap(~dio_label,ncol=1)


#plot dioecious species on top
#plot PCOA points on first two axes coloured by SI
ggplot(data.frame(dataset_pcoa$vectors[dio_label==0,]), aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    aes(size=diaz_cf$Plant_height_m[dio_label==0]),
    col="grey",
    shape=16,
    alpha=0.5,
  ) + 
  geom_point(data=data.frame(dataset_pcoa$vectors[dio_label==1,]),
             aes(x = Axis.1, y = Axis.2,size=diaz_cf$Plant_height_m[dio_label==1]),
             col=harrypotter::hp(4,option="Gryffindor")[2],
             shape=16,
             alpha=0.5,
  ) + 
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2))) +
  theme_bw() +
  theme(legend.position = 'none')



###
# ---- PCoAs with Diaz traits on plots
###

#plot PCOA points on first two axes coloured by plant height
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    aes(col=log(diaz_cf$Plant_height_m)),
    shape=16,
    alpha=0.5,
  ) + 
  scale_color_gradientn(colors=rainbow(5)) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

#plot PCOA points on first two axes coloured by Nitrogen mass
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    aes(col=log(diaz_cf$Nmass_mg_g)),
    shape=16,
    alpha=0.5,
  ) + 
  scale_color_gradientn(colors=rainbow(5)) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

#plot PCOA points on first two axes colour by leaf mass
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    aes(col=log(diaz_cf$Leaf_area_mm2)),
    shape=16,
    alpha=0.5,
  ) + 
  scale_color_gradientn(colors=rainbow(5)) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))



#save pcoa image
#save.image("outputs/8.7_pcoa_diaz_dioecy.Rdata")