rm(list=ls())

#load libraries
library(dplyr)
library(ggplot2)
library(ggalt)
library(cluster)

#load Diaz data
diaz<-read.csv("data/diaz_species_mean_traits.csv")

#read in SI info
si_df<-read.csv("data/species_with_SIstatus.csv")

#combined genus species
si_df$name<-paste(si_df$Genus,si_df$Species,sep=" ")

#SI species in diaz data
diaz_c <- diaz[diaz$Species_name_standardized_against_TPL%in%si_df$name,]
diaz_c<-diaz_c[order(diaz_c$Species_name_standardized_against_TPL),]
str(diaz_c)

#filter by number of observations
#diaz_c<-diaz_c[diaz_c$Number_of_traits_with_values>3,]

#shared species in SI data
si_df<-si_df[si_df$name%in%diaz_c$Species_name_standardized_against_TPL,]
si_df<-si_df[order(si_df$name),]
str(si_df)

#check order
table(si_df$name==diaz_c$Species_name_standardized_against_TPL)

#only those columns for PCOA
diaz_pc<-diaz_c[,c(
  #"Adaptation_to_terrestrial_or_aquatic_habitats",
  #"Woodiness",
  #"Growth_Form",
  #"Succulence",
  #"Nutrition_type_parasitism",
  #"Nutrition_type_carnivory",
  #"Leaf_type",
  "Leaf_area_mm2",
  "Nmass_mg_g",
  "LMA_g_m2",
  "Plant_height_m",
  "Diaspore_mass_mg",
  #"SSD_observed_mg_mm3",
  #"LDMC_g_g",
  #"SSD_imputed_mg_mm3",
  "SSD_combined_mg_mm3"#,
  #"Number_of_traits_with_values"
)
]

#rownames
rownames(diaz_pc)<-diaz_c$Species_name_standardized_against_TPL

#impute missing data
imp <- missForest::missForest(diaz_pc, maxiter = 15, ntree = 100, variablewise = FALSE)
diaz_pc<-imp$ximp

#remove 0s by adding 1
diaz_pc$Diaspore_mass_mg<-diaz_pc$Diaspore_mass_mg + 1

#log transform and scale
diaz_pc<-log(diaz_pc)
diaz_pc<-scale(diaz_pc, center = T, scale = T)

#make data frame
diaz_pc<-as.data.frame(diaz_pc)

#UNCOMMENT TO RUN WHEN NEEDED
#dissimilarity matrix calculation
gower_df <- daisy(diaz_pc,
                  metric = "gower" )

summary(gower_df)

dataset_dist <- stats::as.dist(gower_df)
dataset_pcoa <- ape::pcoa(dataset_dist)

#Recalculate relative eigenvalues by removing negative eigenvalues as in Mouillot et al.  
ev_pcoa <- dataset_pcoa$values$Eigenvalues
ev_pcoa_g0 <- ev_pcoa[ev_pcoa>0]
rel_ev_pcoa_g0 <- ev_pcoa_g0/sum(ev_pcoa_g0)

#plot PCOA points on first two axes coloured by SI
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    aes(col=si_df$SI_status),
    shape=16,
    alpha=0.5,
    size=2
  ) + 
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

#plot PCOA points on first two axes coloured by SI
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.3, y = Axis.4)) +
  geom_point(
    aes(col=si_df$SI_status),
    shape=16,
    alpha=0.25,
    size=2
  ) + 
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[3],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[4],2)))

#check names
table(rownames(dataset_pcoa$vectors)==si_df$name)

#combined dfs
df<-data.frame(dataset_pcoa$vectors,si_df$SI_status)

#plot PCOA points on first two axes coloured by SI
ggplot(data.frame(df), aes(x = Axis.1, y = Axis.2)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", col=NA , n=100, bins=20) +
  scale_fill_distiller(palette = "Greys", direction = 1, guide = "none") +
  geom_point(
    aes(color=si_df.SI_status),
    shape=16,
    alpha=0.25,
    size=2
  ) + 
  theme_bw()+
  theme(legend.position = 'none') +
  facet_wrap(~si_df.SI_status,ncol=1)

#save pcoa image
#save.image("outputs/8.6_pcoa_diaz_si.Rdata")