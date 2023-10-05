rm(list=ls())

#load libraries
library(dplyr)
library(ggplot2)

#set plot margins
par(mar=c(3,3,3,3))

#load formatted DiveRS data
df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

#load Diaz data
diaz<-read.csv("data/diaz_species_mean_traits.csv")

#clean columns
diaz_c<-diaz[,c("Species_name_standardized_against_TPL",
                "Genus",
                "Family",
                "Adaptation_to_terrestrial_or_aquatic_habitats",
                "Woodiness",
                "Growth_Form",
                "Succulence",
                "Nutrition_type_parasitism",
                "Nutrition_type_carnivory",
                "Leaf_type",
                "Leaf_area_mm2",
                "Nmass_mg_g",
                "LMA_g_m2",
                "Plant_height_m",
                "Diaspore_mass_mg",
                "SSD_observed_mg_mm3",
                "LDMC_g_g",
                "SSD_imputed_mg_mm3",
                "SSD_combined_mg_mm3",
                "Number_of_traits_with_values"
)
]

#filter by presence in DiveRS dataset
diaz_cf<-diaz_c[diaz_c$Species_name_standardized_against_TPL%in%rownames(df),]

#which species in DiveRS are not present in Diaz dataset
setdiff(rownames(df),diaz_c$Species_name_standardized_against_TPL)

#filter by number of observations
diaz_cf<-diaz_cf[diaz_cf$Number_of_traits_with_values>3,]

#set row names
rownames(diaz_cf)<-diaz_cf$Species_name_standardized_against_TPL

#only those columns for PCOA
diaz_pcf<-diaz_cf[,c(
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
  #"Plant_height_m", #removed as already in our data
  #"Diaspore_mass_mg",
  #"SSD_observed_mg_mm3",
  #"LDMC_g_g",
  #"SSD_imputed_mg_mm3",
  "SSD_combined_mg_mm3"#,
  #"Number_of_traits_with_values"
)
]

#check for empty species
sort(rowSums(is.na(diaz_pcf)))
head(diaz_pcf)

# #remove 0s by adding 1
# diaz_pcf$Diaspore_mass_mg<-diaz_pcf$Diaspore_mass_mg + 1

#log transform and scale
diaz_pcf<-log(diaz_pcf)
diaz_pcf<-scale(diaz_pcf, center = T, scale = T)

#filter by presence in Diaz dataset
df<-df[rownames(df)%in%rownames(diaz_pcf),]

#combine data sets
rownames(df)==rownames(diaz_pcf)
combined_df<-cbind(df,diaz_pcf)

#plot(combined_df$Diaspore_mass_mg,df$seedMass)
#combined_df<-combined_df[,colnames(combined_df)!="seedMass"]


#dissimilarity matrix calculation
library(cluster)
gower_df <- daisy(combined_df,
                  metric = "gower" )
summary(gower_df)

#NOT RUN:
#Replace NA with 0 to calculate Euclidean distance
#combined_df[is.na(combined_df)] <- 0
#gower_df <- dist(combined_df)

#set as distance
dataset_dist <- stats::as.dist(gower_df)

#run PCoA
dataset_pcoa <- ape::pcoa(dataset_dist)

#Recalculate relative eigenvalues by removing negative eigenvalues as in Mouillot et al.  
ev_pcoa <- dataset_pcoa$values$Eigenvalues
ev_pcoa_g0 <- ev_pcoa[ev_pcoa>0]
rel_ev_pcoa_g0 <- ev_pcoa_g0/sum(ev_pcoa_g0)


#make pcoa vectors a data frame
df_pcoa<-data.frame(dataset_pcoa$vectors)
rownames(df_pcoa)<-rownames(dataset_pcoa$vectors)

# #filter DiveRS dataset by species with Diaz data (enough to calc PCOA)
# df_filt<-df[rownames(df)%in%rownames(combined_df),]
# rownames(diaz_pcf)==rownames(df_filt)
# 
# #save df with only species shared between divers and Diaz et al. for trait space quality comparison
# saveRDS(df_filt, file = here::here("outputs/8.4_df_filt_trans_shared.rds"))

#plot PCOA points on first two axes coloured by woodiness
p1 <- ggplot(df_pcoa, aes(x = Axis.1, y = Axis.2, fill = as.factor(combined_df$FlowerSex))) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.75,
    size=3,
    stroke = 0.5
  ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

#plot PCOA points on first two axes coloured by woodiness
ggplot(df_pcoa, aes(x = Axis.1, y = Axis.2, fill = as.factor(combined_df$Pollination))) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.75,
    size=3,
    stroke = 0.5
  ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

#plot PCOA points on first two axes coloured by woodiness
p2 <- ggplot(df_pcoa, aes(x = Axis.1, y = Axis.2, fill = as.factor(combined_df$Woodiness))) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.75,
    size=3,
    stroke = 0.5
  ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

#plot PCOA points on first two axes coloured by woodiness
ggplot(df_pcoa, aes(x = Axis.1, y = Axis.2, fill = combined_df$seedMass)) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.75,
    size=3,
    stroke = 0.5
  ) +
  scale_fill_gradientn(colours=rainbow(4)) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))


d1 <- ggplot(df_pcoa, aes(x = Axis.1, y = Axis.2, fill = combined_df$Leaf_area_mm2)) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.75,
    size=3,
    stroke = 0.5
  ) + 
  scale_fill_gradientn(colours=rainbow(4)) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))




d2 <- ggplot(df_pcoa, aes(x = Axis.1, y = Axis.2, fill = combined_df$Nmass_mg_g)) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.75,
    size=3,
    stroke = 0.5
  ) + 
  scale_fill_gradientn(colours=rainbow(4)) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

d3 <- ggplot(df_pcoa, aes(x = Axis.1, y = Axis.2, fill = combined_df$LMA_g_m2)) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.75,
    size=3,
    stroke = 0.5
  ) + 
  scale_fill_gradientn(colours=rainbow(4)) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

d4 <- ggplot(df_pcoa, aes(x = Axis.1, y = Axis.2, fill = combined_df$SSD_combined_mg_mm3)) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.75,
    size=3,
    stroke = 0.5
  ) + 
  scale_fill_gradientn(colours=rainbow(4)) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))


library(patchwork)
 (p1 + p2) / (d1 + d2) / (d3 + d4)


#Make data frame of first 9 relative eigenvalues
eig_df<-data.frame(c(1:9),rel_ev_pcoa_g0[1:9])
colnames(eig_df)<-c("pcoa_axis","relative_eigenvalue")
eig_df$pcoa_axis<-as.character(eig_df$pcoa_axis)

par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
boxplot(combined_df$Leaf_area_mm2~combined_df$FlowerSex)
boxplot(combined_df$Nmass_mg_g~combined_df$FlowerSex)
boxplot(combined_df$LMA_g_m2~combined_df$FlowerSex)
boxplot(combined_df$SSD_combined_mg_mm3~combined_df$FlowerSex)


#DiveRS data distance matrix
divers_dist <- daisy(df,
                     metric = "gower" )

#check names
labels(dataset_dist)==labels(divers_dist)

#compare pairwaise distances of the DiveRS and Diaz distance matrices
par(mar=c(5,5,5,5))
plot(dataset_dist,divers_dist,xlim=c(0,1),ylim=c(0,1))


#linear model
summary(lm(dataset_dist~divers_dist))

