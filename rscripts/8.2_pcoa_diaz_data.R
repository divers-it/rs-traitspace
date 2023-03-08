rm(list=ls())
library(dplyr)
library(ggplot2)
library(ggalt)

#load formatted DiveRS data
df<-readRDS(file = here::here("outputs/df_filt_trans.rds"))

#load Diaz data
diaz<-read.csv("data/diaz_species_mean_traits.csv")

#get rid of non-angiosperm
diaz[diaz$Phylogenetic_Group_General=="Angiosperm",]

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

#filter by number of observations
diaz_cf<-diaz_c[diaz_c$Number_of_traits_with_values>3,]

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
                "Plant_height_m",
                "Diaspore_mass_mg",
                #"SSD_observed_mg_mm3",
                #"LDMC_g_g",
                #"SSD_imputed_mg_mm3",
                "SSD_combined_mg_mm3"#,
                #"Number_of_traits_with_values"
)
]

#remove 0s by adding 1
diaz_pcf$Diaspore_mass_mg<-diaz_pcf$Diaspore_mass_mg + 1

#log transform and scale
diaz_pcf<-log(diaz_pcf)
diaz_pcf<-scale(diaz_pcf, center = T, scale = T)

#UNCOMMENT TO RUN WHEN NEEDED
#dissimilarity matrix calculation
library(cluster)
gower_df <- daisy(diaz_pcf,
                  metric = "gower" )

summary(gower_df)

dataset_dist <- stats::as.dist(gower_df)
dataset_pcoa <- ape::pcoa(dataset_dist)

#save/load pcoa image
#save.image("outputs/diaz_pcoa.Rdata")
#load("outputs/diaz_pcoa.Rdata")

#make pcoa vectors a data frame
df_pcoa<-data.frame(dataset_pcoa$vectors)

#check that pcoa is in same order as data
table(rownames(diaz_cf)==rownames(df_pcoa))

#set rownames
rownames(diaz_cf)<-diaz_cf$Species_name_standardized_against_TP
rownames(df_pcoa)<-diaz_cf$Species_name_standardized_against_TP

#make column for inclusion in DiveRS dataset
diaz_cf$divers<-as.numeric(diaz_cf$Species_name_standardized_against_TP%in%rownames(df))

#load substitutes
subs<-readRDS(file = here::here("outputs/sub_genera.rds"))

#make subsitutes 1 in diaz_cf$divers
diaz_cf$divers[which(diaz_cf$Species_name_standardized_against_TPL%in%subs$Species_name_standardized_against_TPL)]<-1

#reorder data frame based in divers inclusion
diaz_cf<-diaz_cf[order(diaz_cf$divers),]
tail(diaz_pcf)

#match order so divers are last
df_pcoa<-df_pcoa[match(diaz_cf$Species_name_standardized_against_TPL, rownames(df_pcoa)),]

#plot PCOA points on first two axes
ggplot(df_pcoa, aes(x = Axis.1, y = Axis.2, fill = as.factor(diaz_cf$divers))) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.75,
    size=2,
    stroke = 0.5
  ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))

ggsave("figures/scatter_pcoa_divers_in_diaz.png",
       width = 30,
       height = 30,
       units = 'cm')

#relative eigenvalues
eig_df<-data.frame(c(1:9),dataset_pcoa$values$Relative_eig[1:9])
colnames(eig_df)<-c("axis","relative_eigenvalue")
eig_df$axis<-as.character(eig_df$axis)

ggplot(eig_df, aes(x=axis, y=relative_eigenvalue)) + 
  geom_bar(stat = "identity")
ggsave("figures/rel_eig_pcoa_diaz.png")

##############################
# only DiveRS
##############################

df_pcoa_od<-df_pcoa[rownames(df_pcoa)%in%rownames(df),]

#add congenerics (this is done after filtering to only those with data for >x traits)
df_pcoa_od<-rbind(df_pcoa_od,df_pcoa[rownames(df_pcoa)%in%subs$Species_name_standardized_against_TPL,])

#plot PCOA points on first two axes
ggplot(df_pcoa_od, aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.75,
    size=2,
    stroke = 0.5
  ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))

#filter DiveRS dataset by species with Diaz data (enough to calc PCOA)
df_filt<-df[rownames(df)%in%rownames(df_pcoa_od),]
rownames(df_pcoa_od)==rownames(df_filt)

#plot PCOA points on first two axes
ggplot(df_pcoa_od, aes(x = Axis.1, y = Axis.2, fill = as.factor(df_filt$Woodiness))) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.75,
    size=2,
    stroke = 0.5
  ) + 
  geom_encircle(aes(fill = as.factor(df_filt$Woodiness)),
                s_shape = 1,
                expand = 0,
                alpha = 0.2,
                color = "black",
                show.legend = FALSE) +
    xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))

