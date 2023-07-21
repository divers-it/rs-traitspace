rm(list=ls())

#load libraries
library(dplyr)
library(ggplot2)
library(ggalt)
library(cluster)

#load most recent PCoA image
load("outputs/diaz_pcoa.Rdata")

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

#make data frame
diaz_pcf<-as.data.frame(diaz_pcf)

##UNCOMMENT TO RUN WHEN NEEDED
##dissimilarity matrix calculation
#gower_df <- daisy(diaz_pcf,
#                  metric = "gower" )
#
#summary(gower_df)
#
#dataset_dist <- stats::as.dist(gower_df)
#dataset_pcoa <- ape::pcoa(dataset_dist)
#
##Recalculate relative eigenvalues by removing negative eigenvalues as in Mouillot et al.  
#ev_pcoa <- dataset_pcoa$values$Eigenvalues
#ev_pcoa_g0 <- ev_pcoa[ev_pcoa>0]
#rel_ev_pcoa_g0 <- ev_pcoa_g0/sum(ev_pcoa_g0)
#
##save pcoa image
#save.image("outputs/diaz_pcoa.Rdata")

#make pcoa vectors a data frame
df_pcoa<-data.frame(dataset_pcoa$vectors)

#check that pcoa is in same order as data
table(rownames(diaz_pcf)==rownames(df_pcoa))

#set rownames
rownames(diaz_pcf)<-diaz_cf$Species_name_standardized_against_TP
rownames(df_pcoa)<-diaz_cf$Species_name_standardized_against_TP

#stats for one species BEFORE reorder
df_pcoa[rownames(df_pcoa)=="Zostera marina",c(1:3)]

#make column for inclusion in DiveRS dataset
diaz_pcf$divers<-as.numeric(rownames(diaz_pcf)%in%rownames(df))

#NOT RUN: Adding congenerics to increase numbers (spaces aren't comparable if added)
#load substitutes
#subs<-readRDS(file = here::here("outputs/sub_genera.rds"))

#make subsitutes 1 in diaz_cf$divers
#diaz_cf$divers[which(diaz_cf$Species_name_standardized_against_TPL%in%subs$Species_name_standardized_against_TPL)]<-1

#reorder data frame based in divers inclusion
diaz_pcf<-diaz_pcf[order(diaz_pcf$divers),]
head(diaz_pcf)
tail(diaz_pcf)

#match order so divers are last
df_pcoa<-df_pcoa[match(rownames(diaz_pcf), rownames(df_pcoa)),]
rownames(diaz_pcf)==rownames(df_pcoa)

#stats for one species AFTER reorder
df_pcoa[rownames(df_pcoa)=="Zostera marina",c(1:3)]

#plot PCOA points on first two axes with DiveRS species as blue triangles
ggplot(df_pcoa, aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    aes(pch = as.factor(diaz_pcf$divers),
        fill = as.factor(diaz_pcf$divers),
        size = as.factor(diaz_pcf$divers)),
    alpha=0.75,
    stroke = 0.5
  ) + scale_shape_manual(values=c(21,24)) +
  scale_size_manual(values=c(2,4)) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

ggsave("figures/scatter_pcoa_divers_in_diaz.png",
       width = 30,
       height = 30,
       units = 'cm')

#Make data frame of first 9 relative eigenvalues
eig_df<-data.frame(c(1:9),rel_ev_pcoa_g0[1:9])
colnames(eig_df)<-c("pcoa_axis","relative_eigenvalue")
eig_df$pcoa_axis<-as.character(eig_df$pcoa_axis)

ggplot(eig_df, aes(x=pcoa_axis, y=relative_eigenvalue)) + 
  geom_bar(stat = "identity")
ggsave("figures/rel_eig_pcoa_diaz.png")
