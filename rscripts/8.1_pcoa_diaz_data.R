rm(list=ls())

# load libraries
library(dplyr)
library(ggplot2)
library(ggalt)
library(cluster)
library(patchwork)

# NOTE: UNCOMMENT TO LOAD WHEN NEEDED
# load most recent PCoA image
# load("outputs/8.1_diaz_pcoa.Rdata")

# load formatted DiveRS data
df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

####
## Format Diaz et al. 2022 data ----
####

# load raw Diaz data
diaz<-read.csv("data/diaz_species_mean_traits.csv")

# get rid of non-angiosperm
diaz[diaz$Phylogenetic_Group_General=="Angiosperm",]

# clean columns
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
        "Number_of_traits_with_values",
        "Phylogenetic_Group_General"
        )
        ]

# filter by number of observations 
diaz_cf<-diaz_c[diaz_c$Number_of_traits_with_values>3,]

# only those columns for PCOA
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

# remove 0s that cause issues when log transforming by adding 1
diaz_pcf$Diaspore_mass_mg<-diaz_pcf$Diaspore_mass_mg + 1

# log transform and scale
diaz_pcf<-log(diaz_pcf)
diaz_pcf<-scale(diaz_pcf, center = T, scale = T)

# make data frame
diaz_pcf<-as.data.frame(diaz_pcf)

# UNCOMMENT TO RUN WHEN NEEDED
# dissimilarity matrix calculation
gower_df <- daisy(diaz_pcf,
                  metric = "gower" )

summary(gower_df)

dataset_dist <- stats::as.dist(gower_df)
dataset_pcoa <- ape::pcoa(dataset_dist)

# Recalculate relative eigenvalues by removing negative eigenvalues as in Mouillot et al.  
ev_pcoa <- dataset_pcoa$values$Eigenvalues
ev_pcoa_g0 <- ev_pcoa[ev_pcoa>0]
rel_ev_pcoa_g0 <- ev_pcoa_g0/sum(ev_pcoa_g0)

# save pcoa image
save.image("outputs/8.1_diaz_pcoa.Rdata")

# Make data frame of first 9 relative eigenvalues and plot
eig_df <- data.frame(c(1:9), rel_ev_pcoa_g0[1:9])
colnames(eig_df) <- c("pcoa_axis", "relative_eigenvalue")
eig_df$pcoa_axis <- as.character(eig_df$pcoa_axis)

ggplot(eig_df, aes(x=pcoa_axis, y=relative_eigenvalue)) + 
  geom_bar(stat = "identity")

# ggsave("figures/8.1_barplot_relative_eigenvalues_pcoa_diaz.png")

# make pcoa vectors a data frame
df_pcoa <- data.frame(dataset_pcoa$vectors)

#check that pcoa is in same order as data
table(rownames(diaz_pcf)==rownames(df_pcoa))

#set rownames
rownames(diaz_pcf)<-diaz_cf$Species_name_standardized_against_TP
rownames(df_pcoa)<-diaz_cf$Species_name_standardized_against_TP

# stats for one species BEFORE reorder
df_pcoa[rownames(df_pcoa)=="Zostera marina",c(1:3)]

# make column for inclusion in DiveRS dataset
diaz_pcf$divers<-as.numeric(rownames(diaz_pcf)%in%rownames(df))

# reorder data frame based in divers inclusion
diaz_pcf <- diaz_pcf[order(diaz_pcf$divers),]
head(diaz_pcf)
tail(diaz_pcf)

# match order so divers are last
df_pcoa <- df_pcoa[match(rownames(diaz_pcf), rownames(df_pcoa)), ]
table(rownames(diaz_pcf) == rownames(df_pcoa))

# stats for one species AFTER reorder
df_pcoa[rownames(df_pcoa) == "Zostera marina", c(1:3)]

####
## Read in Dimensionality (script 9) output ----
#### 

# Prepare Data
load(file = here::here("outputs/9_dimensionality", "res_for_model.RData"))

# get files and names
filenames <- list.files(path = here::here("outputs/9_dimensionality"), pattern = "_res.rds$",
                        full.names = FALSE)

files <- list.files(path = here::here("outputs/9_dimensionality"), pattern = "_res.rds$",
                    full.names = TRUE)

# read in results of all analyses
list_res <- lapply(files, function(x) readRDS(x))

# get data frame of trait space quality metrics
res_for_graph_dim <- data.frame(do.call(rbind, lapply(1:length(list_res),
                                                      function(i) {
                                                        res <- data.frame(list_res[[i]][[3]])
                                                        res$taxa <- gsub("_res\\.rds", "", filenames)[i]
                                                        return(res)
                                                      })))


# make empty columns for species and number of traits
res_for_graph_dim$SP    <- NA 
res_for_graph_dim$trait <- NA 

# fill columns
for (i in 1:nrow(res_for_graph_dim)){ 
  
  res_for_graph_dim$SP[i]    <- res_for_model[rownames(res_for_model) %in% res_for_graph_dim$taxa[i], ]$S
  res_for_graph_dim$trait[i] <- res_for_model[rownames(res_for_model) %in% res_for_graph_dim$taxa[i], ]$Nb_trait
  res_for_graph_dim$elbow[i] <- res_for_model[rownames(res_for_model) %in% res_for_graph_dim$taxa[i], ]$Nb_dim_AUC_elbow
}


# make columns for elbow method and associated AUC values
res_for_graph_dim$selec_elbow_graph <- NA
res_for_graph_dim$AUCwhenelbow      <- NA

# fill columns
for (i in 1:nrow(res_for_graph_dim)) { 
  
  if (res_for_graph_dim$dim[i] == res_for_graph_dim$elbow[i]) {
    
    res_for_graph_dim$selec_elbow_graph[i] <- res_for_graph_dim$dim[i]
    res_for_graph_dim$AUCwhenelbow[i]      <- res_for_graph_dim$AUC[i]
  }
}


# Sort plots by number of species or traits

res_for_graph_dim$taxa <- factor(res_for_graph_dim$taxa, 
                                 levels = unique(res_for_graph_dim$taxa[order(res_for_graph_dim$SP)]), 
                                 ordered = TRUE)

res_for_graph_dim <- res_for_graph_dim[order(res_for_graph_dim$taxa, decreasing = FALSE), ]

####
##  Figure 4a PCoA ----
####

# set colours
mypal <- harrypotter::hp(n = 4,option='slytherin')[c(1,4)]

diaz_pcf$divers[diaz_pcf$divers == 0] <- "Diaz et al. only"
diaz_pcf$divers[diaz_pcf$divers == 1] <- "Shared species"

# plot PCOA points on first two axes with DiveRS species as green circles
s1 <- ggplot(df_pcoa, aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    aes(pch = as.factor(diaz_pcf$divers),
        fill = as.factor(diaz_pcf$divers),
        col = as.factor(diaz_pcf$divers),
        size = as.factor(diaz_pcf$divers),
        alpha = as.factor(diaz_pcf$divers)),
    stroke = 0.5
  ) + 
  scale_fill_manual("Data set",values=mypal) +
  scale_color_manual("Data set",values=mypal) +
  scale_shape_manual("Data set",values=c(21,21)) +
  scale_size_manual("Data set",values=c(2,4)) +
  scale_alpha_manual("Data set",values=c(0.6,0.8)) +
  xlab(paste("PCoA Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("PCoA Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2))) + 
  theme_bw() +
  theme(
    panel.border = element_blank(),
    #panel.grid.major = element_line(colour = "darkgrey"),
    #panel.grid.minor = element_line(colour = "grey"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.85, 0.9),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14))

s1

####
# ---- Figure 4b Dimensions ----
####

# Subset with shared data only
df_shared <- res_for_graph_dim[grep("shared",res_for_graph_dim$taxa),]
df_shared$taxa<-as.character(df_shared$taxa)
df_shared$taxa[grep("Diaz_shared_2022",df_shared$taxa)] <- "Diaz et al. 2022"
df_shared$taxa[grep("DiveRS_shared_2023",df_shared$taxa)] <- "This study"

# set colours
mypal <- harrypotter::hp(n = 4,option='slytherin')[c(2,3)]

# plot
l1 <- ggplot(df_shared, aes(x = dim, y = AUC, colour = taxa)) + 
  stat_summary(fun = "mean", geom = "line", size = 1, alpha = 0.4) +
  stat_summary(fun = "mean", size = 0.88) +
  labs(x = "Number of PCoA axes") + 
  labs(y = "Quality of species trait space (AUC)") +
  facet_wrap(~ taxa,ncol = 2) + theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_color_manual(values=mypal)

# add label boxes about data
l2 <- l1 + geom_label(data = df_shared, 
             aes(label = paste0(taxa, "\n",
                                "No. species = ", SP , "\n",
                                "No. traits = ", trait),
                 y = 0.3, x = 11), size = 4, hjust = 0) +
  
  scale_y_continuous(breaks = seq(0.1, 1, 0.2))

l2

####
## Construct Figure 4 ----
####

patch <- ( s1 / l2 ) + plot_layout(heights=c(2, 1))

patch + plot_annotation(tag_levels = 'a',tag_prefix="(",tag_suffix=")") & 
  theme(plot.tag = element_text(size = 14))

ggsave("figures/figure_4_Diaz_comparison.png",width=10,height=15)
