rm(list=ls())

# Load Project Addins (R Functions and Packages)
# devtools::load_all()

# not loading ggplot?
library(ggplot2)

# load_all doesn't load the files in R/ ?
for (nm in list.files("R", pattern = "[.]R$")) {
   source(file.path("R", nm))
}

# Make output folder
system(paste("mkdir",here::here("outputs/9_dimensionality")))

# Variable to iterate over
percent_list <- seq(0.1, 0.8, by = 0.1)

####
## Run Original DiveRS dataset ----
####

# Read data set
dataset <- readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

# Run dimensionality analyses
run_analysis(dataset, name = "DiveRS_2023")
rm(list = "dataset")

####
## Run DiveRS one-hot dataset ----
####

# read data set
dataset <- readRDS(file = here::here("outputs/one_hot_6_df_filt_trans.rds"))

# run analysis
run_analysis(dataset, name = "DiveRS_one_hot_2023")
rm(list = "dataset")

####
## Run DiveRS imputed dataset with no NA ----
####

# read data set
dataset <- read.csv(file = here::here("outputs/imputed_with_phylo.csv"),row.names = 1,stringsAsFactors = TRUE)

# run analysis
run_analysis(dataset, name = "DiveRS_imputed_2023")
rm(list = "dataset")


####
## Run shared species with Diaz data----
####

# load libraries
library(dplyr)
library(ggplot2)
library(cluster)

# load formatted DiveRS data
df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

####
### Format Diaz et al. 2022 data ----
####

# load Diaz data
diaz<-read.csv("data/diaz_species_mean_traits.csv")

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
                "Number_of_traits_with_values"
)
]

# filter by presence in DiveRS dataset
diaz_cf<-diaz_c[diaz_c$Species_name_standardized_against_TPL%in%rownames(df),]

# which species in DiveRS are not present in Diaz dataset
setdiff(rownames(df),diaz_c$Species_name_standardized_against_TPL)

# filter by number of observations
diaz_cf<-diaz_cf[diaz_cf$Number_of_traits_with_values>3,]

# set row names
rownames(diaz_cf)<-diaz_cf$Species_name_standardized_against_TPL

# only those columns for PCOA
diaz_pcf<-diaz_cf[,c(
  # "Adaptation_to_terrestrial_or_aquatic_habitats",
  # "Woodiness",
  # "Growth_Form",
  # "Succulence",
  # "Nutrition_type_parasitism",
  # "Nutrition_type_carnivory",
  # "Leaf_type",
  "Leaf_area_mm2",
  "Nmass_mg_g",
  "LMA_g_m2",
  "Plant_height_m",
  "Diaspore_mass_mg",
  # "SSD_observed_mg_mm3",
  # "LDMC_g_g",
  # "SSD_imputed_mg_mm3",
  "SSD_combined_mg_mm3"# ,
  # "Number_of_traits_with_values"
)
]

# check for empty species
sort(rowSums(is.na(diaz_pcf)))
head(diaz_pcf)

# remove 0s by adding 1
diaz_pcf$Diaspore_mass_mg<-diaz_pcf$Diaspore_mass_mg + 1

# log transform and scale
diaz_pcf<-log(diaz_pcf)
diaz_pcf<-scale(diaz_pcf, center = T, scale = T)

####
### Make and plot PCOA ----
####

# dissimilarity matrix calculation
gower_df <- daisy(diaz_pcf,
                  metric = "gower" )
summary(gower_df)

# set as distance
dataset_dist <- stats::as.dist(gower_df)

# run PCoA
dataset_pcoa <- ape::pcoa(dataset_dist)

# Recalculate relative eigenvalues by removing negative eigenvalues as in Mouillot et al.  
ev_pcoa <- dataset_pcoa$values$Eigenvalues
ev_pcoa_g0 <- ev_pcoa[ev_pcoa>0]
rel_ev_pcoa_g0 <- ev_pcoa_g0/sum(ev_pcoa_g0)

# set rownames
rownames(diaz_pcf)<-diaz_cf$Species_name_standardized_against_TP

# rename data set
dataset <- diaz_pcf

# run analysis
run_analysis(dataset, name = "Diaz_shared_2022")
rm(list = "dataset")

# make pcoa vectors a data frame
df_pcoa<-data.frame(dataset_pcoa$vectors)
rownames(df_pcoa)<-rownames(dataset_pcoa$vectors)

# filter DiveRS dataset by species with Diaz data
df_filt<-df[rownames(df)%in%rownames(diaz_pcf),]
rownames(diaz_pcf)==rownames(df_filt)

####
## Run shared species with DiveRS data ----
####

# rename data set
dataset <- df_filt

# run analysis
run_analysis(dataset, name = "DiveRS_shared_2023")
rm(list = "dataset")

# plot PCOA points on first two axes coloured by woodiness
ggplot(df_pcoa, aes(x = Axis.1, y = Axis.2, fill = as.factor(df_filt$Woodiness))) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.75,
    size=3,
    stroke = 0.5
  ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

# Make data frame of first 9 relative eigenvalues
eig_df<-data.frame(c(1:9),rel_ev_pcoa_g0[1:9])
colnames(eig_df)<-c("pcoa_axis","relative_eigenvalue")
eig_df$pcoa_axis<-as.character(eig_df$pcoa_axis)

ggplot(eig_df, aes(x=pcoa_axis, y=relative_eigenvalue)) + 
  geom_bar(stat = "identity")

# rename diaz data
diaz_dist<-gower_df

####
## Compare DiveRS and Diaz et al. 2022 ----
####

# DiveRS data distance matrix
divers_dist <- daisy(df_filt,
                     metric = "gower" )

# check names
table(labels(diaz_dist)==labels(divers_dist))

###
# Figure S4: Distances Diaz vs DiveRS ----
###

# compare pairwaise distances of the DiveRS and Diaz distance matrices
png(filename = "figures/figure_S4_distances_diaz_vs_divers.png",width = 750,height = 750)
plot(x=divers_dist,y=diaz_dist,xlim=c(0,1),ylim=c(0,1),xlab="Distance (this study)",ylab="Distance (Diaz et al. 2022)")
abline(0,1,col='red',lty=2,lwd=2)
dev.off()

# mantel test
vegan::mantel(as.matrix(diaz_dist),as.matrix(divers_dist))


## Import Results ----

# This will import all of the results present in the outputs folder
# (from each different dimensionality run)
files <- list.files(path = here::here("outputs/9_dimensionality"), pattern = "_res.rds$",
                    full.names = TRUE)
list_res <- lapply(files, function(x) readRDS(x))

files <- list.files(path = here::here("outputs/9_dimensionality"), pattern = "_cor.rds$",
                    full.names = TRUE)
list_res_cor <- lapply(files, function(x) readRDS(x))

files <- list.files(path = here::here("outputs/9_dimensionality"), pattern = "_single.rds$",
                    full.names = TRUE)
list_res_sigle <- lapply(files, function(x) readRDS(x))

files <- list.files(path = here::here("outputs/9_dimensionality"), pattern = "_pcoa.rds$",
                    full.names = TRUE)
list_res_pcoa <- lapply(files, function(x) readRDS(x))

filenames <- list.files(path = here::here("outputs/9_dimensionality"), pattern = "_pcoa.rds$",
                        full.names = FALSE)

## Create Full Table from all results ----

res_for_model <- t(data.frame(do.call(cbind, lapply(1:length(list_res), 
                                                    function(i) {
                                                      
                                                      res <- data.frame(list_res[[i]][[1]])
                                                      
                                                      # Add information of PCOA Axis and correlation between traits
                                                      res_pcoa1 <- data.frame(list_res_pcoa[[i]]$values[1, 2])
                                                      colnames(res_pcoa1) <- paste0("Relative_eig","_axe1")
                                                      
                                                      res_pcoa2 <- data.frame(list_res_pcoa[[i]]$values[2, 2])
                                                      colnames(res_pcoa2) <- paste0("Relative_eig","_axe2")
                                                      
                                                      res_pcoa3 <- data.frame(list_res_pcoa[[i]]$values[3, 2])
                                                      colnames(res_pcoa3) <- paste0("Relative_eig","_axe3") 
                                                      
                                                      res_cor<- data.frame(list_res_cor[[i]])
                                                      
                                                      # Add information on redundancy per cluster
                                                      
                                                      res_sigle <- t(data.frame(table(list_res_sigle[[i]]$cluster_core))[2:4, 2])
                                                      colnames(res_sigle) <- c("NbS_Cluster1", "NbS_Cluster2", "NbS_Cluster3")
                                                      res <- t(cbind(res, res_pcoa1, res_pcoa2, res_pcoa3, res_cor, res_sigle))
                                                      
                                                      return(res)
                                                    }))))

res_for_model <- as.data.frame(res_for_model)
rownames(res_for_model) <- gsub("_pcoa\\.rds", "", filenames)

# make out of 1
res_for_model$Percentage_lostAUC_depleted0.5  <- res_for_model$Percentage_lostAUC_depleted0.5  / 100
res_for_model$Percentage_lostAUC_depleted0.20 <- res_for_model$Percentage_lostAUC_depleted0.20 / 100

# rename columns
colnames(res_for_model)[14] <- "rowAUClostwhen50percTraitdepleted"
colnames(res_for_model)[15] <- "rowAUClostwhen20percTraitdepleted"

# save results
save(res_for_model, file = here::here("outputs/9_dimensionality", "res_for_model.RData"))

####
## Quality of species trait spaces (AUC) ----
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
                                                        res$taxa <- gsub("_[0-9][0-9]*_res\\.rds", "", filenames)[i]
                                                        return(res)
                                                      })))

rownames(res_for_model) <- gsub("_[0-9][0-9]*", "", rownames(res_for_model))

# make empty columns for species and number of traits
res_for_graph_dim$SP    <- NA 
res_for_graph_dim$trait <- NA 
res_for_graph_dim$elbow <- NA

# remove dataset we don't plot here
res_for_graph_dim <- res_for_graph_dim[res_for_graph_dim$taxa != "proteus",]

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

###
# Figure S5: Dimensionality ----
###

# Plot
p <- ggplot(res_for_graph_dim, aes(x = dim, y = AUC, colour = taxa)) +
  stat_summary(
    fun = "mean",
    geom = "line",
    size = 1,
    alpha = 0.4
  ) +
  stat_summary(fun = "mean", size = 0.88) +
  labs(x = "Number of dimensions(PCoA axes)") +
  labs(y = "Quality of species trait space (AUC)") +
  facet_wrap( ~ taxa, ncol = 6) + theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    panel.background = element_blank(),
    legend.position = "none"
  ) +
  harrypotter::scale_colour_hp_d(option = "LunaLovegood") +
    geom_label(
    data = res_for_graph_dim,
    aes(
      label = paste0("data set = ", taxa, "\n",
                     "# species = ", SP , "\n",
                     "# traits = ", trait),
      y = 0.95,
      x = 0
    ),
    size = 2.1,
    hjust = 0
  ) +
  scale_y_continuous(breaks = seq(0.1, 1, 0.2))


#save plot
grDevices::png(file = here::here("figures", "figure_S5_dimensionality.png"),width = 6250,height=2500,res=500)
print(p)
dev.off()

####
## Not used: Trait omission vs trait space quality ----
####

# Prepare Data

load(file = here::here("outputs/9_dimensionality", "res_for_model.RData"))

# get files and names
filenames <- list.files(path = here::here("outputs/9_dimensionality"), pattern = "_miss.rds$",
                        full.names = FALSE)

files <- list.files(path = here::here("outputs/9_dimensionality"), pattern = "_miss.rds$",
                    full.names = TRUE)

# read in files
list_res <- lapply(files, function(x) readRDS(x))

# make missingness data frame
res_for_graph_miss <- na.omit(data.frame(do.call(rbind, lapply(1:length(list_res),
                                                               function(i) {
                                                                 res <- data.frame(list_res[[i]])
                                                                 res$taxa <- gsub("_miss\\.rds", "", filenames)[i]
                                                                 return(res)
                                                               }))))

# number of species and number of traits columns.
res_for_graph_miss$SP    <- NA 
res_for_graph_miss$trait <- NA 

for (i in 1:nrow(res_for_graph_miss)){ 
  
  res_for_graph_miss$SP[i]    <- res_for_model[rownames(res_for_model) %in% res_for_graph_miss$taxa[i], ]$S
  res_for_graph_miss$trait[i] <- res_for_model[rownames(res_for_model) %in% res_for_graph_miss$taxa[i], ]$Nb_trait
}

# dataset labels
res_for_graph_miss$taxa <- factor(res_for_graph_miss$taxa, 
                                  levels  = unique(res_for_graph_miss$taxa[order(res_for_graph_miss$SP)]), 
                                  ordered = TRUE)

# reorder by data set
res_for_graph_miss <- res_for_graph_miss[order(res_for_graph_miss$taxa, decreasing = FALSE), ]


# plot
p2 <- ggplot(res_for_graph_miss, aes(x = miss_percent * 100, y = AUC,
                                     fill = as.factor(miss_percent * 100))) + 
  geom_boxplot() +
  theme_bw() +
  labs(x = "Traits omission (%)", 
       y = "Quality of species trait space (AUC)", size = 14) +
  harrypotter::scale_fill_hp_d(option = "ronweasley2",direction = -1)+
  facet_wrap(~ taxa, ncol = 6) + 
  theme_bw() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  theme(strip.background = element_blank(),
        strip.text.x     = element_blank(),
        panel.grid.major = element_blank(), 
        panel.background = element_blank(),
        axis.title.x     = element_text( size=14, face="bold"),
        axis.title.y     = element_text( size=14, face="bold"),
        legend.position  = "none") + 
  geom_label(data = res_for_graph_miss, 
             aes(label = paste0("data set = ", taxa),
                 y = 0.95, x = 10), size = 2.1, hjust = 0, fill='white')

# #save figure
# grDevices::png(file = here::here("figures", "9_dimensionality_trait_omission.png"),width = 6250,height=2500,res=500)
# print(p2)
# dev.off()

####
## Not used: Singleton scatterplot ----
####

# prepare data
load(file = here::here("outputs/9_dimensionality", "res_for_model.RData"))

filenames <- list.files(path = here::here("outputs/9_dimensionality"), pattern = "_res.rds$",
                        full.names = FALSE)

files <- list.files(path = here::here("outputs/9_dimensionality"), pattern = "_single.rds$",
                    full.names = TRUE)
list_res_single <- lapply(files, function(x) readRDS(x))

files <- list.files(path = here::here("outputs/9_dimensionality"), pattern = "_pcoa.rds$",
                    full.names = TRUE)
list_res_pcoa  <- lapply(files, function(x) readRDS(x))

# build data frame for singleton
res_for_graph_single <- na.omit(data.frame(do.call(rbind, lapply(1:length(list_res_single),
                                                                 function(i) {
                                                                   res_single <- data.frame(list_res_single[[i]]$cluster_core)
                                                                   
                                                                   # cluster_core = 1 === singleton
                                                                   res_single[res_single[ , 1] >  0, ] <- 2
                                                                   res_single[res_single[ , 1] == 0, ] <- 1
                                                                   res_single[res_single[ , 1] >  1, ] <- 0
                                                                   
                                                                   res_single$taxa <- gsub("_res\\.rds", "", filenames)[i]
                                                                   
                                                                   res_pcoa <- data.frame(list_res_pcoa[[i]]$vectors)[ , c(1:3)]
                                                                   
                                                                   res_cluster <- list_res_single[[i]]$cluster_core
                                                                   res <- cbind(res_single, res_pcoa, res_cluster)
                                                                   colnames(res) <- c("Single","taxa", "Pcoa1", "Pcoa2", "Pcoa3","cluster_ID")
                                                                   
                                                                   return(res)
                                                                 }))))


# number of species and number of traits columns
res_for_graph_single$SP    <- NA 
res_for_graph_single$trait <- NA 

for (i in 1:nrow(res_for_graph_single)) { 
  
  res_for_graph_single$SP[i]    <- res_for_model[rownames(res_for_model) %in% res_for_graph_single$taxa[i], ]$S
  res_for_graph_single$trait[i] <- res_for_model[rownames(res_for_model) %in% res_for_graph_single$taxa[i], ]$Nb_trait
}

# make dataset column and order by 
res_for_graph_single$taxa <- factor(res_for_graph_single$taxa, 
                                    levels  = unique(res_for_graph_single$taxa[order(res_for_graph_single$SP)]), 
                                    ordered = TRUE)
res_for_graph_single <-res_for_graph_single[order(res_for_graph_single$taxa, decreasing = FALSE),]

# modify the PCoA coordinates by jittering
res_for_graph_single$Pcoa1  <- jitter(as.numeric(as.character(res_for_graph_single$Pcoa1)), factor = 50)
res_for_graph_single$Pcoa2  <- jitter(as.numeric(as.character(res_for_graph_single$Pcoa2)), factor = 50)
res_for_graph_single$Pcoa3  <- as.numeric(as.character(res_for_graph_single$Pcoa3))

# make singleton status a factor
res_for_graph_single$Single <- as.factor(as.character(res_for_graph_single$Single))

# NOTE: not sure why this is relevant, seems to make any cluster other than 1 be 0
res_for_graph_single$cluster_ID[res_for_graph_single$cluster_ID != 1] <- 0

# order by species alphabetically
res_for_graph_single<-res_for_graph_single[order(rownames(res_for_graph_single)),]

# cluster_core = 1 === singleton
p3 <- ggplot(res_for_graph_single, aes(x = Pcoa1, y = Pcoa2)) + 
  geom_point(aes(alpha = Single, shape = Single), size = 2) + 
  scale_shape_manual(values = c(4, 16)) +
  scale_alpha_manual(values = c(0.5, 0.8)) +
  ggalt::geom_encircle(s_shape = 1, expand = 0,size = 3, alpha = 0.7, show.legend = FALSE) +
  theme_bw() +
  labs(x = "PCoA axis 1") +
  labs(y = "PCoA axis 2") +
  facet_wrap(~ taxa,ncol = 6, scales = "free") +
  harrypotter::scale_colour_hp_d(option = "LunaLovegood", direction = 1) +
  theme(strip.background = element_blank(),
        strip.text.x     = element_blank(),
        strip.text.y     = element_blank(),
        panel.grid.major = element_blank(), 
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        #legend.position  = "none",
        axis.text.x      = element_blank(),   
        axis.text.y      = element_blank(),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        axis.ticks       = element_blank()) + 
  geom_label(data = res_for_graph_single, 
             aes(label = paste0("data set = ", taxa),
                 y = 0.5, x = -0.5), size = 2.1, hjust = 0)

# save plot
# grDevices::png(file = here::here("figures", "9_dimensionality_singleton.png"),width = 8000,height=2000,res=400) #SAVE A4
# print(p3)
# dev.off()
