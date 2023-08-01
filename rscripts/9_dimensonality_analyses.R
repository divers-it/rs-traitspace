rm(list=ls())

## Load Project Addins (R Functions and Packages) ----
devtools::load_all()

#not loading ggplot?
library(ggplot2)

# Variable to iterate over
percent_list <- seq(0.1, 0.8, by = 0.1)

###
# Original DiveRS dataset
###

# Read Dataset ----
dataset <- readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

# Run dimensionality nalyses ----
run_analysis(dataset, name = "DiveRS_2023")
rm(list = "dataset")

###
# DiveRS one-hot dataset
###

# Read Dataset ----
dataset <- readRDS(file = here::here("outputs/one_hot_6_df_filt_trans.rds"))

# Run Analysis ----
run_analysis(dataset, name = "DiveRS_one_hot_2023")
rm(list = "dataset")

###
# DiveRS imputed dataset with no NA
###

# Read Dataset ----
dataset <- read.csv(file = here::here("outputs/imputed_with_phylo.csv"),row.names = 1,stringsAsFactors = TRUE)

# Run Analysis ----
run_analysis(dataset, name = "DiveRS_imputed_2023")
rm(list = "dataset")

###
# Species from DiveRS data set also in Diaz data set
###

# Read Dataset ----
dataset <- readRDS(file = here::here("outputs/8.4_df_filt_trans_shared.rds"))

# Run Analysis ----
run_analysis(dataset, name = "DiveRS_shared_2023")
rm(list = "dataset")

###
# Species from Diaz data set that are in DiveRS data set
###

#load formatted DiveRS data
df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

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

rownames(diaz_c)<-diaz_c$Species_name_standardized_against_TPL

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

#select only DiveRS species
diaz_pcf_od<-diaz_pcf[rownames(diaz_pcf)%in%rownames(df),]

#remove 0s by adding 1
diaz_pcf_od$Diaspore_mass_mg<-diaz_pcf_od$Diaspore_mass_mg + 1

#log transform and scale
diaz_pcf_od<-log(diaz_pcf_od)
diaz_pcf_od<-scale(diaz_pcf_od, center = T, scale = T)

#rename for analysis
dataset<-data.frame(diaz_pcf_od)

# Run Analysis ----
run_analysis(dataset, name = "Diaz_shared_2022")
rm(list = "dataset")

## Import Results ----

# This will import all of the results present in the outputs folder (from each different dimensionality run)

files <- list.files(path = here::here("outputs"), pattern = "_res.rds$",
                    full.names = TRUE)
list_res <- lapply(files, function(x) readRDS(x))

files <- list.files(path = here::here("outputs"), pattern = "_cor.rds$",
                    full.names = TRUE)
list_res_cor <- lapply(files, function(x) readRDS(x))

files <- list.files(path = here::here("outputs"), pattern = "_single.rds$",
                    full.names = TRUE)
list_res_sigle <- lapply(files, function(x) readRDS(x))

files <- list.files(path = here::here("outputs"), pattern = "_pcoa.rds$",
                    full.names = TRUE)
list_res_pcoa <- lapply(files, function(x) readRDS(x))

filenames <- list.files(path = here::here("outputs"), pattern = "_pcoa.rds$",
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

#make out of 1
res_for_model$Percentage_lostAUC_depleted0.5  <- res_for_model$Percentage_lostAUC_depleted0.5  / 100
res_for_model$Percentage_lostAUC_depleted0.20 <- res_for_model$Percentage_lostAUC_depleted0.20 / 100

#rename columns
colnames(res_for_model)[14] <- "rowAUClostwhen50percTraitdepleted"
colnames(res_for_model)[15] <- "rowAUClostwhen20percTraitdepleted"

#save results
save(res_for_model, file = here::here("outputs", "res_for_model.RData"))

#' FIGURE 2 Quality of species trait spaces (AUC)

## Prepare Data ----

load(file = here::here("outputs", "res_for_model.RData"))

#get files and names
filenames <- list.files(path = here::here("outputs"), pattern = "_res.rds$",
                        full.names = FALSE)

files <- list.files(path = here::here("outputs"), pattern = "_res.rds$",
                    full.names = TRUE)

#read in results of all analyses
list_res <- lapply(files, function(x) readRDS(x))

#get data frame of trait space quality metrics
res_for_graph_dim <- data.frame(do.call(rbind, lapply(1:length(list_res),
                                                      function(i) {
                                                        res <- data.frame(list_res[[i]][[3]])
                                                        res$taxa <- gsub("_res\\.rds", "", filenames)[i]
                                                        return(res)
                                                      })))


#make empty columns for species and number of traits
res_for_graph_dim$SP    <- NA 
res_for_graph_dim$trait <- NA 

#fill columns
for (i in 1:nrow(res_for_graph_dim)){ 
  
  res_for_graph_dim$SP[i]    <- res_for_model[rownames(res_for_model) %in% res_for_graph_dim$taxa[i], ]$S
  res_for_graph_dim$trait[i] <- res_for_model[rownames(res_for_model) %in% res_for_graph_dim$taxa[i], ]$Nb_trait
  res_for_graph_dim$elbow[i] <- res_for_model[rownames(res_for_model) %in% res_for_graph_dim$taxa[i], ]$Nb_dim_AUC_elbow
}


#make columns for elbow method and associated AUC values
res_for_graph_dim$selec_elbow_graph <- NA
res_for_graph_dim$AUCwhenelbow      <- NA

#fill columns
for (i in 1:nrow(res_for_graph_dim)) { 
  
  if (res_for_graph_dim$dim[i] == res_for_graph_dim$elbow[i]) {
    
    res_for_graph_dim$selec_elbow_graph[i] <- res_for_graph_dim$dim[i]
    res_for_graph_dim$AUCwhenelbow[i]      <- res_for_graph_dim$AUC[i]
  }
}


## Sort plots by number of species or traits ----

res_for_graph_dim$taxa <- factor(res_for_graph_dim$taxa, 
                                 levels = unique(res_for_graph_dim$taxa[order(res_for_graph_dim$SP)]), 
                                 ordered = TRUE)

res_for_graph_dim <- res_for_graph_dim[order(res_for_graph_dim$taxa, decreasing = FALSE), ]


## Plot ----

p <- ggplot(res_for_graph_dim, aes(x = dim, y = AUC, colour = taxa)) + 
  stat_summary(fun = "mean", geom = "line", size = 1, alpha = 0.4) +
  stat_summary(fun = "mean", size = 0.88) +
  
  labs(x = "Number of dimensions(PCoA axes)") + 
  labs(y = "Quality of species trait space (AUC)") +
  facet_wrap(~ taxa,ncol = 6) + theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        panel.background = element_blank(),
        legend.position = "none") +
  harrypotter::scale_colour_hp_d(option = "LunaLovegood") +
  
  geom_segment(data = plyr::ddply(res_for_graph_dim, "taxa", dplyr::summarize, wavg = AUC), 
               aes(x = res_for_graph_dim$selec_elbow_graph, 
                   xend = res_for_graph_dim$selec_elbow_graph,
                   y = 0 , yend = res_for_graph_dim$AUCwhenelbow),
               color = "black", linetype = "dotted", size = 1) +
  
  geom_segment(data = plyr::ddply(res_for_graph_dim, "taxa", dplyr::summarize, wavg = AUC), 
               aes(y = res_for_graph_dim$AUCwhenelbow,
                   yend = res_for_graph_dim$AUCwhenelbow ,
                   x = 0 , xend = res_for_graph_dim$selec_elbow_graph),
               color = "black", linetype = "dotted", size = 1) +
  
  geom_point(data = plyr::ddply(res_for_graph_dim, "taxa", dplyr::summarize, wavg = AUC), 
             aes(y = res_for_graph_dim$AUCwhenelbow,
                 x = res_for_graph_dim$selec_elbow_graph),
             color = "black", size = 4, shape = 19) +
  
  geom_label(data = res_for_graph_dim, 
             aes(label = paste0("Elbow-AUC = ", elbow, "\n", 
                                "#S = ", SP , "\n",
                                "#T = ", trait),
                 y = 0.5, x = 10), size = 2.1, hjust = 0) +
  
  geom_label(data = res_for_graph_dim, 
             aes(label = paste0("data set = ", taxa),
                 y = 0.95, x = 0), size = 2.1, hjust = 0) +
  
  scale_y_continuous(breaks = seq(0.1, 1, 0.2))


#save plot
grDevices::png(file = here::here("figures", "9_dimensionality_no_axes.png"),width = 6250,height=2500,res=500)
print(p)
dev.off()

#' FIGURE 5 Trait omission vs trait space quality

## Prepare Data ----

load(file = here::here("outputs", "res_for_model.RData"))

#get files and names
filenames <- list.files(path = here::here("outputs"), pattern = "_miss.rds$",
                        full.names = FALSE)

files <- list.files(path = here::here("outputs"), pattern = "_miss.rds$",
                    full.names = TRUE)

#read in files
list_res <- lapply(files, function(x) readRDS(x))

#make missingness data frame
res_for_graph_miss <- na.omit(data.frame(do.call(rbind, lapply(1:length(list_res),
                                                               function(i) {
                                                                 res <- data.frame(list_res[[i]])
                                                                 res$taxa <- gsub("_miss\\.rds", "", filenames)[i]
                                                                 return(res)
                                                               }))))

#number of species and number of traits columns.
res_for_graph_miss$SP    <- NA 
res_for_graph_miss$trait <- NA 

for (i in 1:nrow(res_for_graph_miss)){ 
  
  res_for_graph_miss$SP[i]    <- res_for_model[rownames(res_for_model) %in% res_for_graph_miss$taxa[i], ]$S
  res_for_graph_miss$trait[i] <- res_for_model[rownames(res_for_model) %in% res_for_graph_miss$taxa[i], ]$Nb_trait
}

#dataset labels
res_for_graph_miss$taxa <- factor(res_for_graph_miss$taxa, 
                                  levels  = unique(res_for_graph_miss$taxa[order(res_for_graph_miss$SP)]), 
                                  ordered = TRUE)

#reorder by data set
res_for_graph_miss <- res_for_graph_miss[order(res_for_graph_miss$taxa, decreasing = FALSE), ]



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

#save figure
grDevices::png(file = here::here("figures", "9_dimensionality_trait_omission.png"),width = 6250,height=2500,res=500)
print(p2)
dev.off()

#' FIGURE 9 Singleton scatterplot

## Prepare Data ----

load(file = here::here("outputs", "res_for_model.RData"))

filenames <- list.files(path = here::here("outputs"), pattern = "_res.rds$",
                        full.names = FALSE)

files <- list.files(path = here::here("outputs"), pattern = "_single.rds$",
                    full.names = TRUE)
list_res_single <- lapply(files, function(x) readRDS(x))

files <- list.files(path = here::here("outputs"), pattern = "_pcoa.rds$",
                    full.names = TRUE)
list_res_pcoa  <- lapply(files, function(x) readRDS(x))

#build data frame for singleton
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


#number of species and number of traits columns
res_for_graph_single$SP    <- NA 
res_for_graph_single$trait <- NA 

for (i in 1:nrow(res_for_graph_single)) { 
  
  res_for_graph_single$SP[i]    <- res_for_model[rownames(res_for_model) %in% res_for_graph_single$taxa[i], ]$S
  res_for_graph_single$trait[i] <- res_for_model[rownames(res_for_model) %in% res_for_graph_single$taxa[i], ]$Nb_trait
}

#make dataset column and order by 
res_for_graph_single$taxa <- factor(res_for_graph_single$taxa, 
                                    levels  = unique(res_for_graph_single$taxa[order(res_for_graph_single$SP)]), 
                                    ordered = TRUE)
res_for_graph_single <-res_for_graph_single[order(res_for_graph_single$taxa, decreasing = FALSE),]

#modify the PCoA coordinates by jittering
res_for_graph_single$Pcoa1  <- jitter(as.numeric(as.character(res_for_graph_single$Pcoa1)), factor = 50)
res_for_graph_single$Pcoa2  <- jitter(as.numeric(as.character(res_for_graph_single$Pcoa2)), factor = 50)
res_for_graph_single$Pcoa3  <- as.numeric(as.character(res_for_graph_single$Pcoa3))

#make singleton status a factor
res_for_graph_single$Single <- as.factor(as.character(res_for_graph_single$Single))

#NOTE: not sure why this is relevant, seems to make any cluster other than 1 be 0
res_for_graph_single$cluster_ID[res_for_graph_single$cluster_ID != 1] <- 0

#order by species alphabetically
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

grDevices::png(file = here::here("figures", "9_dimensionality_singleton.png"),width = 8000,height=2000,res=400) #SAVE A4
print(p3)
dev.off()
