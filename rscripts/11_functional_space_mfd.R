rm(list=ls())

# FROM:
# https://frbcesab.github.io/workshop-free/practice.html

# load packages
library(tidyr)
library(dplyr)
library(mFD)
library(ggplot2)
library(funrar)

# load formatted data
df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

# empty trait code vector
trait_code<-vector()

# read in columns and determine coding (NO ORDINAL YET)
for(i in 1:length(colnames(df))){
  if(is.numeric(df[,i])){
    trait_code[i]<-"Q"
  }
  
  if(is.factor(df[,i])){
    trait_code[i]<-"N"
  }
  
}

# make data frame with trait types
trait_code_df<-data.frame(colnames(df),trait_code)
colnames(trait_code_df)<-c("trait_name","trait_type")

# calculate functional distances using gower
sp_dist <- mFD::funct.dist(
  sp_tr         = df,
  tr_cat        = trait_code_df,
  metric        = "gower",
# scale_euclid  = "scale_center", # already done
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = F)

# do PCoA and evaluate quality of space as axes are added
fspaces_quality <- mFD::quality.fspaces(
  sp_dist             = sp_dist,
  maxdim_pcoa         = 15,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "ward.D")

## Quality metrics of functional spaces ----

# The space with the best quality has the lowest quality metric.
round(fspaces_quality$"quality_fspaces", 4)

# With the mFD package, it is possible to illustrate the quality of PCoA-based multidimensional spaces according
# to deviation between trait-based distances and distances in the functional space

# This function generates a figure with three panels (in rows) for each selected functional space (in columns). 
# Each column represents a functional space, the value of the quality metric is written on the top of each column. 
# The x-axis of all panels represents trait-based distances. The y-axis is different for each row:

# on the first (top) row, the y-axis represents species functional distances in the multidimensional space. 
# Thus, the closer species are to the 1:1 line, the better distances in the functional space fit trait-based ones

# on the second row, the y-axis shows the raw deviation of species distances in the functional space compared to 
# trait-based distances. Thus, the raw deviation reflects the distance to the horizontal line.

# on the third row (bottom), the y-axis shows the absolute or squared deviation of the (“scaled”) distance in 
# the functional space. It is the deviation that is taken into account for computing the quality metric.
mFD::quality.fspaces.plot(
  fspaces_quality            = fspaces_quality,
  quality_metric             = "mad",
  fspaces_plot               = c("pcoa_2d", "pcoa_3d", 
                                 "pcoa_4d", "pcoa_5d", "pcoa_6d", "pcoa_7d"),
  name_file                  = NULL,
  range_dist                 = NULL,
  range_dev                  = NULL,
  range_qdev                 = NULL,
  gradient_deviation         = c(neg = "darkblue", nul = "grey80", pos = "darkred"),
  gradient_deviation_quality = c(low = "yellow", high = "red"),
  x_lab                      = "Trait-based distance")

# ggsave("figures/11_scatterplots_mfd_fspace_quality.png",width=30,height=15,units = "cm")

## Testing the correlation between functional axes and traits ----

# mFD allows to test for correlations between traits and functional axes and then illustrate possible correlations 
# continuous traits = linear model is computed and r2 and associated p-value are returned
# non-continuous traits = Kruskal-Wallis test is computed and eta2 statistic is returned

# a matrix of species coordinates taken from the output
sp_faxes_coord <- fspaces_quality$"details_fspaces"$"sp_pc_coord"

# correlation of continuous traits
df_faxes <- mFD::traits.faxes.cor(
  sp_tr          = df[,c(1:7)], 
  sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4")],
  stop_if_NA = F,
  plot = TRUE)

# Print traits with significant effect
df_faxes$"tr_faxes_stat"[which(df_faxes$"tr_faxes_stat"$"p.value" < 0.05), ]

# plot
df_faxes$"tr_faxes_plot"
# ggsave("figures/11_scatterplots_mfd_quant_traits_vs_axes.png",width=30,height=20,units = "cm")

# correlation of discrete traits 8-14
df_faxes <- mFD::traits.faxes.cor(
  sp_tr          = df[,c(8:14)], 
  sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4")],
  stop_if_NA = F,
  plot = TRUE)

# Print traits with significant effect
df_faxes$"tr_faxes_stat"[which(df_faxes$"tr_faxes_stat"$"p.value" < 0.05), ]

# plot
df_faxes$"tr_faxes_plot"
# ggsave("figures/11_scatterplots_mfd_qual_traits_vs_axes_1.png",width=30,height=20,units = "cm")

# correlation of qualitative traits 15-21
df_faxes <- mFD::traits.faxes.cor(
  sp_tr          = df[,c(15:21)], 
  sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4")],
  stop_if_NA = F,
  plot = TRUE)

# Print traits with significant effect
df_faxes$"tr_faxes_stat"[which(df_faxes$"tr_faxes_stat"$"p.value" < 0.05), ]

# plot
df_faxes$"tr_faxes_plot"
# ggsave("figures/11_scatterplots_mfd_qual_traits_vs_axes_2.png",width=30,height=20,units = "cm")

## Plotting the selected functional space and position of species ----

# get coordinates
sp_faxes_coord <- fspaces_quality$"details_fspaces"$"sp_pc_coord"

# plot
big_plot <- mFD::funct.space.plot(
  sp_faxes_coord  = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4")],
  faxes           = c("PC1", "PC2", "PC3", "PC4"),
  name_file       = NULL,
  faxes_nm        = NULL,
  range_faxes     = c(NA, NA),
  plot_ch         = TRUE,
  plot_vertices   = TRUE,
  plot_sp_nm      = NULL,
  check_input     = TRUE)

big_plot$"patchwork"
# ggsave("figures/11_scatterplots_mfd_functional_space.png",width=30,height=30,units = "cm")


####
## Definitions of FD indices ----
####

# FDis Functional Dispersion: the biomass weighted deviation of species traits values from the center of the functional space filled by the # assemblage i.e. the biomass-weighted mean distance to the biomass-weighted mean trait values of the assemblage.
# 
# FRic Functional Richness: the proportion of functional space filled by species of the studied assemblage, i.e. the volume inside the convex# -hull shaping species. To compute FRic the number of species must be at least higher than the number of functional axis + 1.
# 
# FDiv Functional Divergence: the proportion of the biomass supported by the species with the most extreme functional traits i.e. the ones # located close to the edge of the convex-hull filled by the assemblage.
# 
# FEve Functional Evenness: the regularity of biomass distribution in the functional space using the Minimum Spanning Tree linking all species # present in the assemblage.
# 
# FSpe Functional Specialization: the biomass weighted mean distance to the mean position of species from the global pool (present in all # assemblages).
# 
# FMPD Functional Mean Pairwise Distance: the mean weighted distance between all species pairs.
# 
# FNND Functional Mean Nearest Neighbour Distance: the weighted distance to the nearest neighbor within the assemblage.
# 
# FIde Functional Identity: the mean traits values for the assemblage. FIde is always computed when FDis is computed.
# 
# FOri Functional Originality: the weighted mean distance to the nearest species from the global species pool.

####
## FD indices for clusters ----
####

# read in clustering
clust <- readRDS("outputs/10_clust_num_k_2_7_pam.rds")

# read in robust groups
robust <- readRDS("outputs/10_robust_vect_pam_full.rds")

# check alignment
table(rownames(clust) == names(robust))

# make non-robust NA
clust[is.na(robust),]<-NA

# make df with clustering info
clust_df<-as.data.frame(cbind(rownames(clust),clust$`3clusters`))
colnames(clust_df)<-c("species","cluster")

# recode df into one-hot with species as columns
clust_recode <- clust_df %>% mutate(value = 1)  %>% spread(species, value,  fill = 0) 

#remove NA row
clust_recode <- clust_recode[c(1:3),]

# remove cluster label column and add as name
cw<-clust_recode[,c(2:length(colnames(clust_recode)))]
rownames(cw)<-clust_recode$cluster

# check names
table(colnames(cw)==rownames(sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4")]))

# compute functional indices
alpha_fd_indices <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = as.matrix(sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4")]),
  asb_sp_w         = as.matrix(cw),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

# output indices
fd_ind_values <- alpha_fd_indices$"functional_diversity_indices"
fd_ind_values

fd_ind_table<-round(fd_ind_values[,c(1:9)],3)
rownames(fd_ind_table)<-c("Cluster 1",
                          "Cluster 2",
                          "Cluster 3")

colnames(fd_ind_table) <- c("Species richness",
                            "FDis",
                            "FMPD",
                            "FNND",
                            "FEve",
                            "FRic",
                            "FDiv",
                            "FOri",
                            "FSpe")

write.csv(fd_ind_table,"outputs/11_mfd_indices_PAM_k3.csv")


# information such as coordinates of centroids, distances and identity of the nearest neighbour, 
# distances to the centroid, etc. The user does not have to directly use it but it will be useful 
# if FD indices are then plotted. It can be retrieved through:
details_list <- alpha_fd_indices$"details"
details_list

#### 
## FD indices sexual system ----
#### 

# get flower sex and sexual system columns
ss_df <- data.frame(rownames(df),df[,c("SexualSystem","FlowerSex")])
colnames(ss_df)<-c("species","SexualSystem","FlowerSex")

# convert these to bisexual, monoecy, dioecy
ss_df$new_ss <- rep(NA, length(ss_df[,1]))

### Recode sexual system ----
for (i in 1:length(ss_df[, 1])) {
  
  # bisexual species
  if (grepl("monomorphic", ss_df$SexualSystem[i]) &&
      grepl("bisexual", ss_df$FlowerSex[i], )) {
    
    if(is.na(ss_df$new_ss[i])){
      
      ss_df$new_ss[i] <- "bisexual"
      
    } else {
      
      ss_df$new_ss[i] <- paste(ss_df$new_ss[i],"bisexual",sep="_")
      
    }
    
  }
  
  if (grepl("monomorphic", ss_df$SexualSystem[i]) &&
      grepl("unisexual", ss_df$FlowerSex[i])) {
    
    if(is.na(ss_df$new_ss[i])){
      
      ss_df$new_ss[i] <- "monoecy"
      
    } else {
      
      ss_df$new_ss[i] <- paste(ss_df$new_ss[i],"monoecy",sep="_")
      
    }
    
  }
  
  if (grepl("dimorphic", ss_df$SexualSystem[i]) &&
      grepl("unisexual", ss_df$FlowerSex[i])) {
    
    if(is.na(ss_df$new_ss[i])){
      
      ss_df$new_ss[i] <- "dioecy"
      
    } else {
      
      ss_df$new_ss[i] <- paste(ss_df$new_ss[i],"dioecy",sep="_")
      
    }
    
  }
  
}

ss_df$new_ss

# get rid of other columns
ss_df <- ss_df[,c("species","new_ss")]
ss_df$new_ss <- as.factor(ss_df$new_ss)

# remove missing species
sfc_ss<-sp_faxes_coord[!is.na(ss_df$new_ss),]
ss_df<-ss_df[!is.na(ss_df$new_ss),]
ss_df$new_ss<-droplevels(ss_df$new_ss)

# recode df into one-hot with species as columns
ss_recode <- ss_df %>% mutate(value = 1)  %>% spread(species, value,  fill = 0 ) 

# remove NA column
ss_recode <- ss_recode[1:length(levels(ss_df$new_ss)),]

# remove cluster label column and add as name
ss<-ss_recode[,c(2:length(colnames(ss_recode)))]
rownames(ss)<-as.character(ss_recode$new_ss[1:6])

# check names
table(colnames(ss)==rownames(sfc_ss[ , c("PC1", "PC2", "PC3", "PC4")]))

# remove unwanted sexual systems
ss <- ss[c("bisexual","monoecy","dioecy"),]

# compute 5 functional indices
alpha_fd_indices <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = as.matrix(sfc_ss[ , c("PC1", "PC2", "PC3", "PC4")]),
  asb_sp_w         = as.matrix(ss),
  # ind_vect         = c("fdis", "fric", "fdiv","fspe", "fide"), # output all
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

# output indices
fd_ind_values_ss <- alpha_fd_indices$"functional_diversity_indices"
fd_ind_values_ss

# make table
fd_ind_table_ss<-round(fd_ind_values_ss[,c(1:9)],3)
rownames(fd_ind_table_ss)<-c("monocliny",
                          "monoecy",
                          "dioecy")

colnames(fd_ind_table_ss) <- c("Species richness",
                            "FDis",
                            "FMPD",
                            "FNND",
                            "FEve",
                            "FRic",
                            "FDiv",
                            "FOri",
                            "FSpe")


write.csv(fd_ind_table_ss,"outputs/11_mfd_indices_sexual_system.csv")

### Plot functional indices ----

# colours
pcols <- c(pool = "grey70", 
  asb1 = "skyblue4",
  asb2 = "orange")

# plot functional indices of monoecy and dioecy
plots_alpha_md <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices,
  plot_asb_nm              = c("monoecy", "dioecy"),
  ind_nm                   = c("fdis", "fric", "fdiv", 
                               "fspe", "fide"),
  faxes                    = c("PC1","PC2","PC3"),
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_sp                 = pcols,
  color_vert               = pcols,
  fill_sp = pcols,
  fill_vert = pcols,
  color_ch = pcols,
  fill_ch = pcols,
  alpha_ch                = c(pool = 0.1, asb1 = 0.2, 
                              asb2 = 0.2),
  plot_sp_nm               = NULL,
  save_file                = FALSE,
  check_input              = TRUE)

# colors 
pcols <- c(pool = "grey70", 
  asb1 = "skyblue4",
  asb2 = "#1F968BFF")

# plot functional indices of monoecy and bisexual
plots_alpha_mb <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices,
  plot_asb_nm              = c("monoecy", "bisexual"),
  ind_nm                   = c("fdis", "fric", "fdiv", 
                               "fspe", "fide"),
  faxes                    = c("PC1","PC2","PC3"),
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_sp                 = pcols,
  color_vert               = pcols,
  fill_sp = pcols,
  fill_vert = pcols,
  color_ch = pcols,
  fill_ch = pcols,
  alpha_ch                = c(pool = 0.1, asb1 = 0.2, 
                              asb2 = 0.2),
  plot_sp_nm               = NULL,
  save_file                = FALSE,
  check_input              = TRUE) 

# FDis Functional Dispersion: the biomass weighted deviation of species traits values from the center of the functional space filled by the # assemblage i.e. the biomass-weighted mean distance to the biomass-weighted mean trait values of the assemblage.
plots_alpha_md$"fdis"$"patchwork"
plots_alpha_mb$"fdis"$"patchwork"

# FRic representation: the colored shapes reflect the convex-hull of the studied assemblages
# and the white shape reflects the convex-hull of the global pool of species:
plots_alpha_md$"fric"$"patchwork"
plots_alpha_mb$"fric"$"patchwork"

# FDiv representation: the gravity centers of vertices (i.e. species with the most extreme functional traits) of each 
# assemblages are plotted as a square and a triangle. The two colored circles represent the mean
# distance of species to the gravity center for each assemblage. Species of each assemblage 
# have different size given their relative weight into the assemblage.
plots_alpha_md$"fdiv"$"patchwork"
plots_alpha_mb$"fdiv"$"patchwork"

# FSpe representation: colored traits represent distances of each species from a given assemblage 
# to the center of gravity of the global pool (i.e center of the functional space). the center of
# gravity is plotted with a purple diamond. Species of each assemblage have different size given
# their relative weight into the assemblage.
plots_alpha_md$"fspe"$"patchwork"
plots_alpha_mb$"fspe"$"patchwork"

# FIde representation:colored lines refer to the weighted average position of species of each assemblage
# along each axis. Species of each assemblage have different size given their relative weight
# into the assemblage.
plots_alpha_md$"fide"$"patchwork"
plots_alpha_mb$"fide"$"patchwork"

#### 
## FD indices for mating system ----
#### 

ms_df <- data.frame(rownames(df),df[,"Mating"])
colnames(ms_df)<-c("species","Mating")

# remove missing species
sfc_ms<-sp_faxes_coord[!is.na(ms_df$Mating),]
ms_df<-ms_df[!is.na(ms_df$Mating),]
ms_df$Mating<-droplevels(ms_df$Mating)

# recode df into one-hot with species as columns
ms_recode <- ms_df %>% mutate(value = 1)  %>% spread(species, value,  fill = 0 ) 

# remove NA column
ms_recode <- ms_recode[1:length(levels(ms_df$Mating)),]

# remove cluster label column and add as name
ms<-ms_recode[,c(2:length(colnames(ms_recode)))]
rownames(ms)<-as.character(ms_recode$Mating[1:3])

# check names
table(colnames(ms)==rownames(sfc_ms[ , c("PC1", "PC2", "PC3", "PC4")]))

# compute 5 functional indices
alpha_fd_indices <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = as.matrix(sfc_ms[ , c("PC1", "PC2", "PC3", "PC4")]),
  asb_sp_w         = as.matrix(ms),
  # ind_vect         = c("fdis", "fric", "fdiv","fspe", "fide"), # output all
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

# output indices
fd_ind_values_ms <- alpha_fd_indices$"functional_diversity_indices"
fd_ind_values_ms

# make table
fd_ind_table_ms<-round(fd_ind_values_ms[,c(1:9)],3)
rownames(fd_ind_table_ms)<-c("mixed",
                             "outcrossing",
                             "selfing")

colnames(fd_ind_table_ms) <- c("Species richness",
                               "FDis",
                               "FMPD",
                               "FNND",
                               "FEve",
                               "FRic",
                               "FDiv",
                               "FOri",
                               "FSpe")



write.csv(fd_ind_table_ms,"outputs/11_mfd_indices_mating_system.csv")

### Plot functional indices ----

# colours
pcols <- c(pool = "grey70", 
           asb1 = "darkblue",
           asb2 = "red")

# plot functional indices of monoecy and dioecy
plots_alpha_md <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices,
  plot_asb_nm              = c("selfing", "outcrossing"),
  ind_nm                   = c("fdis", "fric", "fdiv", 
                               "fspe", "fide"),
  faxes                    = c("PC1","PC2","PC3"),
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_sp                 = pcols,
  color_vert               = pcols,
  fill_sp = pcols,
  fill_vert = pcols,
  color_ch = pcols,
  fill_ch = pcols,
  alpha_ch                = c(pool = 0.1, asb1 = 0.2, 
                              asb2 = 0.2),
  plot_sp_nm               = NULL,
  save_file                = FALSE,
  check_input              = TRUE)

# colors 
pcols <- c(pool = "grey70", 
           asb1 = "purple",
           asb2 = "red")

# plot functional indices of monoecy and bisexual
plots_alpha_mb <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices,
  plot_asb_nm              = c("mixed", "outcrossing"),
  ind_nm                   = c("fdis", "fric", "fdiv", 
                               "fspe", "fide"),
  faxes                    = c("PC1","PC2","PC3"),
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  color_sp                 = pcols,
  color_vert               = pcols,
  fill_sp = pcols,
  fill_vert = pcols,
  color_ch = pcols,
  fill_ch = pcols,
  alpha_ch                = c(pool = 0.1, asb1 = 0.2, 
                              asb2 = 0.2),
  plot_sp_nm               = NULL,
  save_file                = FALSE,
  check_input              = TRUE) 

# FDis Functional Dispersion: the biomass weighted deviation of species traits values from the center of the functional space filled by the # assemblage i.e. the biomass-weighted mean distance to the biomass-weighted mean trait values of the assemblage.
plots_alpha_md$"fdis"$"patchwork"
plots_alpha_mb$"fdis"$"patchwork"

# FRic representation: the colored shapes reflect the convex-hull of the studied assemblages
# and the white shape reflects the convex-hull of the global pool of species:
plots_alpha_md$"fric"$"patchwork"
plots_alpha_mb$"fric"$"patchwork"

# FDiv representation: the gravity centers of vertices (i.e. species with the most extreme functional traits) of each 
# assemblages are plotted as a square and a triangle. The two colored circles represent the mean
# distance of species to the gravity center for each assemblage. Species of each assemblage 
# have different size given their relative weight into the assemblage.
plots_alpha_md$"fdiv"$"patchwork"
plots_alpha_mb$"fdiv"$"patchwork"

# FSpe representation: colored traits represent distances of each species from a given assemblage 
# to the center of gravity of the global pool (i.e center of the functional space). the center of
# gravity is plotted with a purple diamond. Species of each assemblage have different size given
# their relative weight into the assemblage.
plots_alpha_md$"fspe"$"patchwork"
plots_alpha_mb$"fspe"$"patchwork"

# FIde representation:colored lines refer to the weighted average position of species of each assemblage
# along each axis. Species of each assemblage have different size given their relative weight
# into the assemblage.
plots_alpha_md$"fide"$"patchwork"
plots_alpha_mb$"fide"$"patchwork"


####
## Functional distinctiveness ----
####

# calculate distinctiveness
sp_di <- distinctiveness_global(sp_dist, di_name = "distinctiveness")

# We get one value of distinctiveness per species. 
# It only considers the functional dissimilarity of all species in the dissimilarity matrix without
# considering their spatial distributions.
summary(sp_di)
sp_di 

# check matching
table(sp_di$species == clust_df$species)

quantile(sp_di$distinctiveness, probs = seq(0, 1, by = 0.1))
subset(sp_di, distinctiveness >= 0.439)

# add clustering
sp_di$cluster <- clust_df$cluster 

# read in robust groups
pam <- readRDS("outputs/10_robust_vect_pam_full.rds")

# check matching
table(sp_di$species == names(pam))

# add robust groups
sp_di$robust <- paste("robust",pam,sep="_")

### Plot distinctiveness per cluster ----

sp_di %>%
  ggplot(aes(x=cluster, y=distinctiveness, fill=cluster)) +
  geom_boxplot() +
  viridis::scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  hrbrthemes::theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A boxplot with jitter") +
  xlab("")

# Most basic violin chart
ggplot(sp_di, aes(x=cluster, y=distinctiveness, fill=cluster)) + # fill=name allow to automatically dedicate a color for each group
  geom_violin()

### Plot distinctiveness per robust group ----

# robust groups separate

sp_di %>%
  ggplot(aes(x=robust, y=distinctiveness, fill=cluster)) +
  geom_boxplot() +
  viridis::scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(aes(fill = cluster), shape = 21,size=1, alpha=0.6, width=0.1) +
  hrbrthemes::theme_ipsum() +
  theme(
    #legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("")

# robust groups grouped by

sp_di %>%
  ggplot(aes(x=cluster, y=distinctiveness, fill=robust)) +
  geom_boxplot() +
  viridis::scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  hrbrthemes::theme_ipsum() +
  theme(
    #legend.position="none",
    plot.title = element_text(size=11)
  ) +
  xlab("")


# Most basic violin chart
ggplot(sp_di, aes(x=robust, y=distinctiveness, fill=robust)) + # fill=name allow to automatically dedicate a color for each group
  geom_violin()

# For the choice or dissimilarity matrix we can use the raw dissimilarity matrix computed directly on raw 
# traits values among species, as we did here. Another option would be to compute a new functional 
# dissimilarity matrix based on the selected functional axes. One advantage of the latter is that it 
# already takes into account the correlation between traits.

# Let’s recompute regional functional distinctiveness based on the four selected functional axes. 
# Because the space comes from a PCA, we can directly use euclidean distance.
new_dissim <- dist(sp_faxes_coord[, c("PC1", "PC2", "PC3", "PC4")])
sp_di_alt <- distinctiveness_global(new_dissim, di_name = "alt_di")

# We can now compare both distinctiveness values.
sp_all_di <- merge(sp_di, sp_di_alt, by = "species")

plot(sp_all_di$distinctiveness, sp_all_di$alt_di)
cor.test(sp_all_di$distinctiveness, sp_all_di$alt_di)

# Both seems very correlated, so in our case using either one should be fine. However, it can be better 
# to use dissimilarity based on a reduced number of well-defined axes because: (1) there are more 
# interpretable thanks to the multivariate analysis, (2) the first ones contain the most information, 
# (3) they explicitly take into account potentially strong correlations between provided traits. 
# We’ll stick here with raw dissimilarity for the sake of simplicity.

####
## Regional scale distinctiveness / uniqueness ----
####

### Clusters ----

# make empty data frame
df_w<-data.frame(matrix(nrow=3, ncol=length(df$Maximumverticalheight)))
colnames(df_w)<-rownames(df)

# put 0/1 rows denoting cluster membership
df_w[1,]<-as.numeric(colnames(df_w)%in%rownames(clust[clust$`3clusters`=="k_3_cluster_1",]))
df_w[2,]<-as.numeric(colnames(df_w)%in%rownames(clust[clust$`3clusters`=="k_3_cluster_2",]))
df_w[3,]<-as.numeric(colnames(df_w)%in%rownames(clust[clust$`3clusters`=="k_3_cluster_3",]))

### Robust groups ----

# make empty data frame
df_w<-data.frame(matrix(nrow=6, ncol=length(df$Maximumverticalheight)))
colnames(df_w)<-rownames(df)

# put 0/1 rows denoting robust group membership
df_w[1,]<-as.numeric(colnames(df_w)%in%names(pam[pam==1]))
df_w[2,]<-as.numeric(colnames(df_w)%in%names(pam[pam==2]))
df_w[3,]<-as.numeric(colnames(df_w)%in%names(pam[pam==3]))
df_w[4,]<-as.numeric(colnames(df_w)%in%names(pam[pam==4]))
df_w[5,]<-as.numeric(colnames(df_w)%in%names(pam[pam==1]))
df_w[6,]<-as.numeric(colnames(df_w)%in%names(pam[pam==1]))

# To compute uniqueness at regional scale we also need the regional level functional dissimilarity matrix 
# with the uniqueness() function, and the site-species matrix:
sp_ui <- uniqueness(
  pres_matrix = as.matrix(df_w),
  as.matrix(sp_dist) # Uncomment for: Gower's distances
  # as.matrix(new_dissim) # Uncomment for: PCoA distances
)

head(sp_ui)
quantile(sp_ui$Ui, probs = seq(0, 1, by = 0.1))

# the most isolated species in the functional space. Meaning that they have the most distant nearest neighbors.
subset(sp_ui, Ui >= 0.21)

### Plot functional distinctiveness / uniqueness ----
# Color species in function of their functional originality.

# Make a summary data.frame
sp_coord_di_ui <- as.data.frame(sp_faxes_coord[, 1:2])
sp_coord_di_ui$species <- rownames(sp_coord_di_ui)
rownames(sp_coord_di_ui) <- NULL
sp_coord_di_ui <- sp_coord_di_ui[, c(3, 1, 2)]
sp_coord_di_ui <- merge(sp_coord_di_ui, sp_di, by = "species")
sp_coord_di_ui <- merge(sp_coord_di_ui, sp_ui, by = "species")

plot_reg_distinctiveness <- ggplot(sp_coord_di_ui, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(aes(color = distinctiveness)) +
  # ggrepel::geom_text_repel(aes(label = species)) +
  scale_color_viridis_c("Functional\nDistinctiveness") +
  theme_bw()

plot_reg_distinctiveness
# ggsave("figures/11_scatterplot_mfd_distinctiveness.png", width= 8, height = 8)

plot_reg_uniqueness <- ggplot(sp_coord_di_ui, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(aes(color = Ui)) +
  # ggrepel::geom_text_repel(aes(label = species)) +
  scale_color_viridis_c("Functional\nUniqueness") +
  theme_bw()

plot_reg_uniqueness
# ggsave("figures/11_scatterplot_mfd_uniqueness.png", width= 8, height = 8)

####
## FD indices table
####

# combine tables
fd_table <- rbind(fd_ind_table, fd_ind_table_ss, fd_ind_table_ms)

# scale variables
for(i in 1:length(colnames(fd_table))){
  
  fd_table[,i] <- scale(fd_table[,i])
  
}

# melt table
fd_table$group <- rownames(fd_table)
fd_table <- reshape2::melt(fd_table, id='group')

# reorder is close to order, but is made to change the order of the factor levels.
fd_table$group <- factor(fd_table$group, levels = c("Cluster 1",
                                          "Cluster 2",
                                          "Cluster 3",
                                          "monocliny",
                                          "monoecy",
                                          "dioecy",
                                          "outcrossing",
                                          "mixed",
                                          "selfing"))


# plot heatmap
ggplot(fd_table, aes(group, variable, fill= value)) + 
  scale_fill_gradient2(low = "#075AFF",
                       mid = "#FFFFCC",
                       high = "#FF0000") +
  theme_minimal() +
  geom_tile()



