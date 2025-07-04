rm(list=ls())

library(ape)
library(cluster)
library(umap)
library(ggplot2)
library(ggpubr)

# This script re-generate initial pcoa plots from gower's distance
# Then it compare them to the same analyses with imputed data and recoding in multistate characters used for phylogenetic analysis
# In addition, the results are compared with non-linear dimensional scaling

# Overall, the differences are rather limited, at least for the main patterns
# The imputed and recoded dataset seems to amplify the differences and the clustering

# Loading and formating datasets
df_final <- readRDS("outputs/6_df_filt_trans.rds")
df_onehot <- readRDS("outputs/one_hot_6_df_filt_trans.rds")
df_imputed <- read.csv("outputs/one_hot_imputed_with_phylo.csv")

# Re-ordering some levels for plotting
df_final$DispersalMode <- factor(
  df_final$DispersalMode,
  levels = c("autonomous",
             "abiotic_autonomous",
             "abiotic",
             "autonomous_biotic",
             "abiotic_autonomous_biotic",
             "abiotic_biotic",
             "biotic")
)
df_final$Pollination <- factor(
  df_final$Pollination,
  levels = c("autonomous",
             "abiotic_autonomous",
             "abiotic",
             "autonomous_biotic",
             "abiotic_biotic",
             "biotic")
)
df_final$Mating <- factor(
  df_final$Mating,
  levels = c("selfing",
             "mixed",
             "outcrossing")
)

# Reformating the imputed dataset
# The initial dataset uses the one-hot encoding
# For discrete traits with more than two states this is not appropriate for evolutionary model fitting
# In addition for "paired traits" due to one-hot encoding, only one need to be simulated
# As the re-coding depend on the traits, this is done manually
# Note also that for traits with more than two states they can be ordered and some transitions can be set to 0 (see below)

# Split the initial matrix into continuous and discrete trait matrices
species_names <- df_imputed$Species
continous_traits <- df_imputed[,c(2:8)]
discrete_traits <- df_imputed[,c(9:43)]

# List of states per trait
lapply(df_final,levels)

## Recoding onehot to multistate ----

# To transform onehot encoded variable into multistate variable
# state 1 = 0, state 2 = 2, polymorphic state = 1
onehot2multi <- function(w,x,y=0,z=0) {
  return(2*x - w*x +5*y -2*w*y - 3*x*y + 9*z - 3*w*z -4*x*z -6*y*z)
}

# Two-state traits
woodiness <- as.factor(onehot2multi(df_imputed$Woodiness_woody,df_imputed$Woodiness_herbaceous))
climbing <- as.factor(onehot2multi(df_imputed$Climbing_non.climbing,df_imputed$Climbing_climbing))
aquatic <- as.factor(onehot2multi(df_imputed$Aquatic_non.aquatic,df_imputed$Aquatic_non.aquatic))
sexualsystem <- as.factor(onehot2multi(df_imputed$SexualSystem_monomorphic,df_imputed$SexualSystem_dimorphic))
lifespan <- as.factor(onehot2multi(df_imputed$Lifespan_long,df_imputed$Lifespan_short))
dispersaldist <- as.factor(onehot2multi(df_imputed$DispersalDist_long,df_imputed$DispersalDist_short))
flowersex <- as.factor(onehot2multi(df_imputed$FlowerSex_bisexual,df_imputed$FlowerSex_unisexual))
showiness <- as.factor(onehot2multi(df_imputed$Showiness_bright,df_imputed$Showiness_dull))

# Three-state traits
pollination <- as.factor(onehot2multi(df_imputed$Pollination_autonomous,df_imputed$Pollination_abiotic,df_imputed$Pollination_biotic))
dispersalmode <- as.factor(onehot2multi(df_imputed$DispersalMode_autonomous,df_imputed$DispersalMode_abiotic,df_imputed$DispersalMode_biotic))
ovaryposition <- as.factor(onehot2multi(df_imputed$OvaryPosition_superior,df_imputed$OvaryPosition_intermediate,df_imputed$OvaryPosition_inferior))
flowersymmetry <- as.factor(onehot2multi(df_imputed$FlowerSymmetry_actinomorphic,df_imputed$FlowerSymmetry_zygomorphic,df_imputed$FlowerSymmetry_other))
floralreward <- as.factor(onehot2multi(df_imputed$FloralReward_none,df_imputed$FloralReward_pollen,df_imputed$FloralReward_nectar,df_imputed$FloralReward_other))

# Special case - polymorphic cases are considered as mixed
matingsystem <- as.factor(onehot2multi(df_imputed$Mating_selfing,df_imputed$Mating_mixed,df_imputed$Mating_outcrossing))
matingsystem <- as.factor(ifelse(matingsystem=="3" | matingsystem=="4",2,matingsystem))

discrete_traits_recoded <- data.frame(
  list(
    "woodiness"=woodiness,
    "climbing"=climbing,
    "aquatic"=aquatic,
    "dioecy"=sexualsystem,
    "flowersex"=flowersex,
    "lifespan"=lifespan,
    "showiness"=showiness,
    "dispersaldist"=dispersaldist,
    "dispersalmode"=dispersalmode,
    "matingsystem"=matingsystem,
    "pollination"=pollination,
    "ovaryposition"=ovaryposition,
    "flowersymetry"=flowersymmetry,
    "reward"=floralreward
  )
)

# Discrete traits as factor for the gower distance analysis
df_recoded <- cbind(
  continous_traits,
  discrete_traits_recoded
)

# save data set for use in simulations
saveRDS(df_recoded,"outputs/13_imputed_recoded.rds")

# Lists of  names
traits_list_final <- names(df_final)
traits_list_onehot <- names(df_onehot)
traits_list_recoded <- names(df_recoded)

## Computing distance matrices ----
dist_final <- daisy(df_final,"gower")
dist_onehot <- daisy(df_onehot,"gower")
dist_recoded <- daisy(df_recoded,"gower")

# Run pcoa and keep the 4 first axes
pcoa_final <- data.frame(pcoa(dist_final)$vector[,c(1:4)])
pcoa_onehot <- data.frame(pcoa(dist_onehot)$vector[,c(1:4)])
pcoa_recoded <- data.frame(pcoa(dist_recoded)$vector[,c(1:4)])

####
## Plot all traits on PCOA ----
####

# Function to generate plot for all trait using
plot_figure_pcoa <- function(df,tr,suffix,label_list) {
  names(df) <- c("dim1","dim2","dim3","dim4")
  label <- label_list[,tr]
  if(is.numeric(label)) {
    # color <- scale_colour_gradientn(colours=terrain.colors(10),name=NULL)
    color <- scale_colour_viridis_c()
  } else
    color <- scale_color_discrete()
  g12 <- ggplot(data=df,aes(x=dim1,y=dim2,col=label)) + geom_point() + color
  g13 <- ggplot(data=df,aes(x=dim1,y=dim3,col=label)) + geom_point() + color
  g14 <- ggplot(data=df,aes(x=dim1,y=dim4,col=label)) + geom_point() + color
  g23 <- ggplot(data=df,aes(x=dim2,y=dim3,col=label)) + geom_point() + color
  g24 <- ggplot(data=df,aes(x=dim2,y=dim4,col=label)) + geom_point() + color
  g34 <- ggplot(data=df,aes(x=dim3,y=dim4,col=label)) + geom_point() + color
  fig <- ggarrange(
    plotlist = list(g12,g13,g14,g23,g24,g34),
    ncol = 2,
    nrow = 3,
    common.legend = T
  )
  annotate_figure(fig,top = text_grob(tr, face = "bold", size = 14))
  ggsave(paste0("figures/pcoa/",tr,"_pcoa_",suffix,".pdf"),width = 7,height = 10)
}

for(tr in traits_list_final) {
  plot_figure_pcoa(df = pcoa_final,tr = tr,suffix = "final",label_list = df_final)
  plot_figure_pcoa(df = pcoa_onehot,tr = tr,suffix = "onehot",label_list = df_final)
  plot_figure_pcoa(df = pcoa_recoded,tr = tr,suffix = "recoded",label_list = df_final)
}

####
## Plot all traits on UMAP ----
####

# The umap method can be applied either on the orginal dataset or on the distance matrix
# We applied it on the distance matrix as it doesn't handle missing data
plot_figure_umap<- function(df,tr,suffix,label_list) {
  names(df) <- c("dim1","dim2")
  label <- label_list[,tr]
  if(is.numeric(label)) {
    #color <- scale_colour_gradientn(colours=terrain.colors(10),name=NULL)
    color <- scale_colour_viridis_c()
  } else
    color <- scale_color_discrete()
  ggplot(data=df,aes(x=dim1,y=dim2,col=label)) + geom_point() + color
  ggsave(paste0("figures/umap/",tr,"_umap_knn",knn,"_",suffix,".pdf"),width = 7,height = 7)
}

# The number of neighbors play a key role in the importance of the small scale versus the large scale
# A there are too many plots so only the final distance matrix is used
# But onehot and recoded matrices can be also run by uncommenting the following code

### Plot all traits on UMAP in folder ----
for(knn in c(10,25,50,100)) {
  custom_config <- umap.defaults
  custom_config$n_components <-  2# number of dimensions targeted
  custom_config$n_neighbors <- knn # number of dimensions targeted
  custom_config$input <- "dist" # The input matrix is a distance matrix
  umap_final <- umap(d = as.matrix(dist_final),config = custom_config)
  df_umap_final <- data.frame(umap_final$layout)
  rownames(df_umap_final) <- species_names
  #umap_onehot <- umap(d = as.matrix(dist_onehot),config = custom_config)
  #df_umap_onehot <- data.frame(umap_onehot$layout)
  #rownames(df_umap_onehot) <- species_names
  #umap_recoded <- umap(d = as.matrix(dist_recoded),config = custom_config)
  #df_umap_recoded <- data.frame(umap_recoded$layout)
  #rownames(df_umap_recoded) <- species_names
  for(tr in traits_list_final) {
    plot_figure_umap(df = df_umap_final,tr = tr,suffix = "umap_final",label_list = df_final)
    #plot_figure_umap(df = df_umap_onehot,tr = tr,suffix = "umap_onehot",label_list = df_final)
    #plot_figure_umap(df = df_umap_recoded,tr = tr,suffix = "umap_recoded",label_list = df_final)
  }
}

####
# Figure S9: UMAP knn10-100 ----
####

#UMAP config
custom_config <- umap.defaults
custom_config$n_components <-  2# number of dimensions targeted
custom_config$input <- "dist" # The input matrix is a distance matrix

#labels
FlowerSex <- df_final$FlowerSex
Woodiness <- df_final$Woodiness

#knn 100
custom_config$n_neighbors <- 100 # number of dimensions targeted
umap_final <- umap(d = as.matrix(dist_final),config = custom_config)
df_umap_final <- data.frame(umap_final$layout)
rownames(df_umap_final) <- species_names
colnames(df_umap_final) <- c("dim1","dim2")

p100<-ggplot(data=df_umap_final,aes(x=dim1,y=dim2,col=FlowerSex,shape=Woodiness)) + geom_point(size=3,alpha=0.6) + scale_color_discrete() +
  theme(legend.position = c(0.2, 0.25), legend.background = element_rect(fill=NA)) + ggtitle("knn = 100")
p100

#knn 50
custom_config$n_neighbors <- 50 # number of dimensions targeted
umap_final <- umap(d = as.matrix(dist_final),config = custom_config)
df_umap_final <- data.frame(umap_final$layout)
rownames(df_umap_final) <- species_names
colnames(df_umap_final) <- c("dim1","dim2")

p50<-ggplot(data=df_umap_final,aes(x=dim1,y=dim2,col=FlowerSex,shape=Woodiness)) + geom_point(size=3,alpha=0.6) + scale_color_discrete() + theme(legend.position = "none") + ggtitle("knn = 50")

#knn 25
custom_config$n_neighbors <- 25 # number of dimensions targeted
umap_final <- umap(d = as.matrix(dist_final),config = custom_config)
df_umap_final <- data.frame(umap_final$layout)
colnames(df_umap_final) <- c("dim1","dim2")

p25<-ggplot(data=df_umap_final,aes(x=dim1,y=dim2,col=FlowerSex,shape=Woodiness)) + geom_point(size=3,alpha=0.6) + scale_color_discrete() + theme(legend.position = "none") + ggtitle("knn = 25")

#knn 10
custom_config$n_neighbors <- 10 # number of dimensions targeted
umap_final <- umap(d = as.matrix(dist_final),config = custom_config)
df_umap_final <- data.frame(umap_final$layout)
colnames(df_umap_final) <- c("dim1","dim2")

p10<-ggplot(data=df_umap_final,aes(x=dim1,y=dim2,col=FlowerSex,shape=Woodiness)) + geom_point(size=3,alpha=0.6) + scale_color_discrete() + theme(legend.position = "none") + ggtitle("knn = 10")

# make combined plot
( p100 + p50 ) /
  ( p25 + p10 )

ggsave("figures/figure_S10_scatterplots_umap_knn10-100.png",width = 15, height=15)
