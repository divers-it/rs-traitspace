rm(list=ls())

#FROM:
#https://frbcesab.github.io/workshop-free/practice.html

#load packages
library(tidyr)
library(dplyr)
library(mFD)
library(ggplot2)
library(funrar)

#load formatted data
df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

#empty trait code vector
trait_code<-vector()

#read in columns and determine coding (NO ORDINAL YET)
for(i in 1:length(colnames(df))){
  if(is.numeric(df[,i])){
    trait_code[i]<-"Q"
  }
  
  if(is.factor(df[,i])){
    trait_code[i]<-"N"
  }
  
}

#make data frame with trait types
trait_code_df<-data.frame(colnames(df),trait_code)
colnames(trait_code_df)<-c("trait_name","trait_type")

#calculate functional distances using gower
sp_dist <- mFD::funct.dist(
  sp_tr         = df,
  tr_cat        = trait_code_df,
  metric        = "gower",
#  scale_euclid  = "scale_center", #already done
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = F)

#do PCoA and evaluate quality of space as axes are added
fspaces_quality <- mFD::quality.fspaces(
  sp_dist             = sp_dist,
  maxdim_pcoa         = 15,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "ward.D")

## Quality metrics of functional spaces ----

#The space with the best quality has the lowest quality metric.
round(fspaces_quality$"quality_fspaces", 3)
#write.csv(round(fspaces_quality$"quality_fspaces", 3),"outputs/mfd_qual.csv")

#With the mFD package, it is possible to illustrate the quality of PCoA-based multidimensional spaces according
#to deviation between trait-based distances and distances in the functional space

#This function generates a figure with three panels (in rows) for each selected functional space (in columns). 
#Each column represents a functional space, the value of the quality metric is written on the top of each column. 
#The x-axis of all panels represents trait-based distances. The y-axis is different for each row:

#on the first (top) row, the y-axis represents species functional distances in the multidimensional space. 
#Thus, the closer species are to the 1:1 line, the better distances in the functional space fit trait-based ones

#on the second row, the y-axis shows the raw deviation of species distances in the functional space compared to 
#trait-based distances. Thus, the raw deviation reflects the distance to the horizontal line.

#on the third row (bottom), the y-axis shows the absolute or squared deviation of the (“scaled”) distance in 
#the functional space. It is the deviation that is taken into account for computing the quality metric.
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

ggsave("figures/10_scatterplots_mfd_fspace_quality.png",width=30,height=15,units = "cm")

## Testing the correlation between functional axes and traits ----

#mFD allows to test for correlations between traits and functional axes and then illustrate possible correlations 
#continuous traits = linear model is computed and r2 and associated p-value are returned
#non-continuous traits = Kruskal-Wallis test is computed and eta2 statistic is returned

#a matrix of species coordinates taken from the output
sp_faxes_coord <- fspaces_quality$"details_fspaces"$"sp_pc_coord"

#correlation of continuous traits
df_faxes <- mFD::traits.faxes.cor(
  sp_tr          = df[,c(1:7)], 
  sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4")],
  stop_if_NA = F,
  plot = TRUE)

# Print traits with significant effect
df_faxes$"tr_faxes_stat"[which(df_faxes$"tr_faxes_stat"$"p.value" < 0.05), ]

#plot
df_faxes$"tr_faxes_plot"
ggsave("figures/10_scatterplots_mfd_quant_traits_vs_axes.png",width=30,height=20,units = "cm")

#correlation of discrete traits 8-14
df_faxes <- mFD::traits.faxes.cor(
  sp_tr          = df[,c(8:14)], 
  sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4")],
  stop_if_NA = F,
  plot = TRUE)

# Print traits with significant effect
df_faxes$"tr_faxes_stat"[which(df_faxes$"tr_faxes_stat"$"p.value" < 0.05), ]

#plot
df_faxes$"tr_faxes_plot"
ggsave("figures/10_scatterplots_mfd_qual_traits_vs_axes_1.png",width=30,height=20,units = "cm")

#correlation of qualitative traits 15-21
df_faxes <- mFD::traits.faxes.cor(
  sp_tr          = df[,c(15:21)], 
  sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4")],
  stop_if_NA = F,
  plot = TRUE)

# Print traits with significant effect
df_faxes$"tr_faxes_stat"[which(df_faxes$"tr_faxes_stat"$"p.value" < 0.05), ]

#plot
df_faxes$"tr_faxes_plot"
ggsave("figures/10_scatterplots_mfd_qual_traits_vs_axes_2.png",width=30,height=20,units = "cm")

## Plotting the selected functional space and position of species ----

#get coordinates
sp_faxes_coord <- fspaces_quality$"details_fspaces"$"sp_pc_coord"

#plot
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
ggsave("figures/10_scatterplots_mfd_functional_space.png",width=30,height=30,units = "cm")

## Computing and plotting alpha FD indices ----

#read in clustering
clust_ward<-readRDS("outputs/11_clust_num_k_2_7_ward.rds")

#make df with clustering info
clust_ward_df<-as.data.frame(cbind(rownames(clust_ward),clust_ward$`3clusters`))
colnames(clust_ward_df)<-c("species","cluster")

#recode df into one-hot with species as columns
clust_ward_recode <- clust_ward_df %>% mutate(value = 1)  %>% spread(species, value,  fill = 0 ) 

#remove cluster label column and add as name
cw<-clust_ward_recode[,c(2:length(colnames(clust_ward_recode)))]
rownames(cw)<-clust_ward_recode$cluster

#check names
table(colnames(cw)==rownames(sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4")]))

#compute 5 functional indices
alpha_fd_indices <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = as.matrix(sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4")]),
  asb_sp_w         = as.matrix(cw),
  #ind_vect         = c("fdis", "fric", "fdiv","fspe", "fide"), #output all
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

#output indices
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


write.csv(fd_ind_table,"outputs/10_mfd_indices_wardD2_k3.csv")

#FDis Functional Dispersion: the biomass weighted deviation of species traits values from the center of the functional space filled by the #assemblage i.e. the biomass-weighted mean distance to the biomass-weighted mean trait values of the assemblage.
#
#FRic Functional Richness: the proportion of functional space filled by species of the studied assemblage, i.e. the volume inside the convex#-hull shaping species. To compute FRic the number of species must be at least higher than the number of functional axis + 1.
#
#FDiv Functional Divergence: the proportion of the biomass supported by the species with the most extreme functional traits i.e. the ones #located close to the edge of the convex-hull filled by the assemblage.
#
#FEve Functional Evenness: the regularity of biomass distribution in the functional space using the Minimum Spanning Tree linking all species #present in the assemblage.
#
#FSpe Functional Specialization: the biomass weighted mean distance to the mean position of species from the global pool (present in all #assemblages).
#
#FMPD Functional Mean Pairwise Distance: the mean weighted distance between all species pairs.
#
#FNND Functional Mean Nearest Neighbour Distance: the weighted distance to the nearest neighbor within the assemblage.
#
#FIde Functional Identity: the mean traits values for the assemblage. FIde is always computed when FDis is computed.
#
#FOri Functional Originality: the weighted mean distance to the nearest species from the global species pool.

#information such as coordinates of centroids, distances and identity of the nearest neighbour, 
#distances to the centroid, etc. The user does not have to directly use it but it will be useful 
#if FD indices are then plotted. It can be retrieved through:
details_list <- alpha_fd_indices$"details"
details_list

#plot functional indices of clusters 1 and 2
plots_alpha <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices,
  plot_asb_nm              = c("k_3_cluster_1", "k_3_cluster_2"),
  ind_nm                   = c("fdis", "fric", "fdiv", 
                               "fspe", "fide"),
  faxes                    = NULL,
  faxes_nm                 = NULL,
  range_faxes              = c(NA, NA),
  plot_sp_nm               = NULL,
  save_file                = FALSE,
  check_input              = TRUE) 


#FRic representation: the colored shapes reflect the convex-hull of the studied assemblages
#and the white shape reflects the convex-hull of the global pool of species:

plots_alpha$"fric"$"patchwork"
ggsave("figures/10_scatterplots_mfd_funct_rich_wardD2_k3_clust1_clust2.png",width=30,height=30,units = "cm")

#FDiv representation: the gravity centers of vertices (i.e. species with the most extreme functional traits) of each 
#assemblages are plotted as a square and a triangle. The two colored circles represent the mean
#distance of species to the gravity center for each assemblage. Species of each assemblage 
#have different size given their relative weight into the assemblage.
plots_alpha$"fdiv"$"patchwork"
ggsave("figures/10_scatterplots_mfd_funct_div_wardD2_k3_clust1_clust2.png",width=30,height=30,units = "cm")

#FSpe representation: colored traits represent distances of each species from a given assemblage 
#to the center of gravity of the global pool (i.e center of the functional space). the center of
#gravity is plotted with a purple diamond. Species of each assemblage have different size given
#their relative weight into the assemblage.
plots_alpha$"fdis"$"patchwork"
ggsave("figures/10_scatterplots_mfd_funct_disp_wardD2_k3_clust1_clust2.png",width=30,height=30,units = "cm")

#FIde representation:colored lines refer to the weighted average position of species of each assemblage
#along each axis. Species of each assemblage have different size given their relative weight
#into the assemblage.
plots_alpha$"fide"$"patchwork"
ggsave("figures/10_scatterplots_mfd_funct_ident_wardD2_k3_clust1_clust2.png",width=30,height=30,units = "cm")

## Functional originality at regional scale ----
sp_di <- distinctiveness_global(sp_dist, di_name = "distinctiveness")

#We get one value of distinctiveness per species. 
#It only considers the functional dissimilarity of all species in the dissimilarity matrix without
#considering their spatial distributions. We get one value of distinctiveness per species.
summary(sp_di)

quantile(sp_di$distinctiveness, probs = seq(0, 1, by = 0.1))
subset(sp_di, distinctiveness >= 0.439)

#For the choice or dissimilarity matrix we can use the raw dissimilarity matrix computed directly on raw 
#traits values among species, as we did here. Another option would be to compute a new functional 
#dissimilarity matrix based on the selected functional axes. One advantage of the latter is that it 
#already takes into account the correlation between traits.

#Let’s recompute regional functional distinctiveness based on the four selected functional axes. 
#Because the space comes from a PCA, we can directly use euclidean distance.

new_dissim <- dist(sp_faxes_coord[, c("PC1", "PC2", "PC3", "PC4")])

sp_di_alt <- distinctiveness_global(new_dissim, di_name = "alt_di")

#We can now compare both distinctiveness values.
sp_all_di <- merge(sp_di, sp_di_alt, by = "species")

plot(sp_all_di$distinctiveness, sp_all_di$alt_di)

cor.test(sp_all_di$distinctiveness, sp_all_di$alt_di)

#Both seems very correlated, so in our case using either one should be fine. However, it can be better 
#to use dissimilarity based on a reduced number of well-defined axes because: (1) there are more 
#interpretable thanks to the multivariate analysis, (2) the first ones contain the most information, 
#(3) they explicitly take into account potentially strong correlations between provided traits. 
#We’ll stick here with raw dissimilarity for the sake of simplicity.

##make empty data frame
df_w<-data.frame(matrix(nrow=3, ncol=length(df$Maximumverticalheight)))
colnames(df_w)<-rownames(df)
#
##put 0/1 rows denoting group membership
df_w[1,]<-as.numeric(colnames(df_w)%in%rownames(clust_ward[clust_ward$`3clusters`=="k_3_cluster_1",]))
df_w[2,]<-as.numeric(colnames(df_w)%in%rownames(clust_ward[clust_ward$`3clusters`=="k_3_cluster_2",]))
df_w[3,]<-as.numeric(colnames(df_w)%in%rownames(clust_ward[clust_ward$`3clusters`=="k_3_cluster_3",]))

#To compute uniqueness at regional scale we also need the regional level functional dissimilarity matrix 
#with the uniqueness() function, and the site-species matrix:
sp_ui <- uniqueness(
  pres_matrix = as.matrix(df_w),
  as.matrix(sp_dist)
)

head(sp_ui)
quantile(sp_ui$Ui, probs = seq(0, 1, by = 0.1))

#the most isolated species in the functional space. Meaning that they have the most distant nearest neighbors.
subset(sp_ui, Ui >= 0.21)

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
  #ggrepel::geom_text_repel(aes(label = species)) +
  scale_color_viridis_c("Functional\nDistinctiveness") +
  theme_bw()

plot_reg_uniqueness <- ggplot(sp_coord_di_ui, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_point(aes(color = Ui)) +
  #ggrepel::geom_text_repel(aes(label = species)) +
  scale_color_viridis_c("Functional\nUniqueness") +
  theme_bw()

plot_reg_distinctiveness
ggsave("figures/10_scatterplot_mfd_distinctiveness.png", width= 8, height = 8)

plot_reg_uniqueness
ggsave("figures/10_scatterplot_mfd_uniqueness.png", width= 8, height = 8)

####
## NOT RUN: phylo stuff
####
#
#library(ape)
#phy<-read.tree("outputs/pruned_tree.tre")
#
#pdf("figures/10_phy.pdf",height = 20, width = 20)
#plot(phy,cex=0.5)
#nodelabels(cex=0.5)
#dev.off()
#
##split species into monocots, magnoliids and dicots
##NODES NEED TO BE CHANGED FOR NEW TREE
#
#monocots_phy<-extract.clade(phy,node=609)
#monocots_phy$tip.label
#mag_phy<-extract.clade(phy,node=665)
#mag_phy$tip.label
#dicot_phy<-extract.clade(phy,node=347)
#dicot_phy$tip.label
#
##species not in previous three categories
#other<-setdiff(sort(phy$tip.label),sort(c(monocots_phy$tip.label,
#                                    mag_phy$tip.label,
#                                    dicot_phy$tip.label)))
#
##make empty data frame
#df_w<-data.frame(matrix(, nrow=4, ncol=length(phy$tip.label)))
#colnames(df_w)<-sort(phy$tip.label)
#
##put 0/1 rows denoting group membership
#df_w[1,]<-as.numeric(colnames(df_w)%in%monocots_phy$tip.label)
#df_w[2,]<-as.numeric(colnames(df_w)%in%mag_phy$tip.label)
#df_w[3,]<-as.numeric(colnames(df_w)%in%dicot_phy$tip.label)
#df_w[4,]<-as.numeric(colnames(df_w)%in%other)
#
##rename rows
#rownames(df_w)<-c("monocots","mag","dicots","other")
#
##remove underscore to match df names
#colnames(df_w)<-gsub("_"," ",sort(phy$tip.label))
#
##remove species not present in df
#df_w<-df_w[, colnames(df_w) %in% rownames(df)]
#
##compute 5 functional indices
#alpha_fd_indices <- mFD::alpha.fd.multidim(
#  sp_faxes_coord   = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4")],
#  asb_sp_w         = as.matrix(df_w),
#  ind_vect         = c("fdis", "fric", "fdiv", 
#                       "fspe", "fide"),
#  scaling          = TRUE,
#  check_input      = TRUE,
#  details_returned = TRUE)
