#FROM:
#https://frbcesab.github.io/workshop-free/practice.html

#load formatted data
df<-readRDS(file = here::here("outputs/df_filt_trans.rds"))

#load Diaz data
diaz<-read.csv("data/diaz_species_mean_traits.csv")

#load substitutes
subs<-readRDS(file = here::here("outputs/sub_genera.rds"))

rownames(diaz)<-diaz$Species_name_standardized_against_TPL

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

#add congenerics (this is done after filtering to only those with data for >x traits)
diaz_pcf_ds<-diaz_pcf[rownames(diaz_pcf)%in%c(subs$Species_name_standardized_against_TPL,rownames(df)),]

#empty trait code vector
trait_code<-vector()

#read in columns and determine coding (NO ORDINAL YET)
for(i in 1:length(colnames(diaz_pcf_ds))){
  if(is.numeric(diaz_pcf_ds[,i])){
    trait_code[i]<-"Q"
  }
  
  if(is.factor(diaz_pcf_ds[,i])){
    trait_code[i]<-"N"
  }
  
}

#make data frame
trait_code_df<-data.frame(colnames(diaz_pcf_ds),trait_code)
colnames(trait_code_df)<-c("trait_name","trait_type")

library(mFD)

sp_dist <- mFD::funct.dist(
  sp_tr         = data.frame(diaz_pcf_ds),
  tr_cat        = trait_code_df,
  metric        = "gower",
  scale_euclid  = "scale_center",
  ordinal_var   = "classic",
  weight_type   = "equal",
  stop_if_NA    = F)


fspaces_quality <- mFD::quality.fspaces(
  sp_dist             = sp_dist,
  maxdim_pcoa         = 15,
  deviation_weighting = "absolute",
  fdist_scaling       = FALSE,
  fdendro             = "average")

## Quality metrics of functional spaces ----

#The space with the best quality has the lowest quality metric.
round(fspaces_quality$"quality_fspaces", 3)
write.csv(round(fspaces_quality$"quality_fspaces", 3),"outputs/mfd_qual_diaz.csv")

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

png("figures/mfd_quality_diaz.png",height=1000,width=3000,res=200)
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
dev.off()

## Testing the correlation between functional axes and traits ----

#mFD allows to test for correlations between traits and functional axes and then illustrate possible correlations 
#continuous traits = linear model is computed and r2 and associated p-value are returned
#non-continuous traits = Kruskal-Wallis test is computed and eta2 statistic is returned

#a matrix of species coordinates taken from the output
sp_faxes_coord <- fspaces_quality$"details_fspaces"$"sp_pc_coord"

#correlation of continuous traits
df_faxes <- mFD::traits.faxes.cor(
  sp_tr          = data.frame(diaz_pcf_ds), 
  sp_faxes_coord = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4")],
  stop_if_NA = F,
  plot = TRUE)

# Print traits with significant effect
df_faxes$"tr_faxes_stat"[which(df_faxes$"tr_faxes_stat"$"p.value" < 0.05), ]

#plot
png("figures/mfd_traits_vs_axes_quant_diaz.png",width=2500,height = 2000,res=200)
df_faxes$"tr_faxes_plot"
dev.off()

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

png("figures/mfd_functional_space_diaz.png",width=2000,height = 2000, res=200)
big_plot$"patchwork"
dev.off()

## Computing and plotting alpha FD indices ----

library(ape)
phy<-read.tree("outputs/pruned_tree.tre")

#pdf("figures/phy.pdf",height = 20, width = 20)
#plot(phy,cex=0.5)
#nodelabels(cex=0.5)
#dev.off()

#split species into monocots, magnoliids and dicots
monocots_phy<-extract.clade(phy,node=609)
monocots_phy$tip.label
mag_phy<-extract.clade(phy,node=665)
mag_phy$tip.label
dicot_phy<-extract.clade(phy,node=347)
dicot_phy$tip.label

#species not in previous three categories
other<-setdiff(sort(phy$tip.label),sort(c(monocots_phy$tip.label,
                                    mag_phy$tip.label,
                                    dicot_phy$tip.label)))

#make empty data frame
df_w<-data.frame(matrix(, nrow=4, ncol=length(phy$tip.label)))
colnames(df_w)<-sort(phy$tip.label)

#put 0/1 rows denoting group membership
df_w[1,]<-as.numeric(colnames(df_w)%in%monocots_phy$tip.label)
df_w[2,]<-as.numeric(colnames(df_w)%in%mag_phy$tip.label)
df_w[3,]<-as.numeric(colnames(df_w)%in%dicot_phy$tip.label)
df_w[4,]<-as.numeric(colnames(df_w)%in%other)

#rename rows
rownames(df_w)<-c("monocots","mag","dicots","other")

#remove underscore to match df names
colnames(df_w)<-gsub("_"," ",sort(phy$tip.label))

#remove species not present in df
df_w<-df_w[, colnames(df_w) %in% rownames(df)]

#compute 5 functional indices
alpha_fd_indices <- mFD::alpha.fd.multidim(
  sp_faxes_coord   = sp_faxes_coord[ , c("PC1", "PC2", "PC3", "PC4")],
  asb_sp_w         = as.matrix(df_w),
  ind_vect         = c("fdis", "fric", "fdiv", 
                       "fspe", "fide"),
  scaling          = TRUE,
  check_input      = TRUE,
  details_returned = TRUE)

#output indices
fd_ind_values <- alpha_fd_indices$"functional_diversity_indices"
fd_ind_values
write.csv(fd_ind_values,"outputs/mfd_ind_values.csv")

#information such as coordinates of centroids, distances and identity of the nearest neighbour, 
#distances to the centroid, etc. The user does not have to directly use it but it will be useful 
#if FD indices are then plotted. It can be retrieved through:
details_list <- alpha_fd_indices$"details"
details_list

#plot functional indices
plots_alpha <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices,
  plot_asb_nm              = c("monocots", "dicots"),
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

png("figures/mfd_funct_rich_monocot_dicot.png",height = 1500, width = 1500,res=150)
plots_alpha$"fric"$"patchwork"
dev.off()

#FDiv representation: the gravity centers of vertices (i.e. species with the most extreme functional traits) of each 
#assemblages are plotted as a square and a triangle. The two colored circles represent the mean
#distance of species to the gravity center for each assemblage. Species of each assemblage 
#have different size given their relative weight into the assemblage.
png("figures/mfd_funct_div_monocot_dicot.png",height = 1500, width = 1500,res=150)
plots_alpha$"fdiv"$"patchwork"
dev.off()

#FSpe representation: colored traits represent distances of each species from a given assemblage 
#to the center of gravity of the global pool (i.e center of the functional space). the center of
#gravity is plotted with a purple diamond. Species of each assemblage have different size given
#their relative weight into the assemblage.

png("figures/mfd_funct_disp_monocot_dicot.png",height = 1500, width = 1500, res=150)
plots_alpha$"fdis"$"patchwork"
dev.off()

#FIde representation:colored lines refer to the weighted average position of species of each assemblage
#along each axis. Species of each assemblage have different size given their relative weight
#into the assemblage.

pngg("figures/mfd_funct_ident_monocot_dicot.png",,height = 1500, width = 1500, res=150)
plots_alpha$"fide"$"patchwork"
dev.off()


## Functional originality at regional scale ----
library(funrar)

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
#interpretable thanks to the multivariate analysis, (2) the first one contain de most information, 
#(3) they explicitly take into account potentially strong correlations between provided traits. 
#We’ll stick here with raw dissimilarity for the sake of simplicity.

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

#We already saw in part 1 how to visualize the functional space using the pre-made functions of mFD. Here we will use our own functions to be able to color species in function of their functional originality.

# Make a summary data.frame
sp_coord_di_ui <- as.data.frame(sp_faxes_coord[, 1:2])
sp_coord_di_ui$species <- rownames(sp_coord_di_ui)
rownames(sp_coord_di_ui) <- NULL
sp_coord_di_ui <- sp_coord_di_ui[, c(3, 1, 2)]
sp_coord_di_ui <- merge(sp_coord_di_ui, sp_di, by = "species")
sp_coord_di_ui <- merge(sp_coord_di_ui, sp_ui, by = "species")
library("ggplot2")

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
ggsave("figures/mfd_distinctiveness.png", width= 8, height = 8)

plot_reg_uniqueness
ggsave("figures/mfd_uniqueness.png", width= 8, height = 8)
