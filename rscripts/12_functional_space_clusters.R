#FROM:
#https://frbcesab.github.io/workshop-free/practice.html

#load formatted data
df<-readRDS(file = here::here("outputs/df_filt_trans.rds"))

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

#make data frame
trait_code_df<-data.frame(colnames(df),trait_code)
colnames(trait_code_df)<-c("trait_name","trait_type")

library(mFD)

sp_dist <- mFD::funct.dist(
  sp_tr         = df,
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

#a matrix of species coordinates taken from the output
sp_faxes_coord <- fspaces_quality$"details_fspaces"$"sp_pc_coord"

#load formatted data
df_w<-readRDS(file = here::here("outputs/clust.num.stack.Rds"))

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
write.csv(fd_ind_values,"outputs/mfd_ind_values_clusters_w5.csv")

#information such as coordinates of centroids, distances and identity of the nearest neighbour, 
#distances to the centroid, etc. The user does not have to directly use it but it will be useful 
#if FD indices are then plotted. It can be retrieved through:
details_list <- alpha_fd_indices$"details"
details_list

#plot functional indices
plots_alpha <- mFD::alpha.multidim.plot(
  output_alpha_fd_multidim = alpha_fd_indices,
  plot_asb_nm              = c("cluster1", "cluster2"),
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

png("figures/mfd_funct_ident_monocot_dicot.png",,height = 1500, width = 1500, res=150)
plots_alpha$"fide"$"patchwork"
dev.off()
