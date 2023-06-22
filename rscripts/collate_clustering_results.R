###
# Make table of robust groups / clustering from different approaches
###

#load formatted data
df <- readRDS(file = here::here("outputs/df_filt_trans.rds"))

### Ward

ward <- readRDS(file = here::here("outputs/clust_num_k_2_7_ward.rds"))

# no robust groups as hierarchical
# k = 3 seemed the best 
rob_df<-data.frame(ward$`3clusters`)

#get rid of prefix
rob_df$ward..3clusters. <- as.numeric(gsub("k_3_cluster_","",rob_df$ward..3clusters.))

#fix names
rownames(rob_df)<-rownames(ward)
colnames(rob_df)[1]<-"ward"

### k prototype

kpro <- readRDS(file = here::here("outputs/robust_vect_kpro_full.rds"))

#check names
rownames(rob_df)==names(kpro)

#merge
rob_df<-cbind(rob_df,kpro)

### PAM

pam <- readRDS(file = here::here("outputs/robust_vect_pam_full.rds"))

#check names
rownames(rob_df)==names(pam)

#merge
rob_df<-cbind(rob_df,pam)

### Density clustering

dens <- readRDS(file = here::here("outputs/robust_vect_dens_full.rds"))

#check names
rownames(rob_df)==names(dens)

#merge
rob_df<-cbind(rob_df,dens)

rob_df
