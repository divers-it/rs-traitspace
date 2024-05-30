###
# ---- Make table of robust groups / clustering from different approaches ----
###

#load formatted data
df <- readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

### Ward

ward <- readRDS(file = here::here("outputs/11_clust_num_k_2_7_ward.rds"))

# no robust groups as hierarchical
# k = 3 seemed the best 
rob_df<-data.frame(ward$`3clusters`)

#get rid of prefix
rob_df$ward..3clusters. <- as.numeric(gsub("k_3_cluster_","",rob_df$ward..3clusters.))

#fix names
rownames(rob_df)<-rownames(ward)
colnames(rob_df)[1]<-"ward_cluster"

### k prototype

kpro_robust <- readRDS(file = here::here("outputs/11.1_robust_vect_kpro_full.rds"))

#check names
rownames(rob_df)==names(kpro_robust)

#merge
rob_df<-cbind(rob_df,kpro_robust)

### PAM

pam_robust <- readRDS(file = here::here("outputs/11.2_robust_vect_pam_full.rds"))

#check names
rownames(rob_df)==names(pam_robust)

#merge
rob_df<-cbind(rob_df,pam_robust)

# ### LCM
# 
# lcm <- readRDS(file = here::here("outputs/11.3_clust_vect_lcm.rds"))
# 
# #check names
# gsub(" ","_",rownames(rob_df))==names(lcm)
# 
# #merge
# rob_df<-cbind(rob_df,lcm)

### Density clustering
# 
# dens <- readRDS(file = here::here("outputs/11.4_robust_vect_dens_full.rds"))
# 
# #check names
# rownames(rob_df)==names(dens)
# 
# #merge
# rob_df<-cbind(rob_df,dens)
# rob_df

###
# One-hot
###

### Ward

ward <- readRDS(file = here::here("outputs/one_hot_11_clust_num_k_2_7_ward.rds"))

rownames(rob_df)==rownames(ward)

# no robust groups as hierarchical
# k = 5 seemed the best 
rob_df<-cbind(rob_df,ward$`5clusters`)

#get rid of prefix
rob_df$`ward$\`5clusters\`` <- as.numeric(gsub("k_5_cluster_","",rob_df$`ward$\`5clusters\``))

#fix names
colnames(rob_df)[length(colnames(rob_df))]<-"one_hot_ward_cluster"

rob_df

### k prototype

one_hot_kpro_robust <- readRDS(file = here::here("outputs/one_hot_11.1_robust_vect_kpro_full.rds"))

#check names
sort(rownames(rob_df))==sort(names(one_hot_kpro_robust))

#merge
rob_df<-cbind(rob_df,one_hot_kpro_robust)

### PAM

one_hot_pam_robust <- readRDS(file = here::here("outputs/one_hot_11.2_robust_vect_pam_full.rds"))

#check names
rownames(rob_df)==names(one_hot_pam_robust)

#merge
rob_df<-cbind(rob_df,one_hot_pam_robust)

# Clustering methods not used in manuscript

# ### LCM
# 
# one_hot_lcm <- readRDS(file = here::here("outputs/one_hot_11.3_clust_vect_lcm.rds"))
# 
# #check names
# gsub(" ","_",rownames(rob_df))==names(one_hot_lcm)
# 
# #merge
# rob_df<-cbind(rob_df,one_hot_lcm )

# ### Density clustering
# 
# one_hot_dens <- readRDS(file = here::here("outputs/one_hot_11.4_robust_vect_dens_full.rds"))
# 
# #check names
# rownames(rob_df)==names(one_hot_dens)
# 
# #merge
# rob_df<-cbind(rob_df,one_hot_dens)

#compare normal to one hot
table(rob_df$ward,rob_df$one_hot_ward_cluster)
table(rob_df$pam,rob_df$one_hot_pam_robust)
table(rob_df$kpro,rob_df$one_hot_kpro_robust)

#output csv
rownames(rob_df)==rownames(df)
write.csv(cbind(rob_df),"outputs/collate_clustering_results.csv")
