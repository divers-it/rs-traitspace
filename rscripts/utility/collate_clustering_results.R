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
colnames(rob_df)[1]<-"ward"

### k prototype

kpro <- readRDS(file = here::here("outputs/11.1_robust_vect_kpro_full.rds"))

#check names
rownames(rob_df)==names(kpro)

#merge
rob_df<-cbind(rob_df,kpro)

### PAM

pam <- readRDS(file = here::here("outputs/11.2_robust_vect_pam_full.rds"))

#check names
rownames(rob_df)==names(pam)

#merge
rob_df<-cbind(rob_df,pam)

### LCM

lcm <- readRDS(file = here::here("outputs/11.3_clust_vect_lcm.rds"))

#check names
gsub(" ","_",rownames(rob_df))==names(lcm)

#merge
rob_df<-cbind(rob_df,lcm)

### Density clustering

dens <- readRDS(file = here::here("outputs/11.4_robust_vect_dens_full.rds"))

#check names
rownames(rob_df)==names(dens)

#merge
rob_df<-cbind(rob_df,dens)
rob_df

###
# One-hot
###

### Ward

ward <- readRDS(file = here::here("outputs/one_hot_11_clust_num_k_2_7_ward.rds"))

rownames(rob_df)==rownames(ward)

# no robust groups as hierarchical
# k = 4 seemed the best 
rob_df<-cbind(rob_df,ward$`4clusters`)

#get rid of prefix
rob_df$`ward$\`4clusters\`` <- as.numeric(gsub("k_4_cluster_","",rob_df$`ward$\`4clusters\``))

#fix names
colnames(rob_df)[length(colnames(rob_df))]<-"one_hot_ward"

rob_df

### k prototype

one_hot_kpro <- readRDS(file = here::here("outputs/one_hot_11.1_robust_vect_kpro_full.rds"))

#check names
rownames(rob_df)==names(one_hot_kpro)

#merge
rob_df<-cbind(rob_df,one_hot_kpro)

### PAM

one_hot_pam <- readRDS(file = here::here("outputs/one_hot_11.2_robust_vect_pam_full.rds"))

#check names
rownames(rob_df)==names(one_hot_pam)

#merge
rob_df<-cbind(rob_df,one_hot_pam)

### LCM

one_hot_lcm <- readRDS(file = here::here("outputs/one_hot_11.3_clust_vect_lcm.rds"))

#check names
gsub(" ","_",rownames(rob_df))==names(one_hot_lcm)

#merge
rob_df<-cbind(rob_df,one_hot_lcm )

### Density clustering

one_hot_dens <- readRDS(file = here::here("outputs/one_hot_11.4_robust_vect_dens_full.rds"))

#check names
rownames(rob_df)==names(one_hot_dens)

#merge
rob_df<-cbind(rob_df,one_hot_dens)

#compare normal to one hot
table(rob_df$ward,rob_df$one_hot_ward)
table(rob_df$pam,rob_df$one_hot_pam)
table(rob_df$kpro,rob_df$one_hot_kpro)
table(rob_df$lcm,rob_df$one_hot_lcm)

#output csv
rownames(rob_df)==rownames(df)
write.csv(cbind(rob_df,df),"outputs/collate_clustering_results.csv")
