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

###
# One-hot
###

### Ward

ward <- readRDS(file = here::here("outputs/clust_num_k_2_7_ward_one_hot.rds"))

rownames(rob_df)==rownames(ward)

# no robust groups as hierarchical
# k = 4 seemed the best 
rob_df<-cbind(rob_df,ward$`4clusters`)

#get rid of prefix
rob_df$`ward$\`4clusters\`` <- as.numeric(gsub("k_4_cluster_","",rob_df$`ward$\`4clusters\``))

#fix names
colnames(rob_df)[length(colnames(rob_df))]<-"ward_one_hot"

rob_df

### k prototype

kpro_one_hot <- readRDS(file = here::here("outputs/robust_vect_kpro_full_one_hot.rds"))

#check names
rownames(rob_df)==names(kpro_one_hot)

#merge
rob_df<-cbind(rob_df,kpro_one_hot)

### PAM

pam_one_hot <- readRDS(file = here::here("outputs/robust_vect_pam_full_one_hot.rds"))

#check names
rownames(rob_df)==names(pam_one_hot)

#merge
rob_df<-cbind(rob_df,pam_one_hot)

### Density clustering

dens_one_hot <- readRDS(file = here::here("outputs/robust_vect_dens_full_one_hot.rds"))

#check names
rownames(rob_df)==names(dens_one_hot)

#merge
rob_df<-cbind(rob_df,dens_one_hot)

#compare normal to one hot
table(rob_df$ward,rob_df$ward_one_hot)
table(rob_df$pam,rob_df$pam_one_hot)
table(rob_df$kpro,rob_df$kpro_one_hot)

#output csv
rownames(rob_df)==rownames(df)
write.csv(cbind(rob_df,df),"outputs/collate_clustering_results.csv")
