# Make table of robust groups / clustering from different approaches

####
## Original encoding ----
####

# load formatted data
df <- readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

### Hierarchical Ward.D2 ----
ward <- readRDS(file = here::here("outputs/10.2_clust_num_k_2_7_ward.rds"))

# no robust groups as hierarchical
# k = 3 seemed the best 
rob_df<-data.frame(ward$`3clusters`)

# get rid of prefix
rob_df$ward..3clusters. <- as.numeric(gsub("k_3_cluster_","",rob_df$ward..3clusters.))

# fix names
rownames(rob_df)<-rownames(ward)
colnames(rob_df)[1]<-"ward_cluster"

### k-prototypes ----
kpro_robust <- readRDS(file = here::here("outputs/10.1_robust_vect_kpro_full.rds"))

# check names
table(rownames(rob_df)==names(kpro_robust))

# merge
rob_df<-cbind(rob_df,kpro_robust)

### PAM ----
pam_robust <- readRDS(file = here::here("outputs/10_robust_vect_pam_full.rds"))

# check names
table(rownames(rob_df)==names(pam_robust))

# merge
rob_df<-cbind(rob_df,pam_robust)

#### 
## One-hot encoding ----
#### 

### Hierarchical ward.D2 ----
ward <- readRDS(file = here::here("outputs/one_hot_10.2_clust_num_k_2_7_ward.rds"))

# check names
table(rownames(rob_df)==rownames(ward))

# no robust groups as hierarchical
# k = 5 seemed the best 
rob_df<-cbind(rob_df,ward$`5clusters`)

# get rid of prefix
rob_df$`ward$\`5clusters\`` <- as.numeric(gsub("k_5_cluster_","",rob_df$`ward$\`5clusters\``))

# fix names
colnames(rob_df)[length(colnames(rob_df))]<-"one_hot_ward_cluster"

### k-prototypes ----
one_hot_kpro_robust <- readRDS(file = here::here("outputs/one_hot_10.1_robust_vect_kpro_full.rds"))

# check names
sort(rownames(rob_df))==sort(names(one_hot_kpro_robust))

# merge
rob_df<-cbind(rob_df,one_hot_kpro_robust)

### PAM ----
one_hot_pam_robust <- readRDS(file = here::here("outputs/one_hot_10_robust_vect_pam_full.rds"))

# check names
table(rownames(rob_df)==names(one_hot_pam_robust))

## merge results ----
rob_df<-cbind(rob_df,one_hot_pam_robust)

# compare normal to one hot
table(rob_df$ward_cluster,rob_df$one_hot_ward_cluster)
table(rob_df$pam_robust,rob_df$one_hot_pam_robust)
table(rob_df$kpro_robust,rob_df$one_hot_kpro_robust)

# compare PAM to others
table(rob_df$pam_robust,rob_df$ward_cluster)
table(rob_df$pam_robust,rob_df$kpro_robust)

# output csv
table(rownames(rob_df)==rownames(df))
write.csv(cbind(rob_df),"outputs/collate_clustering_results.csv")
