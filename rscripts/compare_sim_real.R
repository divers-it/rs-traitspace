
library(GGally)

sim_list <- readRDS("outputs/phylo_simulated_datasets.rds")

#reorder
for(i in 1:length(sim_list)){
  sim_list[[i]]<-sim_list[[i]][order(rownames(sim_list[[i]])),]
}

df_sim <- as.data.frame(unclass(sim_list[[1]]), stringsAsFactors = TRUE)
rownames(df_sim)<-rownames(sim_list[[1]])

# Loading and formating datasets
df_final <- readRDS("outputs/6_df_filt_trans.rds")
df_onehot <- readRDS("outputs/one_hot_6_df_filt_trans.rds")
df_imputed <- read.csv("outputs/imputed_with_phylo.csv",row.names = 1,stringsAsFactors = TRUE)

#dissimilarity matrix calculation
gower_final <- daisy(df_final,
                          metric = "gower" )

gower_onehot <- daisy(df_onehot,
                     metric = "gower" )

gower_imputed <- daisy(df_imputed,
                     metric = "gower" )

gower_sim <- daisy(df_sim,
                       metric = "gower" )

gsub("_"," ",rownames(df_sim))==rownames(df_final)

plot(gower_sim~gower_onehot)

pairs(c(gower_final,gower_onehot,gower_imputed))


pairs(data.frame(list("original"=as.vector(gower_final),
    "onehot"=as.vector(gower_onehot),
    "imputed"=as.vector(gower_imputed),
    "simulated"=as.vector(gower_sim))),pch=16,cex=0.5
)


##plot pairs of variables (removing some with high cat count)
png("figures/ggpairs_sim.png",width=4000,height=4000,res=100)
ggpairs(df_sim,cardinality_threshold=16) 
dev.off()
