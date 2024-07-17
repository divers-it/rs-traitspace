# https://rpubs.com/pjmurphy/269609

#load data set
df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

# split data set into quant and qualitative columns
df_num <- df[sapply(df,is.numeric)]
df_cat <- df[sapply(df,is.factor)]

### hierarchical clustering for variables ----

# do clustering and plot tree
tree <- ClustOfVar::hclustvar(df_num, df_cat)
plot(tree)

# estimate number of clusters using stability of partitions
stab <- ClustOfVar::stability(tree,
                              B=50, # number of bootstraps
                              graph = FALSE
                              )

# plot stability results
plot(stab)
boxplot(stab$matCR[,1:7])

### K-means clustering for variables ----

# # make new df with unique colnames
# kdf_cat <- df_cat
# 
# #make sure that categories are not the same
# kdf_cat$col_Lifespan <- as.factor(paste("life_", kdf_cat$col_Lifespan,sep=""))
# kdf_cat$col_Pollination <- as.factor(paste("poll_", kdf_cat$col_Pollination,sep=""))
# kdf_cat$col_FlowerSymmetry <- as.factor(paste("symm_", kdf_cat$col_FlowerSymmetry,sep=""))
# str(kdf_cat)
# 
# 
# k.means <- ClustOfVar::kmeansvar(X.quanti = df_num, X.quali = kdf_cat, init = 4)
# summary(k.means)
# 
# # NOTE: There is an issue where NAs are treated as separate levels
# # so datasets with NAs will not run as package thinks these are shared 
# # states across traits: 
# vect.all.levels<-as.character(unlist(apply(kdf_cat,2,unique)))
# vect.all.levels.unique<-unique(vect.all.levels)
# test.name.categ<-length(vect.all.levels)==length(vect.all.levels.unique) 
# test.name.categ
