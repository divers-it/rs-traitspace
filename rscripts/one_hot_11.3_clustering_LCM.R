rm(list = ls())

#load packages
library(VarSelLCM)
library(ggplot2)
library(cluster)

#load data with NAs imputed
#NOTE: Run "impute_missing_data.R" to update missing data
df2<-read.csv("outputs/one_hot_imputed_with_phylo.csv",row.names = 1,stringsAsFactors = TRUE)

#NOTE: This could be modified but some columns that are quantitative and difficult
#to interpret / neutral in their expected effect on strategies appear to bias this
#method much more than the others
#remove biasing quantitative columns
df2<-subset(df2, select = -c(Number.of.fertile.stamens,Fusion.of.ovaries,Number.of.ovules.per.carpel))

# Cluster analysis without variable selection
res_without <- VarSelCluster(df2, gvals = 1:5, vbleSelec = FALSE, crit.varsel = "BIC")

# Cluster analysis with variable selection (with parallelisation)
res_with <- VarSelCluster(df2, gvals = 1:5, crit.varsel = "BIC")

#Comparison of the BIC for both models: variable selection can improve the BIC
#lower is better
BIC(res_without)
BIC(res_with)

# Estimated partition
fitted(res_with)

# Estimated probabilities of classification
head(fitted(res_with, type="probability"))

# Summary of the best model
summary(res_with)

# Discriminative power of the variables. The greater this index, the more the variable distinguishes the clusters.
# Too many for colours but this is coded into plot method
plot(res_with)

# Summary of qualitative variable
plot(res_with, y="Woodiness_herbaceous")

# Summary of quantitative variable
plot(res_with, y="Maximum.vertical.height")

# More detailed output
print(res_with)

# Print model parameter
coef(res_with)

#make clustering vector
lcm_res<-as.factor(fitted(res_with))
names(lcm_res)<-rownames(res_with@data@dataContinuous@notNA)
saveRDS(lcm_res, file = here::here("outputs/one_hot_11.3_clust_vect_lcm.rds"))

#dissimilarity matrix
gower_df <- daisy(df2,
                  metric = "gower" )
summary(gower_df)
dataset_dist <- stats::as.dist(gower_df)

#pcoa
dataset_pcoa <- ape::pcoa(dataset_dist)

#Recalculate relative eigenvalues by removing negative eigenvalues as in Mouillot et al.  
ev_pcoa <- dataset_pcoa$values$Eigenvalues
ev_pcoa_g0 <- ev_pcoa[ev_pcoa>0]
rel_ev_pcoa_g0 <- ev_pcoa_g0/sum(ev_pcoa_g0)

#check names
rownames(res_with@data@dataContinuous@notNA)==rownames(dataset_pcoa$vectors)

#plot points on first two axes, coloured by cluster
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2, fill = as.factor(fitted(res_with)))) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.5,
    size=3,
    stroke = 0.5
  ) + 
  stat_ellipse(geom = "polygon",
               aes(fill =  as.factor(fitted(res_with))), 
               alpha = 0.25) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

ggsave("figures/one_hot_11.3_scatterplot_pcoa_lcm_k3_coloured_by_cluster.png",width = 12,height=10)
