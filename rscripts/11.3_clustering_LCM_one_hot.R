library(VarSelLCM)

#load formatted data
df<-readRDS(file = here::here("outputs/df_filt_trans.rds"))

#load data with NAs imputed
df2<-read.csv("outputs/imputed_with_phylo.csv",row.names = 1)

#remove biasing quantitative columns
df2<-subset(df2, select = -c(Number.of.fertile.stamens,Fusion.of.ovaries,Number.of.ovules.per.carpel))

# Cluster analysis without variable selection
res_without <- VarSelCluster(df2, gvals = 1:5, vbleSelec = FALSE, crit.varsel = "BIC")

# Cluster analysis with variable selection (with parallelisation)
res_with <- VarSelCluster(df2, gvals = 1:5, crit.varsel = "BIC")

#Comparison of the BIC for both models: variable selection permits to improve the BIC
BIC(res_without)
BIC(res_with)

# Estimated partition
fitted(res_with)

# Estimated probabilities of classification
head(fitted(res_with, type="probability"))

# Summary of the best model
summary(res_with)

# Discriminative power of the variables. The greater this index, the more the variable distinguishes the clusters.
plot(res_with)

# Summary of categorical variable
plot(res_with, y="Woodiness_herbaceous")

# More detailed output
print(res_with)

# Print model parameter
coef(res_with)

# Probabilities of classification for new observations 
# predict(res_with, newdata = x[1:3,])

# Imputation by posterior mean for the first observation
#not.imputed <- x[1,]
#imputed <- VarSelImputation(res_with, x[1,], method = "sampling")
#rbind(not.imputed, imputed)

#dissimilarity matrix calc - weights?
library(cluster)
gower_df <- daisy(df2,
                  metric = "gower" )

summary(gower_df)

dataset_dist <- stats::as.dist(gower_df)

dataset_pcoa <- ape::pcoa(dataset_dist)


#plot points on first two axes, coloured by cluster
library(ggplot2)

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
  xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))

ggsave("figures/pcoa_LCM_k4.png",width = 12,height=10)
