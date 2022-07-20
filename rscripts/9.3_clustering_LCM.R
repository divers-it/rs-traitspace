library(VarSelLCM)

#load formatted data
df<-readRDS(file = here::here("outputs/df_filt_trans.rds"))

# Remove traits with too much NA ----
df <- df[ , (colSums(is.na(df)) < length(df[,1])*0.6)]
str(df)

# Remove line with too much NA ----
df <- df[(rowSums(is.na(df)) < length(df[1,])*0.5), ]

#remove mating system
#df<-subset(df, select=-c(sexmorphs))

#numeric columns only
nums <- unlist(lapply(df, is.numeric))
facts <- unlist(lapply(df, is.factor))

df_nums<-scale(df[ , nums])

df2<-cbind(df_nums,df[ , facts])

# Cluster analysis without variable selection
res_without <- VarSelCluster(df2, gvals = 1:4, vbleSelec = FALSE, crit.varsel = "BIC")

# Cluster analysis with variable selection (with parallelisation)
res_with <- VarSelCluster(df2, gvals = 1:4, crit.varsel = "BIC")

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

# Boxplot for the continuous variable MaxHeartRate
plot(x=res_with, y="Numberoffertilestamens")

#Empirical and theoretical distributions of the most discriminative variable (to check that the distribution is well-fitted)
plot(res_with, y="Numberoffertilestamens", type="cdf")

# Summary of categorical variable
plot(res_with, y="Woodiness")

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


par(mfrow=c(1,1))
plot(dataset_pcoa$vectors[,1]~dataset_pcoa$vectors[,2],col=as.factor(fitted(res_with)),pch=16,cex=2)

