rm(list = ls())
library(densityClust)
library(cluster)

#load formatted data
df <- readRDS(file = here::here("outputs/df_filt_trans.rds"))

#calc gower distance
gower_df <- daisy(df,
                  metric = "gower" )

summary(gower_df)

#
protClust <- densityClust(gower_df, gaussian=TRUE)

# Inspect clustering attributes to define thresholds
plot(protClust) 

#
protClust <- findClusters(protClust, rho = 15,delta = 0.15,verbose = TRUE,plot = TRUE)

plotMDS(protClust)


split(rownames(df), protClust$clusters)
