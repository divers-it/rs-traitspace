rm(list = ls())

#load packages
library(dplyr)
library(gridExtra)
library(Rtsne)
library(ggplot2)
library(cluster)
library(RColorBrewer)
library(ggrepel)
library(networkD3)
library(patchwork)

# load formatted data
df <- readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

# build distance matrix
gower_df <- daisy(df,
                  metric = "gower" )
summary(gower_df)

## Select number of clusters (k) ----

# Silhouette Width to select the optimal number of clusters
# The silhouette width is one of the very popular choices when it comes to selecting the optimal number of clusters. 
# It measures the similarity of each point to its cluster, and compares that to the similarity of the point with the closest neighboring cluster. 
# This metric ranges between -1 to 1, where a higher value implies better similarity of the points to their clusters. 
# Therefore, a higher value of the Silhouette Width is desirable. 
# We calculate this metric for a range of cluster numbers and find where it is maximized. 
# The following code shows the implementation in R:

#empty vector
silhouette <- c()

#run PAM with different values of K from 2-10 and calculate silhouette width
for(i in 2:10){
  pam_clusters <- pam(as.matrix(gower_df),
                      diss = TRUE,
                      k = i)
  silhouette <- c(silhouette ,pam_clusters$silinfo$avg.width)
}

### Plot silhouette width ----
plot(2:10, silhouette,
     xlab = "Clusters",
     ylab = "Silhouette Width")
lines(2:10, silhouette)

j <- 1

# min clusters
mic <- 6
mac <- 10

# make df with cluster membership for each value of k from 2-7
for(i in mic:mac){
  if(i == mic){
    pam.gower = pam(gower_df, diss = TRUE, k = i)
    pam_df<-data.frame(pam.gower$clustering)
    colnames(pam_df)[j] <- paste(i,"clusters",sep="")
  } else {
    pam.gower = pam(gower_df, diss = TRUE, k = i)
    pam_df<-cbind(pam_df,pam.gower$clustering)
    colnames(pam_df)[j] <- paste(i,"clusters",sep="")
  }
  
  j <- j + 1
}

pam_df

clust_df <-as.data.frame(pam_df)
rownames(clust_df)<-names(pam.gower$clustering)

####
## Identify robust groups ----
####
# Robust groups are those that consistently stay together as values of k change

# make data frame of combo frequencies
combos <- as.data.frame(table(clust_df))

# remove no existant combos
combos <- combos[combos$Freq > 0, ]

# order
combos <- combos[order(combos$Freq, decreasing = T), ]

# change to strings
combos <-data.frame(lapply(combos, as.character), stringsAsFactors = FALSE)
plot(combos$Freq,
     main = paste(colnames(combos)[1],"to",colnames(combos)[length(combos[1,])-1]))

# empty list
robust<-list()

# empty vector
robust_vect_pam<-rep(NA,length(rownames(df)))
names(robust_vect_pam)<-rownames(df)

combos

# loop through ordered table to extract robust groups
for(i in 1:length(combos$Freq[as.numeric(combos$Freq)>0])){
              foo<-rownames(clust_df[clust_df[, 1] == combos[i, 1] & 
                                     clust_df[, 2] == combos[i, 2] &
                                     clust_df[, 3] == combos[i, 3] &
                                    clust_df[, 4] == combos[i, 4] &
                                     clust_df[, 5] == combos[i, 5]# &
                                     #clust_df[, 6] == combos[i, 6] &
                                       #clust_df[, 7] == combos[i, 7] &
                                       #clust_df[, 8] == combos[i, 8] &
                                       #clust_df[, 9] == combos[i, 9]
                                     ,])
  
  robust[[i]]<-foo
  
  robust_vect_pam[foo]<-i
  
}

# robust groups
robust

# keep 80% of the species in robust clusters ; others are NA
sum = 0

for (i in 1:max(robust_vect_pam)) {
  
  if (sum > 0.8 * length(robust_vect_pam)) {
    robust_vect_pam[robust_vect_pam == i] <- NA
    
  }
  
  sum = sum + length(robust[[i]])
  
}

table(robust_vect_pam)
