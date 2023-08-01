rm(list=ls())

# Library
library(ggplot2)
library(dplyr)
library(forcats)
library(cluster)

#load data set
df <- readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

#new df for editing
df2<-df

#dissimilarity matrix calculation of matrix
gower_df2 <- daisy(df2,
                   metric = "gower")

#make into distance object
dataset_dist2 <- stats::as.dist(gower_df2)


#vector for order of columns dropped
cols_dropped<-matrix(nrow=length(colnames(df))-1,ncol=2)

for (i in 1:(length(colnames(df))-1)) {
  
  auc_changes <- vector()
  
  for (j in 1:length(colnames(df2))) {
    df3 <- df2[, -j]
    gower_df3 <- daisy(as.data.frame(df3), metric = "gower")
    dataset_dist3 <- stats::as.dist(gower_df3)
    Co_rank   <-
      coRanking::coranking(dataset_dist2, dataset_dist3, input_Xi = "dist")
    NX        <- coRanking::R_NX(Co_rank)
    AUC       <- coRanking::AUC_ln_K(NX)
    auc_changes[j] <- 1-AUC
    
  }
  
  names(auc_changes)<-colnames(df2)
  
  #remove auc value with least change
  df2 <- df2[, !colnames(df2) %in% names(sort(auc_changes)[1])]
  
  #store name of column dropped
  cols_dropped[i,1]<-names(sort(auc_changes)[1])
  cols_dropped[i,2]<-sort(auc_changes)[1]

  
}

cols_dropped<-as.data.frame(cols_dropped)
colnames(cols_dropped)<-c('trait','AUC_decrease')

#add row for last trait remaining
remaining_trait<-setdiff(colnames(df),cols_dropped$trait)
cols_dropped<-rbind(cols_dropped,c(remaining_trait,1.0))

cols_dropped$AUC_decrease<-as.numeric(cols_dropped$AUC_decrease)

cols_dropped

cols_dropped %>%
  mutate(trait = fct_reorder(trait, desc(AUC_decrease))) %>%
  ggplot( aes(x=trait, y=AUC_decrease)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  ylab("Cumulative AUC decrease relative to original trait distance matrix") +
  theme_bw()

#NOTE: WHAT DOES THIS WARNING MEAN?
#Warning messages:
#1: In rankmatrix(Xi, input = input_Xi, use = use) :
#  0 outside of diagonal in distance matrix


