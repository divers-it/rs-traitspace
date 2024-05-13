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

####
# ---- AUC ----
####

#empty list
list_auc_changes<-list()

for (i in 1:(length(colnames(df))-1)) {
  
  auc_changes <- vector()
  
  for (j in 1:length(colnames(df2))) {
    df3 <- df2[, -j]
    gower_df3 <- daisy(as.data.frame(df3), metric = "gower")
    dataset_dist3 <- stats::as.dist(gower_df3)
    
    
    #This method was used in Dimensionality (Mouillot et al. 2021)
    #Supposed to be used for measuring the effect of dimensionality reduction
    #e.g. distance matrix vs pcoa
    #but here comparing distance matrices with variables dropped
    #unsure whether this is entirely appropriate
    Co_rank   <-
      coRanking::coranking(dataset_dist2, dataset_dist3, input_Xi = "dist")
    NX        <- coRanking::R_NX(Co_rank)
    AUC       <- coRanking::AUC_ln_K(NX)
    auc_changes[j] <- 1-AUC
    
  }
  
  names(auc_changes)<-colnames(df2)
  
  list_auc_changes[[i]]<-auc_changes
  
  #remove auc value with least change
  df2 <- df2[, !colnames(df2) %in% names(sort(auc_changes)[1])]
  
  #store name of column dropped
  cols_dropped[i,1]<-names(sort(auc_changes)[1])
  cols_dropped[i,2]<-sort(auc_changes)[1]

  
}

#examine an AUC change comparison
list_auc_changes[[1]]

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

####
# ---- Correlation ----
####

#Has problems when only a few traits with NAs in distance matrices

#new df for editing
df2<-df

#dissimilarity matrix calculation of matrix
gower_df2 <- daisy(df2,
                   metric = "gower")

#make into distance object
dataset_dist2 <- stats::as.dist(gower_df2)

#vector for order of columns dropped
cols_dropped<-matrix(nrow=length(colnames(df))-1,ncol=2)

list_cors<-list()

for (i in 1:(length(colnames(df))-1)) {
  
  cors <- vector()
  
  for (j in 1:length(colnames(df2))) {
    df3 <- df2[, -j]
    gower_df3 <- daisy(as.data.frame(df3), metric = "gower")
    dataset_dist3 <- stats::as.dist(gower_df3)
    
    #same as cor if not using p val
    #cors[j]<-vegan::mantel(dataset_dist2,dataset_dist3)
    
    cors[j]<-cor(dataset_dist2,dataset_dist3)[1]
    
    }
  
  names(cors)<-colnames(df2)
  
  list_cors[[i]]<-cors
  
  #remove auc value with least change
  df2 <- df2[, !colnames(df2) %in% names(sort(cors,decreasing=TRUE)[1])]
  
  #store name of column dropped
  cols_dropped[i,1]<-names(sort(cors,decreasing=TRUE)[1])
  cols_dropped[i,2]<-sort(cors,decreasing=TRUE)[1]
  
}

cols_dropped<-na.omit(as.data.frame(cols_dropped))
colnames(cols_dropped)<-c('trait','Cor')
cols_dropped

#add row for last trait remaining
remaining_traits<-setdiff(colnames(df),cols_dropped$trait)

remaining_traits_table<-cbind(remaining_traits,rep(0,length(remaining_traits)))
colnames(remaining_traits_table)<-c('trait','Cor')

cols_dropped<-rbind(cols_dropped,remaining_traits_table)

cols_dropped$Cor<-as.numeric(cols_dropped$Cor)

cols_dropped %>%
  mutate(trait = fct_reorder(trait, rev(desc(Cor)))) %>%
  ggplot( aes(x=trait, y=Cor)) +
  geom_bar(stat="identity", fill="#f68060", alpha=.6, width=.4) +
  coord_flip() +
  xlab("") +
  ylab("Correlation when removed to original distance matrix") +
  theme_bw()

####
# ---- Redo PCoA ----
####

#load data set
df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

#dissimilarity matrix calculation
gower_df <- daisy(df,
                  metric = "gower" )
summary(gower_df)

#make into distance object
dataset_dist <- stats::as.dist(gower_df)

#run PCoA on distance matrix
dataset_pcoa <- ape::pcoa(dataset_dist)

#plot PCOA points on first two axes
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    color=1,
    shape=16,
    alpha=0.75,
    size=2,
    stroke = 0.5
  )

df2 <- subset(df, select = -c(seedMass,
                              Numberofovulesperfunctionalcarpel,
                              flowerSize,
                              Numberofstructuralcarpels,
                              Maximumverticalheight,
                              Numberoffertilestamens,
                              Fusionofovaries#,
                              #Woodiness,
                              #FlowerSex
))

# df2 <- subset(df, select = c(Woodiness,
#                              Pollination,
#                              FlowerSex,
#                              DispersalMode,
#                              Lifespan,
#                              Showiness,
#                              seedMass))

#dissimilarity matrix calculation
gower_df2 <- daisy(df2,
                  metric = "gower" )
summary(gower_df2)

#make into distance object
dataset_dist2 <- stats::as.dist(gower_df2)

#run PCoA on distance matrix
dataset_pcoa2 <- ape::pcoa(dataset_dist2)

#plot PCOA points on first two axes
ggplot(data.frame(dataset_pcoa2$vectors), aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    color=1,
    shape=16,
    alpha=0.75,
    size=2,
    stroke = 0.5
  )

