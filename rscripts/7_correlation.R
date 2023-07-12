rm(list=ls())

#load libraries
library(dplyr)
library(ggmosaic)
library(tibble)
library(factoextra)
library(dplyr)
library(cluster)
library(GGally)
 
#load formatted data
dataset<-readRDS(file = here::here("outputs/df_filt_trans.rds"))

# Correlation between Traits ----

#make empty matrix for correlation results
dataset_cor <- matrix(0, ncol(dataset), ncol(dataset))

#do kendall rank correlation between all traits
for (i in 1:ncol(dataset)) {

  for (j in i:ncol(dataset)) {

    dataset_cor[i, j] <- stats::cor(
      x      = rank(dataset[ , i]),
      y      = rank(dataset[ , j]),
      method = "kendall"
    )
  }
}

#transpose matrix to lower triangle
dataset_cor[lower.tri(dataset_cor)] <- t(dataset_cor)[lower.tri(dataset_cor)]

#make diagonal NA
diag(dataset_cor) <- NA

#make df for labelling
dataset_cor_df<-as.data.frame(dataset_cor)
rownames(dataset_cor_df)<-colnames(dataset)
colnames(dataset_cor_df)<-colnames(dataset)

view(dataset_cor_df)

#correlation stats
dataset_cor_summ <- data.frame(
  mean_cor = mean(abs(dataset_cor), na.rm = TRUE),
  sd_cor   = sd(abs(dataset_cor),   na.rm = TRUE),
  max_cor  = max(abs(dataset_cor),  na.rm = TRUE),
  min_cor  = min(abs(dataset_cor),  na.rm = TRUE)
)

dataset_cor_summ


#Alternative correlation on numeric traits

#make vectors to split numeric and factor columns
nums <- unlist(lapply(dataset, is.numeric))

#examine data distribution
dataset_num<-(dataset[ , nums])

#make empty matrix for correlation results
dataset_cor_num <- matrix(0, ncol(dataset), ncol(dataset))

#do correlation between all traits
for (i in 1:ncol(dataset_num)) {
  
  for (j in i:ncol(dataset_num)) {
    
    dataset_cor_num[i, j] <- stats::cor(
      x      = rank(dataset_num[ , i]),
      y      = rank(dataset_num[ , j]),
      method = "pearson"
    )
  }
}

#transpose matrix to lower triangle
dataset_cor_num[lower.tri(dataset_cor_num)] <- t(dataset_cor_num)[lower.tri(dataset_cor_num)]

#make diagonal NA
diag(dataset_cor_num) <- NA

#correlation stats for numeric traits
dataset_cor_num_summ <- data.frame(
  mean_cor = mean(abs(dataset_cor_num), na.rm = TRUE),
  sd_cor   = sd(abs(dataset_cor_num),   na.rm = TRUE),
  max_cor  = max(abs(dataset_cor_num),  na.rm = TRUE),
  min_cor  = min(abs(dataset_cor_num),  na.rm = TRUE)
)
dataset_cor_num_summ

###
# Uncomment to run (takes quite a long time)
###

##plot pairs of variables (removing some with high cat count)
#png("figures/ggpairs.png",width=4000,height=4000,res=100)
#ggpairs(df,cardinality_threshold=16) 
#dev.off()

#dissimilarity matrix calc
par(mfrow=c(1,1))

gower_df <- daisy(dataset,
                  metric = "gower" )

summary(gower_df)

#visualize the distances between species
png("figures/distance_heatmap.png")
fviz_dist(dist.obj = gower_df,
          order = TRUE, show_labels = F)
dev.off()

#One-liner to look at frequency of trait combinations
combo_df <-
  dataset %>% group_by(SexualSystem, FlowerSex, .drop = FALSE) %>%
  summarize(count = n())
combo_df

#mosaic plot to reperesent combinations visually
ggplot(data = dataset) +
  geom_mosaic(aes(x = product(Mating,Pollination,Lifespan), fill=Mating), na.rm=TRUE) + 
  labs(x = "Measurement_1", y = "Measurement_2")

