 # Read in trait data ----
library(dplyr)
library(ggmosaic)
library(tibble)

#load formatted data
df<-readRDS(file = here::here("outputs/df_filt_trans.rds"))

# Correlation between Traits ----

dataset<-df

dataset_cor <- matrix(0, ncol(dataset), ncol(dataset))

for (i in 1:ncol(dataset)) {

  for (j in i:ncol(dataset)) {

    dataset_cor[i, j] <- stats::cor(
      x      = rank(dataset[ , i]),
      y      = rank(dataset[ , j]),
      method = "kendall"
    )
  }
}

dataset_cor[lower.tri(dataset_cor)] <- t(
  dataset_cor)[lower.tri(dataset_cor)]

diag(dataset_cor) <- NA

dataset_cor <- data.frame(
  mean_cor = mean(abs(dataset_cor), na.rm = TRUE),
  sd_cor   = sd(abs(dataset_cor),   na.rm = TRUE),
  max_cor  = max(abs(dataset_cor),  na.rm = TRUE),
  min_cor  = min(abs(dataset_cor),  na.rm = TRUE)
)

dataset_cor

#mosaic plot
ggplot(data = df) +
  geom_mosaic(aes(x = product(Mating,Pollination,Lifespan), fill=Mating), na.rm=TRUE) + 
  labs(x = "Measurement_1", y = "Measurement_2")

#plot pairs of variables (removing some with high cat count)
library(GGally)
#pdf("figures/ggpairs.pdf",height=30,width=30)
png("figures/ggpairs.png",width=4000,height=4000,res=100)
ggpairs(df,cardinality_threshold=16) 
dev.off()

#dissimilarity matrix calc
par(mfrow=c(1,1))
library(cluster)
gower_df <- daisy(df,
                  metric = "gower" )

summary(gower_df)

#introduced NAs - need to solve better
#gower_df[is.na(gower_df)]<-0

#visualize the distances between species
#Need R > 4.0
library(factoextra)
png("figures/proteus_fviz_dist.png")
fviz_dist(dist.obj = gower_df,
          order = TRUE, show_labels = F)
dev.off()
