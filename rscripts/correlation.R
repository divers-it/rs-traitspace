 # Read in trait data ----
library(dplyr)
library(ggmosaic)

#load formatted data
df<-readRDS(file = here::here("outputs/df_filt.rds"))

#numeric columns only
nums <- unlist(lapply(df, is.numeric))
facts <- unlist(lapply(df, is.factor))
  
df_nums <-df[ , nums]

df2<-cbind(df[ , nums],df[ , facts])

str(df2)

# Correlation between Traits ----

dataset<-df2

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

view(dataset_cor)

#mosaic plot
ggplot(data = df2) +
  geom_mosaic(aes(x = product(Woodiness,SexualSystem), fill=Woodiness), na.rm=TRUE) + 
  labs(x = "Measurement_1", y = "Measurement_2")

ggplot(data = df2) +
  geom_mosaic(aes(x = product(Woodiness,Lifespan), fill=Woodiness), na.rm=TRUE) + 
  labs(x = "Measurement_1", y = "Measurement_2")

ggplot(data = df2) +
  geom_mosaic(aes(x = product(Woodiness,Mating), fill=Woodiness), na.rm=TRUE) + 
  labs(x = "Measurement_1", y = "Measurement_2")

ggplot(data = df2) +
  geom_mosaic(aes(x = product(Woodiness,Pollination), fill=Woodiness), na.rm=TRUE) + 
  labs(x = "Measurement_1", y = "Measurement_2")

ggplot(data = df2) +
  geom_mosaic(aes(x = product(Mating,Pollination,Lifespan), fill=Mating), na.rm=TRUE) + 
  labs(x = "Measurement_1", y = "Measurement_2")

 
#plot all numeric variables against one another
pairs(log(df_nums))


#plot pairs of variables (removing some with high cat count)
library(GGally)

pdf("figures/ggpairs.pdf",height=30,width=30)
ggpairs(df2[,c(1:11)]) 
dev.off()
