rm(ls=list())
library(dplyr)

#load formatted data
df<-readRDS(file = here::here("outputs/df_filt.rds"))

#character to factor
df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)],
                                       as.factor)

#integer to numeric
df[sapply(df, is.integer)] <- lapply(df[sapply(df, is.integer)],
                                     as.numeric)

# Remove traits with too much NA ----
df <- df[ , (colSums(is.na(df)) < length(df[,1])*0.6)]
str(df)

# Remove line with too much NA ----
df <- df[(rowSums(is.na(df)) < length(df[1,])*0.5), ]

#remove reproductive traits
df2<-subset(df, select=-c(SexualSystem,Mating,FlowerSex))

#check structure
str(df2)

#scale numeric
nums <- unlist(lapply(df2, is.numeric))
facts <- unlist(lapply(df2, is.factor))
df2<-cbind(scale(df2[ , nums]),df2[ , facts])
str(df2)

par(mfrow=c(1,1))

#dissimilarity matrix calculation
library(cluster)
gower_df <- daisy(df2,
                  metric = "gower" )

summary(gower_df)

dataset_dist <- stats::as.dist(gower_df)
dataset_pcoa <- ape::pcoa(dataset_dist)

#relative eigenvalues
eig_df<-data.frame(c(1:9),dataset_pcoa$values$Relative_eig[1:9])
colnames(eig_df)<-c("axis","relative_eigenvalue")
eig_df$axis<-as.character(eig_df$axis)

ggplot(eig_df, aes(x=axis, y=relative_eigenvalue)) + 
  geom_bar(stat = "identity")
ggsave("figures/rel_eig_pcoa_no_reproductive.pdf")

#plot points on first two axes
p1 <- ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2, fill = as.factor(df$SexualSystem))) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.5,
    size=4,
    stroke = 0.5
  ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))

p2 <- ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2, fill = as.factor(df$Mating))) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.5,
    size=4,
    stroke = 0.5
  ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))

p3 <- ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2, fill = as.factor(df$FlowerSex))) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.5,
    size=4,
    stroke = 0.5
  ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))


###
# Combined
###
library(patchwork)

p1 / p2 / p3
ggsave("figures/scatter_pcoa_no_reproductive.pdf",
       width = 20,
       height = 40,
       units = 'cm')
