library(dplyr)

#load formatted data
df<-readRDS(file = here::here("outputs/df_filt.rds"))

#limit missing data
library(visdat)
vis_miss(df)

#check structure
str(df)

#add names
rownames(df)<-df$X

#remove species column
df<-subset(df, select=-c(X))

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

#numeric columns only
nums <- unlist(lapply(df2, is.numeric))
facts <- unlist(lapply(df2, is.factor))

df2<-cbind(df2[ , nums],df2[ , facts])

str(df2)

#centre and scale
#df_nums<-scale(df_nums)

#transform?

pdf("figures/proteus_trait_hists.pdf")
par(mfrow=c(3,3))
#look at hists
for(i in 1:6){
  hist(df2[,i],main=colnames(df2)[i])
}
dev.off()


pdf("figures/proteus_trait_hists_transformed.pdf")
par(mfrow=c(3,3))
#look at log10 hists
for(i in 1:6){
  hist(log(df2[,i]),main=colnames(df2)[i])
}
dev.off()

#do log transformations
for(i in 1:6){
  df2[,i]<-log(df2[,i])
  #df2[,i]<-scale(df2[,i]) #if we want to scale
}


par(mfrow=c(1,1))

#dissimilarity matrix calc - weights?
library(cluster)
gower_df <- daisy(df2,
                  metric = "gower" )

summary(gower_df)

dataset_dist <- stats::as.dist(gower_df)

dataset_pcoa <- ape::pcoa(dataset_dist)

plot(dataset_pcoa$vectors[,1]~dataset_pcoa$vectors[,2],col=as.factor(df$SexualSystem))


