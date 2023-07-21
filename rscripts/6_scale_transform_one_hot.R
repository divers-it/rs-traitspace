rm(list=ls())

#Load libraries
library(dplyr)

#load formatted data
df<-readRDS(file = here::here("outputs/df_filt_one_hot.rds"))

#make vectors to split numeric and factor columns
nums <- unlist(lapply(df, is.numeric))
facts <- unlist(lapply(df, is.factor))

#examine data distribution
boxplot(df[ , nums])

#combine to ensure correct order
df2<-cbind(df[ , nums],df[ , facts])

#plot histograms of quantitative variables
pdf("figures/proteus_trait_hists_one_hot.pdf")
par(mfrow=c(3,3))
#look at hists
for(i in 1:7){
  hist(df2[,i],main=colnames(df2)[i])
}
dev.off()

#plot histograms of logged (log10) variables
pdf("figures/proteus_trait_hists_transformed_one_hot.pdf")
par(mfrow=c(3,3))

for(i in 1:7){
  hist(log(df2[,i]),main=colnames(df2)[i])
}
dev.off()

#do log transformations
#NOTE: not logging ovaries and others that don't work
for(i in c(1,2,3,4,6,7)){
  df2[,i]<-log(df2[,i])
}

#scale and centre numeric traits
df2<-cbind(scale(df2[ , 1:7],center = T, scale = T),df[ , facts])

#plot histograms of logged (log10) variables
pdf("figures/proteus_trait_hists_transformed_scaled_one_hot.pdf")
par(mfrow=c(3,3))

for(i in 1:7){
  hist(df2[,i],main=colnames(df2)[i])
}
dev.off()

# Save scaled and transformed dataset
saveRDS(df2, file = here::here("outputs/df_filt_trans_one_hot.rds"))
