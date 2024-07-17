rm(list=ls())

#Load libraries
library(dplyr)

#load formatted data
df<-readRDS(file = here::here("outputs/5_df_filt.rds"))

#make vectors to split numeric and factor columns
nums <- unlist(lapply(df, is.numeric))
facts <- unlist(lapply(df, is.factor))

#examine data distribution
boxplot(df[ , nums])

#combine to ensure correct order
df2<-cbind(df[ , nums],df[ , facts])

#plot histograms of quantitative variables
pdf("figures/6_histograms.pdf")
par(mfrow=c(3,3))
#look at hists
for(i in 1:7){
  hist(df2[,i],main=colnames(df2)[i])
}
dev.off()

#plot histograms of logged (log10) variables
pdf("figures/6_histograms_transformed.pdf")
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

# Logit for proportional data ----

hist(df2[,'Fusionofovaries'])
hist(car::logit(df2[,'Fusionofovaries']))

# NOTE: doesn't seem to do much to the distribution
# df2[,'Fusionofovaries'] <- car::logit(df2[,'Fusionofovaries'])

# scale and centre numeric traits ----

df2<-cbind(scale(df2[ , 1:7],center = T, scale = T),df[ , facts])

#plot histograms of logged (log10) variables
pdf("figures/6_histograms_transformed_scaled.pdf")
par(mfrow=c(3,3))

for(i in 1:7){
  hist(df2[,i],main=colnames(df2)[i])
}
dev.off()

# Save scaled and transformed dataset
saveRDS(df2, file = here::here("outputs/6_df_filt_trans.rds"))
