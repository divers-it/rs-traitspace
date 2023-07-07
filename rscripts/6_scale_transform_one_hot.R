rm(list=ls())

#Load libraries
library(dplyr)

#load formatted data
df<-readRDS(file = here::here("outputs/df_filt_one_hot.rds"))

#numeric columns only
nums <- unlist(lapply(df, is.numeric))
facts <- unlist(lapply(df, is.factor))

str(df)

boxplot(df[ , nums])

#centring introduced negative numbers that cant be logged
boxplot(scale(df[ , nums],center = F))

#scale and combine
df2<-cbind(scale(df[ , nums],center = F),df[ , facts])

nums <- unlist(lapply(df2, is.numeric))
facts <- unlist(lapply(df2, is.factor))

#plot histograms of quantitative variables
pdf("figures/proteus_trait_hists_one_hot.pdf")
par(mfrow=c(3,ceiling(sum(nums==TRUE)/3)))
#look at hists
for(i in 1:length(df2)){
  if(nums[i]){
  hist(df2[,i],main=colnames(df2)[i])
}
}
dev.off()

#plot histograms of logged (log10) variables
pdf("figures/proteus_trait_hists_transformed_one_hot.pdf")
par(mfrow=c(3,ceiling(sum(nums==TRUE)/3)))

for(i in 1:length(df2)){
  if(nums[i]){
  hist(log(df2[,i]),main=colnames(df2)[i])
}
}
dev.off()

#do log transformations
#not logging ovaries outcrossing rate
for(i in 1:length(df2)){
  if(nums[i]){
	if(colnames(df2)[i] != "Outcrossing.rate" & colnames(df2)[i] != "Fusion.of.ovaries"){
  df2[,i]<-log(df2[,i])
}
}
}

# Save scaled and transformed dataset
saveRDS(df2, file = here::here("outputs/df_filt_trans_one_hot.rds"))
