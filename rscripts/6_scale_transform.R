library(dplyr)

#load formatted data
df<-readRDS(file = here::here("outputs/df_filt.rds"))

#remove mating system
#df<-subset(df, select=-c(sexmorphs))

#numeric columns only
nums <- unlist(lapply(df, is.numeric))
facts <- unlist(lapply(df, is.factor))

str(df)

boxplot(df[ , nums])

#centring introduced negative numbers that cant be logged
boxplot(scale(df[ , nums],center = F))

#scale and combine
df2<-cbind(scale(df[ , nums],center = F),df[ , facts])

#plot histograms of quantitative variables
pdf("figures/proteus_trait_hists.pdf")
par(mfrow=c(3,3))
#look at hists
for(i in 1:6){
  hist(df2[,i],main=colnames(df2)[i])
}
dev.off()

#plot histograms of logged (log10) variables
pdf("figures/proteus_trait_hists_transformed.pdf")
par(mfrow=c(3,3))

for(i in 1:6){
  hist(log(df2[,i]),main=colnames(df2)[i])
}
dev.off()

#do log transformations
#not logging ovaries and others that dont work
#for(i in c(1,2,3,5,6)){
for(i in c(1,2,3,4,6)){
  df2[,i]<-log(df2[,i])
}

# Save scaled and transformed dataset
saveRDS(df2, file = here::here("outputs/df_filt_trans.rds"))
