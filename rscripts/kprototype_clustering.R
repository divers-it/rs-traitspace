library(dplyr)
library(clustMixType)
library(wesanderson)

#load formatted data
df<-readRDS(file = here::here("outputs/df_filt.rds"))

#remove mating system
#df<-subset(df, select=-c(sexmorphs))

#numeric columns only
nums <- unlist(lapply(df, is.numeric))
facts <- unlist(lapply(df, is.factor))

df2<-cbind(df[ , nums],df[ , facts])

str(df2)

ss<-vector()

for(i in 2:15){
  kproto_out<-kproto(df2, k=i, lambda = NULL, iter.max = 1000, nstart = 10, na.rm = F)
  ss[i]<-kproto_out$tot.withinss
}

#plot total ss
plot(ss,type='b')
# k = 4 or 5

#rerun with chosen value of k
kproto_out<-kproto(df2, k=3, lambda = NULL, iter.max = 1000, nstart = 10, na.rm = F)

#relationships of dataset properties to clusters
par(mfrow=c(2,2))
clprofiles(kproto_out, df2, col = wes_palette("Royal1", 5, type = "continuous"))


## --------- PCOA scatterplot with cluster annotation ---------

#convert to distance matrix
dataset_dist <- stats::as.dist(gower_df)

#run PCOA
dataset_pcoa <- ape::pcoa(dataset_dist)

#plot points on first two axes, coloured by cluster
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.3, color = as.factor(kproto_out$cluster))) +
  geom_point() +
  stat_ellipse(geom = "polygon",
               aes(fill = as.factor(kproto_out$cluster)), 
               alpha = 0.25) +
  xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))

