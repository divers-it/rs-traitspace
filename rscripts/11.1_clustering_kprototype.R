library(dplyr)
library(clustMixType)
library(wesanderson)

#load formatted data
df<-readRDS(file = here::here("outputs/df_filt_trans.rds"))

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
# k = 4 or 6

#rerun with chosen value of k
kproto_out<-kproto(df2, k=4, lambda = NULL, iter.max = 1000, nstart = 10, na.rm = F)

#relationships of dataset properties to clusters


#removes asking for plot
source("R/myclprofiles.R")

png("figures/kproto_cluster_characteristics_quant.png",height = 2000, width=1500,res=200)
par(mfrow=c(3,2))
myclprofiles(kproto_out, df2[,1:6], col = wes_palette("Royal1", 4, type = "continuous"))
dev.off()

png("figures/kproto_cluster_characteristics_qual1.png",height = 2000, width=1500,res=200)
par(mfrow=c(3,2))
myclprofiles(kproto_out, df2[,7:12], col = wes_palette("Royal1", 4, type = "discrete"))
dev.off()

png("figures/kproto_cluster_characteristics_qual2.png",height = 2000, width=1500,res=200)
par(mfrow=c(3,2))
myclprofiles(kproto_out, df2[,13:18], col = wes_palette("Royal1", 4, type = "discrete"))
dev.off()

## --------- PCOA scatterplot with cluster annotation ---------

par(mfrow=c(1,1))

#convert to distance matrix
dataset_dist <- stats::as.dist(gower_df)

#run PCOA
dataset_pcoa <- ape::pcoa(dataset_dist)

#plot points on first two axes, coloured by cluster
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2, fill = as.factor(kproto_out$cluster))) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.5,
    size=3,
    stroke = 0.5
  ) + 
  stat_ellipse(geom = "polygon",
               aes(fill = as.factor(kproto_out$cluster)), 
               alpha = 0.25) +
  xlab(paste("Axis 1: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(dataset_pcoa$values$Relative_eig[2],2)))

ggsave("figures/pcoa_kproto_k4.png",width = 12,height=10)
