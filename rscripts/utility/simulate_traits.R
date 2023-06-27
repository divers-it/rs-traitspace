library(phytools)

#read in phylogenetic tree
phy<-read.tree("outputs/pruned_tree.tre")

#load formatted data
df<-readRDS(file = here::here("outputs/df_filt_trans.rds"))

#http://www.phytools.org/static.help/fastBM.html
#simulate quantitative trait using BM model and bounded by values from real data
x<-fastBM(phy,bounds = c(min(na.omit(df$Maximumverticalheight)),max(na.omit(df$Maximumverticalheight)))) # Brownian motion
x

#Simulate categorical variable

#transition matrix
Q <- matrix(c(-0.001,0.001,
              0.001,-0.001),ncol=2)

M <- sim.history(phy,Q,anc='1')
plot(M)
