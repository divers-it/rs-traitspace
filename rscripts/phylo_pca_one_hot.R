library(phytools)

df<-read.csv("outputs/imputed_with_phylo.csv", row.names = 1,stringsAsFactors = T)
str(df)

#examine data
head(df[,c(1:6)])

#SCALE? Need to transform to avoid biasing dataset over binary data?
#for(i in 1:6){ df[,i]<-((df[,i] - mean(df[,i]))/(sd(df[,i])))} #own implementation
df[,1:6]<-scale(df[,1:6],center = TRUE, scale = TRUE)

#check transformationg
head(df[,c(1:6)])

#read in phylogenetic tree
phy<-read.tree("outputs/pruned_tree.tre")
#plot(phy,cex=0.5)

#in dataset but not in phylo
setdiff(rownames(df),phy$tip.label)

#in phylo but not in dataset
setdiff(phy$tip.label,(rownames(df)))

#drop tips not in dataset
phy<-drop.tip(phy,setdiff(phy$tip.label,(rownames(df))))

#make neat
phy<-ladderize(phy)

#adephylo pPCA
library(adephylo)

#combine dataset and phylo
vir_4d <- phylobase::phylo4d(phy,df) 

#run pPCA
vir_pca <- adephylo::ppca(vir_4d, center = TRUE, scale = TRUE, scannf = FALSE, nfposi = 2, method = "Abouheif")

#plot pPCA
plot(vir_pca)

#phytools pPCA
ppca<-phyl.pca(phy,df)
print(ppca)

#plot pPCA
plot(ppca)
biplot(ppca)

#library(ggbiplot)
#ggbiplot(ppca, obs.scale = 1, var.scale = 1) +
#  scale_color_discrete(name = '') +
#  theme(legend.direction = 'horizontal', legend.position = 'top')
