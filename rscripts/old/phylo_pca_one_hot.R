rm(list = ls())

library(phytools)
library(ggplot2)
library(ggrepel)
library(patchwork)

df<-read.csv("outputs/one_hot_imputed_with_phylo.csv", row.names = 1,stringsAsFactors = T)
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

screeplot(vir_pca)

ggplot() + 
  geom_point(data=vir_pca$li,aes(x=PC1,y=PC2),shape=21,fill="grey",alpha=0.4,size=6)

p1 <- ggplot() + 
  geom_point(data=vir_pca$li,aes(x=PC1,y=PC2,fill=as.factor(vir_pca$tab$Woodiness_herbaceous)),shape=21,alpha=0.4,size=6) + 
  xlim(min(vir_pca$li$PC1)-1,max(vir_pca$li$PC1)+1) +
  ylim(min(vir_pca$li$PC2)-1,max(vir_pca$li$PC2)+1) +
  scale_fill_discrete(labels=c('0', '1'),name="Herbaceous") +
  geom_segment(data=vir_pca$c1,
               aes(x=0,y=0,xend=PA1*10,yend=PA2*10), col="grey30",
               #aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2,col=signal),
               #aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2,col=log(transition_rates)),
               arrow=arrow(length = unit(0.4, "cm")),
               alpha=0.7,
               linewidth=0.85,
               lineend='round',
               linejoin='round')  +
  geom_text_repel(data=vir_pca$c1,
                  aes(x=PA1*10,y=PA2*10,label=rownames(vir_pca$c1)),size=4) +
  theme_bw() +
  theme(legend.position = c(0.1,0.9))

p2 <- ggplot() + 
  geom_point(data=vir_pca$li,aes(x=PC1,y=PC2,fill=as.factor(vir_pca$tab$FlowerSex_unisexual)),shape=21,alpha=0.4,size=6) + 
  xlim(min(vir_pca$li$PC1)-1,max(vir_pca$li$PC1)+1) +
  ylim(min(vir_pca$li$PC2)-1,max(vir_pca$li$PC2)+1) +
  scale_fill_discrete(labels=c('0', '1'), name="Unisexual flowers") +
  geom_segment(data=vir_pca$c1,
               aes(x=0,y=0,xend=PA1*10,yend=PA2*10), col="grey30",
               #aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2,col=signal),
               #aes(x=0,y=0,xend=Dim1/2,yend=Dim2/2,col=log(transition_rates)),
               arrow=arrow(length = unit(0.4, "cm")),
               alpha=0.7,
               linewidth=0.85,
               lineend='round',
               linejoin='round')  +
  geom_text_repel(data=vir_pca$c1,
                  aes(x=PA1*10,y=PA2*10,label=rownames(vir_pca$c1)),size=4) +
  theme_bw() +
  theme(legend.position = c(0.1,0.9))


p1 + p2 

ggsave("figures/scatterplots_phylo_pca_loadings.png",width=20,height=12)

#plot pPCA
plot(vir_pca)

#phytools pPCA
ppca<-phyl.pca(phy,df)
print(ppca)

#plot pPCA
plot(ppca)
biplot(ppca, var.axes=TRUE)

s<-data.frame(ppca$S)
s$PC1

plot(s$PC1,s$PC2)

library(ggbiplot)
ggbiplot(ppca, obs.scale = 1, var.scale = 1) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')
