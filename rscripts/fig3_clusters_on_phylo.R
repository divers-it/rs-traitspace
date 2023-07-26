par(mar=c(1,1,1,1))

library(ape)

#read in clustering results
df<-read.csv("outputs/collate_clustering_results.csv",row.names=1)

#insert '_' into rownames to match phylo
rownames(df)<-gsub(" ","_",rownames(df))

#read in phylogenetic tree
phy<-read.tree("outputs/pruned_tree.tre")

#in dataset but not in phylo
setdiff(rownames(df),phy$tip.label)

#in phylo but not in dataset
setdiff(phy$tip.label,rownames(df))

#in phylo but not in dataset
setdiff(phy$tip.label,rownames(df))

#drop tips
phy<-drop.tip(phy,setdiff(phy$tip.label,rownames(df)))

#check order
df<-df[match(phy$tip.label,rownames(df)),]
phy$tip.label==rownames(df)

#colours
palette(brewer.pal(3,"Set1"))

png("figures/phylo_ward.png",width=750,height=750)
#plot phylogeny and example trait
plot(phy, show.tip.label = FALSE
     ,type='fan'
     )
tiplabels(pch = 16, col = as.factor(df$ward), cex = 1.2)
dev.off()
