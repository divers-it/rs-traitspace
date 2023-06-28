par(mar=c(4,4,4,4))

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


#plot phylogeny and example trait
plot(phy, show.tip.label = FALSE
     ,type='fan'
     )
tiplabels(pch = 16, col = as.factor(df$pam_one_hot), cex = 0.75)
