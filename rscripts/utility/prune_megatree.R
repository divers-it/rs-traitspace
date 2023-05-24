library(ape)
library(V.PhyloMaker)
library(Taxonstand)

qry_div=read.csv("data/qryDiveRS_Data_2023-05-21.csv")
spec_list=as.array(row.names(table(qry_div$NTaxDat)))
tpl_list= TPL(spec_list)
spec_list=tpl_list[,c("Taxon","New.Genus","Family")]
maker_tree=phylo.maker(spec_list)
plot(maker_tree$scenario.3)
write.tree(maker_tree$scenario.3,file="outputs/pruned_tree.tre")
