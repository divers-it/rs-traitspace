#Function by Hervé Sauquet (2016)
#Writes custom labels (clade names) into a tree, and saves the labelled tree as a new file (optional, e.g., for import in Mesquite)
#This is different from the native nodelabels() function of ape, which merely adds labels to a plotted tree

writeClades <- function(tree, mapclades, mapnodes, outfilename = "")
{
	nrtips <- length(tree$tip.label)
	usrnodelabels <- rep("",nrtips-1)
	usrnodelabels[mapnodes-nrtips] <- mapclades
	tree$node.label <- usrnodelabels
	if (outfilename != "") {write.tree(tree,file=outfilename)}
	return(tree)
}

#Optional: write leaf collections in a text file to be pasted in BayesTraits script
#	logfilename <- "BayesTraitsNodes.txt"
#	sink(logfilename)
#	for (i in 1:nrmapnodes)
#		{
#			leaves <- gsub(",","",toString(sort(tips(tree,mapnodes[i]))))
#			cat(sep=" ","AddNode",mapclades[i],leaves,"\n")
#		}
#	sink()