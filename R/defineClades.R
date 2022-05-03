#Function by Hervé Sauquet (2016)

#Define a number of key clades for mapping or synthesizing results in ASR analyses
#The clades are defined in a separate csv file, where the first column is the name of the clade, and the next columns are the clade definitions (i.e., two or more tips, of which the MRCA defines the clade)
#Returns a vector of clade names (mapclades) and a vector corresponding node numbers in the given tree (mapnodes)
#Note that these numbers are conditional on a particular presentation of the tree (i.e., node numbers change if tree ladderized differently)

defineClades <- function(tree, csvfilename, silent=TRUE)
{
	if (silent==FALSE) { cat("Reading clade definition file (", csvfilename, ")...\n\n", sep="") }
	cladedef <- read.csv(csvfilename,header=FALSE,sep=",")
	mapclades <- as.character(cladedef[,1])
	nrmapnodes <- length(mapclades)
	mapnodes <- vector()
	mapsizes <- vector()
	for (i in 1:nrmapnodes)
		{
			deftiplist <- as.character(unlist(cladedef[i,-1]))
			deftiplist <- deftiplist[deftiplist!=""]
			mrca <- getMRCA(tree, deftiplist)
			if (is.null(mrca)) { mrca <- 0 }
			mapnodes[i] <- mrca
			mapsizes[i] <- length(tips(tree,mapnodes[i]))
			if (silent==FALSE) { cat(mapclades[i]," defined as MRCA of ", length(deftiplist), " taxa (", paste(deftiplist, collapse=", "), ")\n", sep="") }
		}
	mapclades=mapclades[mapnodes!=0]
	mapnodes=mapnodes[mapnodes!=0]
	mapsizes=mapsizes[mapnodes!=0]
	list(mapclades=mapclades,mapnodes=mapnodes,mapsizes=mapsizes)
}