#Function by Hervé Sauquet (2016)
#(adapted from readData.R and defineClades.R)

#Prunes tips from a single tree (typically the MCC from a BEAST analysis, in Newick format) or a sample of trees (typically a posterior sample of trees from a BEAST analysis, in NEXUS format)

#At the time this function was written, the ingroup extraction code was superfluous because the ingroup had already been extracted from all trees in the posterior (using Mesquite). However, this could be handy in the future (to bypass Mesquite).
#Currently set to prune all tips from one or more designated clades, but could easily be modified to prune tips randomly or from a given list of tips.

pruneTrees <- function(datafolder, treefile, ingroupdef, cladedefsfile, prunecladelist, outfilename)
{

	#Add data subfolder path to data file names
		treefile <- paste(datafolder,treefile,sep="")
		cladedefsfile <- paste(datafolder,cladedefsfile,sep="")
	
	#Read (collection of) tree(s)
		if (grepl("MCC", treefile)==TRUE) {
			singletree <- TRUE
			tree <- read.tree(treefile)
		} else {
			singletree <- FALSE
			multitree <- read.nexus(treefile)
			nrtrees <- length(multitree)
			tree <- multitree[[1]]
		}
		
	#Extract ingroup subtree (i.e., remove all outgroups)
		if (ingroupdef != "") {
			ingroupnode <- getMRCA(tree, ingroupdef)
			tree <- extract.clade(tree, ingroupnode, root.edge=0)
		}

	#Collect the list of tips to prune
		cladedef <- read.csv(cladedefsfile,header=FALSE,sep=",")
		mapclades <- as.character(cladedef[,1])
		nrpruneclades <- length(prunecladelist)
		tiplist <- list()
		prunetiplist <- c()
		mapnodes <- vector()
		mapsizes <- vector()
		for (i in 1:nrpruneclades)
			{
				deftiplist <- as.character(unlist(cladedef[match(prunecladelist,mapclades)[i],-1]))
				deftiplist <- deftiplist[deftiplist!=""]
				mapnodes[i] <- getMRCA(tree, deftiplist)
				tiplist[[i]] <- tips(tree,mapnodes[i])
				prunetiplist <- c(prunetiplist,tiplist[[i]])
				mapsizes[i] <- length(tiplist[[i]])
			}

	#Prune tips, tree by tree, then write (collection of) pruned tree(s) in an output file
		if (singletree==TRUE) {
			prunedtree <- drop.tip(tree,prunetiplist)
			write.tree(prunedtree, file=outfilename)
		} else {
			prunedtree <- list()
			for (i in 1:nrtrees)
				{
					prunedtree[[i]] <- drop.tip(multitree[[i]],prunetiplist)
				}
			write.nexus(prunedtree, file=outfilename)
		}

}