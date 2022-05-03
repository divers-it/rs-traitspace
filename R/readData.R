#Function by Hervé Sauquet (2016, 2019)

#Description of inputs:
	# datafolder: path to data subfolder
	# treefile: core tree file with branch length for analyses (chronogram or phylogram)
	# maptreefile: optional secondary tree file with transformed branch lengths (cladogram) for mapping purposes only
	# datafile: core data file (matrix of one or more discrete characters in csv format)
	# charnamefile: file with names of characters and character states, for display purposes only
	# ingroupdef: ingroup definition (using two or more tips); by default, this script will prune the outgroups; to keep all taxa, set as ""
	# cladedefsfile

readData <- function(datafolder, treefile, maptreefile="", datafile, charnamefile, ingroupdef, cladedefsfile, trimtotree=TRUE, stretchtotree=FALSE)
{

	#Prepare a summary table of data read

		dataread = c("Tree file", treefile)
		dataread = rbind(dataread, c("Data file", datafile))
		dataread = rbind(dataread, c("Character names file", charnamefile))
		dataread = rbind(dataread, c("Clade definition file", cladedefsfile))
		rownames(dataread)=NULL

	#Add data subfolder path to data file names
	
		treefile <- paste(datafolder,treefile,sep="")
		if (maptreefile != "") maptreefile <- paste(datafolder,maptreefile,sep="")
		datafile <- paste(datafolder,datafile,sep="")
		charnamefile <- paste(datafolder,charnamefile,sep="")
		cladedefsfile <- paste(datafolder,cladedefsfile,sep="")
	
	#Read collection of trees and prepare master tree
	
		multitree <- read.tree(treefile)
		targettree <- multitree #replace with multitree[[1]] if source tree file is true multitree
		tree <- ladderize(targettree, TRUE) #Ladderize target tree (for clarity) and store in single-tree object "tree"
		treeUP <- ladderize(targettree, FALSE) #Ladderize target tree (for clarity) and store in single-tree object "tree"
		dataread <- rbind(dataread, c("Nr of taxa in tree", length(tree$tip.label)))
		
		#(Optional) cladogram version of this tree (created with Mesquite) for the purpose of graphical mapping only (use FALSE to ladderize correctly if tree will be drawn upright):
			if (maptreefile != "") {
				maptree <- ladderize(read.tree(maptreefile), FALSE)
			} else {
				maptree <- ""
			}
		
		#Extract ingroup subtree (i.e., remove all outgroups)
			if (ingroupdef != "") {
				ingroupnode <- getMRCA(tree, ingroupdef)
				tree <- extract.clade(tree, ingroupnode, root.edge=0)
				treeUP <- extract.clade(treeUP, ingroupnode, root.edge=0)
				if (maptree != "") maptree <- extract.clade(maptree, ingroupnode, root.edge=0)
				dataread <- rbind(dataread, c("Nr of taxa in tree (outgroups pruned)", length(tree$tip.label)))
			}

	#Read and prepare matrix
	
		#The lines below gets the matrix exactly right in the format required by the rayDISC function, with the first column (not rownames) containing the names of taxa and the following columns the trait data:
			charmatrix <- read.csv(datafile,header=FALSE,sep=",")
			if (is.na(charmatrix[1,length(charmatrix)])) charmatrix <- charmatrix[-length(charmatrix)] #drop the last column (artificially read because of the extra comma at the end of each row, due to the way the Output to R function of PROTEUS is currently programmed)
			colnames(charmatrix) <- c("Taxon",c(1:(length(charmatrix)-1))) #rename columns (for aesthetic purposes only)
			dataread <- rbind(dataread, c("Nr of taxa in matrix", nrow(charmatrix)))
			dataread <- rbind(dataread, c("Nr of characters in matrix", length(charmatrix) - 1))
			nrchars <- length(charmatrix) - 1
			charmatrix[charmatrix==""] <- "?" #replace empty cells with missing data (this should never be necessary with data exported from PROTEUS except in rare circumstances where a state has been deleted after the last data optimization)
			if (trimtotree==TRUE) charmatrix <- charmatrix[charmatrix[,1] %in% tree$tip.label,] #deletes from matrix taxa not sampled in the tree
			if (stretchtotree==TRUE) { #adds missing taxa rows to matrix (required for MP function); this code should work even when there are no missing taxa to add
				missingtaxa <- sort(tree$tip.label[! (tree$tip.label %in% as.vector(charmatrix$Taxon))])
				missingdata <- rep("?", length(missingtaxa))
				missingdata <- matrix(ncol = nrchars, nrow = length(missingtaxa))
				missingdata[] <- "?"
				missingrows <- cbind(data.frame(missingtaxa), missingdata)
				colnames(missingrows)[1] <- "Taxon"
				charmatrix <- rbind(charmatrix, missingrows)
			}
		
		charnames <- as.matrix(read.csv(charnamefile,header=TRUE,sep=";"))

	#Designate nodes to magnify (where to show pie charts)
	
		cladedefs <- defineClades(tree,cladedefsfile)
		cladedefsUP <- defineClades(treeUP,cladedefsfile)
		mapclades <- cladedefs$mapclades
		mapnodes <- cladedefs$mapnodes #list for vertical or fan tree (ladderize TRUE, plot direction rightwards)
		mapnodesUP <- cladedefsUP$mapnodes #list for horizontal tree (ladderize FALSE, plot direction upwards)
		dataread <- rbind(dataread, c("Nr of clades to map", length(mapclades)))
	
	list(tree=tree, treeUP=treeUP, maptree=maptree, charmatrix=charmatrix, charnames=charnames, mapclades=mapclades, mapnodes=mapnodes, mapnodesUP=mapnodesUP, dataread=dataread)
}