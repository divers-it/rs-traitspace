#Function by Hervé Sauquet (2016)

normModelstr <- function(modelstrcol)
{
	modelstrcol <- factor(modelstrcol)
	modelstrlist <- levels(modelstrcol)
	nrmodels <- length(modelstrlist)
	newmodelstrlist <- NULL
	for (i in 1:nrmodels) {
		modelstr <- gsub("'", "", modelstrlist[i])
		modelpar <- unlist(strsplit(modelstr, " "))
		freemodelpar <- modelpar[modelpar!="Z"]
		par_levels <- levels(factor(freemodelpar))
		if (length(par_levels)>1) {
			par_inorder <- as.character(unique(freemodelpar))
			transmatrix <- matrix(c(par_inorder, par_levels), ncol=2)
			transmatrix <- transmatrix[order(transmatrix[,1]),]
			modelpar <- factor(modelpar)
			if ("Z" %in% modelpar) { addZ <- "Z" } else { addZ <- NULL }
			levels(modelpar) <- c(transmatrix[,2], addZ)
			newmodelstr <- paste(paste(modelpar, collapse=" "), " ", sep="")
		} else {
			newmodelstr <- modelstr
		}
		newmodelstrlist <- c(newmodelstrlist, newmodelstr)
	}
	levels(modelstrcol) <- newmodelstrlist
	return(modelstrcol)
}