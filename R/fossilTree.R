#Function by Hervé Sauquet (2017)

#Important note: all of the settings below (page size, margins, justifications, cex values, offset values) are optimized for a given number of taxa,
#here 792 taxa as in eFLOWER. Adapting this script for other tree sizes will require adjustments of these parameters.

fossilTree <- function(shape="long", tree, fossiltaxon, nrchars, edgescores, colorgradient, mapclades, mapnodes, subfolder="", graphpar=c(80,8,20,0.05))
{

	#Open a PDF graphical device to plot tree
		filename <- paste(subfolder, "Phyloscan_Paleo-eFLOWER_", fossiltaxon, "_", shape, ".pdf", sep="")
		if (shape=="long") {
			pdf(file=filename,width=graphpar[1],height=graphpar[2],useDingbats=FALSE)
			par(oma=c(0,0,0,0),mar=c(1,graphpar[3],0,0),xpd=TRUE)
		} else {
			pdf(file=filename,width=11,height=8,useDingbats=FALSE)
			layout(matrix(c(0,0,2,2,0,0,2,2,rep(1,16)),4,6))
			par(oma=c(5,5,5,3),mar=c(0,0,0,0),xpd=TRUE)
		}
		
	#Plot the tree
		if (shape=="long") {
			plot(tree, type="phylogram", direction="upwards", root.edge=FALSE, edge.width=3, edge.color=colorgradient[edgescores-min(edgescores)+1], cex=0.5, label.offset=3.5)
		} else {
			plot(tree, type="fan", rotate.tree=90, root.edge=TRUE, edge.width=1.5, edge.color=colorgradient[edgescores-min(edgescores)+1], cex=0.2, label.offset=1)
		}

	#Plot clade names
		if (shape=="long") {
			nodelabels(text=mapclades, node=mapnodes, frame="none", adj=c(-0.1,1.6), cex=0.8)
		} else {
		}

	#Plot title
		printfossiltaxon <- gsub("X ", "+", gsub("_", " ", fossiltaxon))
		title <- paste("Paleo-eFLOWER phyloscanning using parsimony\n", printfossiltaxon, " (", nrchars, " characters)", sep="")
		subtitle <- format(Sys.time(), "%A, %d %b %Y (%H:%M)")
		if (shape=="long") {
			text(x=-65,y=160,labels=title,adj=0,cex=1.2,font=1)
		} else {
			# par(oma=c(5,2,5,2),mar=c(0,0,0,0))
			title(main=title, cex.main=2, outer=TRUE)
		}

	#Plot legend
		if (shape=="long") {
			par(new=T)
			par(oma=c(0,0,0,0),mar=c(12,12,12,370),xpd=TRUE)
		} else {
		}
		# hist(edgescores, breaks=seq(min(edgescores),max(edgescores)+1,1), labels=TRUE, col=colorgradient, main=NULL, )
		edgecount <- count(as.factor(edgescores))
		xnames <- as.character(edgecount$x)
		xnames[1] <- paste(edgecount$x[1], "\n (MP)", sep="")
		xnames[2] <- paste(edgecount$x[2], "\n (MP+1)", sep="")
		if (length(edgecount$x)>2) xnames[3] <- paste(edgecount$x[3], "\n (MP+2)", sep="")
		mybar <- barplot(edgecount$freq, names.arg=xnames, col=colorgradient, ylim=c(0,1000))
		# text(mybar, edgecount$freq+20 , paste("n = ",edgecount$freq,sep="") ,cex=1) 
		text(mybar, edgecount$freq+40, edgecount$freq, cex=2) 

	#Close PDF device
		cat(sep="","\n","Phyloscan plot output to ",filename,"\n\n")
		dev.off()
}