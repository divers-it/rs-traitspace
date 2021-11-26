#Function by Hervé Sauquet (2016)

#Plots a tree and reconstructed ancestral states

#Important note: all of the settings below (page size, margins, justifications, cex values, offset values) are optimized for a given number of taxa,
#here 792 taxa as in eFLOWER. Adapting this script for other tree sizes will require adjustments of these parameters.

fabTree <- function(shape="fan", tree, ansAndmore, char, co, subco, circo, mapclades, mapnodes, mapnodesUP, sumtable, bestfit = FALSE, charID, subfolder="", method="ML", graphpar=c(80,8,20,0.05,0.025,3,3.5,-65))
{

	#Open a PDF graphical device to plot circular tree and progress pie chart
		charIDpath <- paste(subfolder, charID, sep="")
		if (shape=="fan") {
			filename <- paste(charIDpath, "_fantree_", method, ".pdf", sep="")
			pdf(file=filename,width=8,height=9,useDingbats=FALSE)
			par(oma=c(0,0,5,0),mar=c(0,0,0,0),xpd=TRUE)
		} else {
			if (bestfit==TRUE) {
				filename <- paste(charIDpath, "_fulltree_", method, ".pdf", sep="")
			} else {
				filename <- paste(charIDpath,"_", ansAndmore$mkmodel, ".pdf", sep="")
			}
			pdf(file=filename,width=graphpar[1],height=graphpar[2],useDingbats=FALSE)
			par(oma=c(0,0,0,0),mar=c(1,graphpar[3],0,0),xpd=TRUE)
		}
		
	#Plot the tree
		if (shape=="fan") {
#			plot(tree, type="fan", rotate.tree=90, root.edge=TRUE, edge.width=1.5, edge.color=ansAndmore$edgeco, cex=0.2, label.offset=1)
			plot(tree, type="fan", rotate.tree=90, root.edge=TRUE, edge.width=3, edge.color=ansAndmore$edgeco, cex=0.4, label.offset=2)
		} else {
#			plot(tree, type="phylogram", direction="upwards", root.edge=FALSE, edge.width=graphpar[6], edge.color=ansAndmore$edgecoUP, cex=0.5, label.offset=graphpar[7])
			plot(tree, type="phylogram", direction="upwards", root.edge=FALSE, edge.width=graphpar[6], edge.color=ansAndmore$edgecoUP, cex=0.65, label.offset=graphpar[7])
		}

	#Plot tip states
		tips <- tree$tip.label
		if (shape=="fan") {
#			tiplabels(pch=19, col=co[char$charData[tips]+1], cex=0.18)
			tiplabels(pch=19, col=co[char$charData[tips]+1], cex=0.6)
		} else {
			tiplabels(pch=21, col=circo[char$charData[tips]+1], bg=co[char$charData[tips]+1], adj=c(0.5,2.5), cex=1)
		}

	#Plot polymorphic tip states
		if (length(char$poltips)>0) {
			if (shape=="fan") {
				#for now, do not bother plotting polymorphic tip states on fan trees
			} else {
				tiplabels(tip=char$poltips[,1],pie=char$poltips[,2:(char$nrstates+1)], piecol=co, adj=c(0.5,2.5), cex=graphpar[5])
			}
		}

	#Plot node state probabilities
		if (shape=="fan") {
			nodelabels(node=mapnodes, pie=ansAndmore$mapstates, piecol=subco, cex=0.5) #plot selected nodes only
		} else {
			nodelabels(pie=ansAndmore$ancstates, piecol=subco, cex=graphpar[4]) #plot all nodes
		}

	#Plot clade names
		if (shape=="fan") {
		} else {
			nodelabels(text=mapclades, node=mapnodesUP, frame="none", adj=c(-0.1,1.6), cex=0.8)
		}

	#Plot legend
		if (shape=="fan") {
			legend("topleft",char$states,pt.bg=subco,pch=21,bty="n",cex=1,adj=0,pt.cex=2)
		} else {
#			legend(x=graphpar[8],y=140,char$states,pt.bg=subco,pch=21,bty="n",cex=1,adj=0,pt.cex=2)
#		  legend(x=graphpar[8],y=65,char$states,pt.bg=subco,pch=21,bty="n",cex=1,adj=0,pt.cex=2)
		  legend(x=graphpar[8],y=180,char$states,pt.bg=subco,pch=21,bty="n",cex=1,adj=0,pt.cex=2)
		}

	#Plot title
		if (method=="ML") {
			title <- paste("ML ancestral state reconstruction using rayDISC (R:corHMM)\n", char$charname,", ",ansAndmore$mkmodel," model",sep="")
		} else {
			title <- paste("MP ancestral state reconstruction using ancestral.pars\n", "(R:phangorn)\n", char$charname, ", ", ansAndmore$ans$nrsteps, " steps", sep="")
		}
		subtitle <- format(Sys.time(), "%A, %d %b %Y (%H:%M)")
		if (shape=="fan") {
			title(main=title, cex.main=1, sub=subtitle, cex.sub=0.5, outer=TRUE)
		} else {
#			text(x=graphpar[8],y=160,labels=title,adj=0,cex=1.2,font=1)
#		  text(x=graphpar[8],y=70,labels=title,adj=0,cex=1.2,font=1)
		  text(x=graphpar[8],y=210,labels=title,adj=0,cex=1.2,font=1)
#		  addtable2plot(x=-40,y=137,yjust=0,table=ansAndmore$bestcladestates,cex=1)
#		  addtable2plot(x=-45,y=55,yjust=0,table=ansAndmore$bestcladestates,cex=1)
		  addtable2plot(x=-35,y=137,yjust=0,table=ansAndmore$bestcladestates,cex=1)
		  if ((bestfit==TRUE) & (method=="ML")) {
#				addtable2plot(x=-65,y=0,table=sumtable,cex=1)
				addtable2plot(x=-40,y=0,table=sumtable,cex=1)
			}
		}

	#Close PDF device
		cat(sep="","\n","ASR plot output to ",filename,"\n\n")
		dev.off()
}