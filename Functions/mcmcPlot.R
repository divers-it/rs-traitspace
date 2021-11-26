#Function by Hervé Sauquet (2016)

mcmcPlot <- function(charID, subfolder="", iter="", par, parname)
{
	#Open a PDF graphical device to plot circular tree and progress pie chart
		charIDpath <- paste(subfolder, charID, sep="")
		if (iter=="") {
			filename <- paste(charIDpath, "_", parname, "_dist.pdf", sep="")
		} else {
			filename <- paste(charIDpath, "_", parname, "_trace.pdf", sep="")
		}
		pdf(file=filename,width=11,height=8,useDingbats=FALSE)
		
	#Plot the parameter
		if (iter=="") {
			hist(par, col="blue", breaks=50, main="", xlab=parname)
		} else {
			plot(iter, par, type="p", pch=20, cex=0.5, col="blue", xlab="Iteration", ylab=parname)
		}

	#Close PDF device
		cat(sep="","\n","ASR plot output to ",filename,"\n\n")
		dev.off()
}
