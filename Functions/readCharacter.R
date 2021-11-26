#Function by Hervé Sauquet (2016)
#Reads character data, character name, and character state names for a given character (charindex) from a matrix of data (charmatrix) and character and character state names (charnames)
#The list of tips in correct order (tips) is required to prepare the list of polymorphic tips (for mapping purposes)
#The color vector (co) is required to be subset in the event that absent states are detected
#Returns a whole bunch of objects used in the ASR script:
#	charData is a vector of data (in numeric format, missing as NA) tagged with taxon names; it is used in the ASR and fabtree functions for tip branch and terminal node coloring (for mapping purposes) only, not as source data
#	charname is a single-element vector containing the fully spelled character name
#	states is a vector of character state names
#	nrstates is a single-element vector containing the number of character states present in the dataset (i.e., original number of character states minus number of absent states)
#	dropdiag is a vector of diagonal references used to extract estimated transition rates from the transition rate matrix output by the rayDISC function
#	ratecolnames is a vector of transition rate names used to produce meaningful column names in summary table outputs
#	poltips is a matrix of polymorphic tips (1st column: tip number; next columns: state pseudoprobabilities)
#	subco is a vector of colors subset from the original color vector provided
#Note that none of this code is required to perform the actual ASR optimization or even to map simple node pie charts with ancestral states onto a plotted tree

readCharacter <- function(charindex, charmatrix, charnames, tips, co)
{
	#Read character data from matrix
		mxcharData <- charmatrix[charindex+1]
		rownames(mxcharData) <- unlist(charmatrix[1])
	#Transform character data into an appropriate object for preparing edge colors and plotting monomorphic tip states
		#I have no idea why this conversion is necessary but this is the recipe that works
		#The downside is that all missing and polymorphic data alike are transformed into NA in the process
		charData <- c(apply(mxcharData,1,as.numeric)) #used only for plotting terminal states on tree
	#Read character, character number (original), and character state names
		charname <- charnames[charindex,1]
		a <- unlist(strsplit(charname,"[.]"))
		charnr <- a[1] #tree string left of "."
		states <- as.vector(charnames[charindex,])
		states <- states[-1]
		states <- states[states != ""]
		states <- states[!is.na(states)]
		nrstates <- length(states)
	#Prepare a matrix of "probabilities" for polymorphic tip data in order to plot them as pie charts
		listpoltax <- grep("&",mxcharData[,1])
		poltips <- matrix(data=0,nrow=length(listpoltax),ncol=nrstates+1)
		if (length(listpoltax)>0) {
			listpoltax <- rownames(mxcharData)[listpoltax]
			polmxcharData <- cbind(listpoltax,gsub("&","",mxcharData[listpoltax,]))
			for (i in 1:length(listpoltax))
				{
					nrst <- nchar(polmxcharData[i,2])
					for (j in 1:nrst)
						{
							poltips[i,as.numeric(substr(polmxcharData[i,2],j,j))+2] <- 1/nrst
						}
					if (listpoltax[i] %in% tips) poltips[i,1] <- which(tips==listpoltax[i])
				}
			colnames(poltips) <- c("tip",0:(nrstates-1))
			if (length(which(poltips[,1]==0))>0) poltips <- poltips[-which(poltips[,1]==0),] #deletes taxa with no match in tree
		}
	#Check for absent states
		datastates <- levels(factor(mxcharData[,1]))
		datastates <- datastates[datastates != "?"]
		datastates <- datastates[datastates != ""] #although empty cells have already been replaced with missing data above (after reading the matrix), apparently the levels of this column have conserved the empty string as a possible value
		if (length(listpoltax)>0) datastates <- datastates[-grep("&",datastates)]
		datastates[datastates == "A"] <- 10
		datastates[datastates == "B"] <- 11
		datastates[datastates == "C"] <- 12
		datastates[datastates == "D"] <- 13
		datastates[datastates == "E"] <- 14
		datastates <- as.numeric(datastates)
		if (length(datastates) != nrstates) { #absent states detected!
			states <- states[datastates+1]
			nrstates <- length(states)
			subco <- co[datastates+1]
			cat(sep="","\n","WARNING: some of the states listed for the character are not sampled in the dataset. ASR will be performed ignoring these states.\n")
			flush.console() #print immediately to console to follow progress (otherwise delayed until function or loop finished)
		} else {
			subco <- co
		}
	#Prepare nrstates-dependent objects
		dropdiag <- c(1)
		for (i in 2:nrstates)
			{
				dropdiag <- c(dropdiag,(nrstates+1)*(i-1)+1)
			}
		ratecolnames <- c()
		for (i in 1:nrstates)
			{
				for (j in 1:nrstates)
					{
						ratecolnames <- c(ratecolnames,paste("q",i-1,j-1,sep=""))
					}
			}
		ratecolnames <- ratecolnames[-dropdiag] #drop diagonal (NA) rates
	list(charData=charData,charname=charname,charnr=charnr,states=states,nrstates=nrstates,dropdiag=dropdiag,ratecolnames=ratecolnames,poltips=poltips,subco=subco)
}