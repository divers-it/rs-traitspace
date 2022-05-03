#Function by Hervé Sauquet (2017)

#Primer (set-up mode)
	# charindexlist = c(1,3)
	# charmatrix = usrdata$charmatrix
	# charnames = usrdata$charnames
	# tips = usrdata$treeUP$tip.label
	# co = usrcolors$co

	readcorrCharacter <- function(charindexlist, charmatrix, charnames, tips, co)
{
	#Read character data from matrix
		# mxcharData <- charmatrix[charindex+1]
		# rownames(mxcharData) <- unlist(charmatrix[1])
	#Transform character data into an appropriate object for preparing edge colors and plotting monomorphic tip states
		#I have no idea why this conversion is necessary but this is the recipe that works
		#The downside is that all missing and polymorphic data alike are transformed into NA in the process
		# charData <- c(apply(mxcharData,1,as.numeric)) #used only for plotting terminal states on tree
	#Read character, character number (original), and character state names
		charname1 <- charnames[charindexlist[1],1]
		charname2 <- charnames[charindexlist[2],1]
		corrcharname <- paste(charname1, charname2, sep=" / ")
		a1 <- unlist(strsplit(charname1,"[.]"))
		charnr1 <- a1[1] #tree string left of "."
		a2 <- unlist(strsplit(charname2,"[.]"))
		charnr2 <- a2[1] #tree string left of "."
		corrcharnr <- paste(charnr1, charnr2, sep=" / ")
		states1 <- as.vector(charnames[charindexlist[1],])
		states1 <- states1[-1]
		states1 <- states1[states1 != ""]
		states1 <- states1[!is.na(states1)]
		nrstates1 <- length(states1)
		states2 <- as.vector(charnames[charindexlist[2],])
		states2 <- states2[-1]
		states2 <- states2[states2 != ""]
		states2 <- states2[!is.na(states2)]
		nrstates2 <- length(states2)
		nrcorrstates <- nrstates1*nrstates2
		corrstates <- c()
		corrsymbols <- c()
		pagelindex <- c(1:(nrcorrstates))
		for (i in 1:nrstates1) {
			for (j in 1:nrstates2) {
				corrstates <- c(corrstates, paste(states1[i], states2[j], sep=" / "))
				corrsymbols <- c(corrsymbols, paste(i-1, j-1, sep=","))
			}
		}
	# #Prepare a matrix of "probabilities" for polymorphic tip data in order to plot them as pie charts
		# listpoltax <- grep("&",mxcharData[,1])
		# poltips <- matrix(data=0,nrow=length(listpoltax),ncol=nrstates+1)
		# if (length(listpoltax)>0) {
			# listpoltax <- rownames(mxcharData)[listpoltax]
			# polmxcharData <- cbind(listpoltax,gsub("&","",mxcharData[listpoltax,]))
			# for (i in 1:length(listpoltax))
				# {
					# nrst <- nchar(polmxcharData[i,2])
					# for (j in 1:nrst)
						# {
							# poltips[i,as.numeric(substr(polmxcharData[i,2],j,j))+2] <- 1/nrst
						# }
					# if (listpoltax[i] %in% tips) poltips[i,1] <- which(tips==listpoltax[i])
				# }
			# colnames(poltips) <- c("tip",0:(nrstates-1))
			# if (length(which(poltips[,1]==0))>0) poltips <- poltips[-which(poltips[,1]==0),] #deletes taxa with no match in tree
		# }
	# #Check for absent states
		# datastates <- levels(factor(mxcharData[,1]))
		# datastates <- datastates[datastates != "?"]
		# datastates <- datastates[datastates != ""] #although empty cells have already been replaced with missing data above (after reading the matrix), apparently the levels of this column have conserved the empty string as a possible value
		# if (length(listpoltax)>0) datastates <- datastates[-grep("&",datastates)]
		# datastates <- as.numeric(datastates)
		# if (length(datastates) != nrstates) { #absent states detected!
			# states <- states[datastates+1]
			# nrstates <- length(states)
			# subco <- co[datastates+1]
			# cat(sep="","\n","WARNING: some of the states listed for the character are not sampled in the dataset. ASR will be performed ignoring these states.\n")
			# flush.console() #print immediately to console to follow progress (otherwise delayed until function or loop finished)
		# } else {
			# subco <- co
		# }
	# #Prepare nrstates-dependent objects
		# dropdiag <- c(1)
		# for (i in 2:nrstates)
			# {
				# dropdiag <- c(dropdiag,(nrstates+1)*(i-1)+1)
			# }
		# ratecolnames <- c()
		# for (i in 1:nrstates)
			# {
				# for (j in 1:nrstates)
					# {
						# ratecolnames <- c(ratecolnames,paste("q",i-1,j-1,sep=""))
					# }
			# }
		# ratecolnames <- ratecolnames[-dropdiag] #drop diagonal (NA) rates
	list(charname=corrcharname,charnr=corrcharnr,states=corrstates,nrstates=nrcorrstates,corrsymbols=corrsymbols)
}