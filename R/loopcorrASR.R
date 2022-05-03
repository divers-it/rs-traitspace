#Function by Hervé Sauquet (2017)

#Primer (set-up mode)
	# methodoption="MLrjMCMC"
	# modeloption="comprehensive"
	# rootmode="both"
	# bct=0.5
	# colorscheme="bright"
	# mapallmodels=FALSE
	# addfantree=FALSE
	# writedata=FALSE
		# # charlist=c(2,21)
		# charlist=c(3:5)
		# # charlist=c(1,2)
	 # # # charlist=c(1,3,4)
	# testmode=FALSE
	# runIDprefix=runIDprefix
	# resformat="Excel"
		# # BayesTraitsfolder=paste("./BayesTraitsOutput", "/BayesTraits_", series, suffix, "/", sep="")
		# BayesTraitsfolder=paste("./BayesTraits_", dsversion, series, suffix, "/", sep="") 
	# graphpar=c(80,8,20,0.05)
	# prepareBayesTraits=""

loopcorrASR <- function(usrdata, methodoption="ML", modeloption="comprehensive", rootmode="both", bct=0.5, colorscheme="bright", mapallmodels=FALSE, addfantree=FALSE, writedata=FALSE, charlist="", testmode=TRUE, runIDprefix="", resformat="Excel", BayesTraitsfolder="./", graphpar, prepareBayesTraits="")
{

	#Set character state colors
		usrcolors <- setColors(colorscheme)
	
	#Set up log and output files and objects
		uniquedatetime <- format(Sys.time(), "%Y-%m-%d_%H-%M")
		runID <- paste(runIDprefix,"_", uniquedatetime, sep="")
		runfolder <- paste("./", runID, "/", sep="")
		# treefolder <- paste(runfolder, "TreesByCharacter/", sep="")
		dir.create(runfolder)
		# dir.create(treefolder)
		logfilename <- paste(runfolder, runIDprefix, ".out", sep="")
		if (resformat=="Excel") {
			resfilename <- paste(runfolder, runIDprefix, ".xlsx", sep="")
			wb <- createWorkbook()
			sheet <- createSheet(wb, sheetName=runID)
		} else {
			statsfolder <- paste(runfolder, "StatsByCharacter/", sep="")
			dir.create(statsfolder)
			sumresfilename <- paste(runfolder, runIDprefix, ".csv", sep="")
		}

	#Start logging (to both console and log file)
		sink(file=logfilename, split=TRUE)
		print(capture.output(sessionInfo()), quote=FALSE)
		cat("\n\n")
		print(usrdata$dataread, quote=FALSE)
		cat("\n")
	
	#Loop through all pairs of characters and perform a number of actions for each of them
		if (charlist=="") charlist <- c(1:nrow(usrdata$charnames))
		nrchars <- length(charlist)
		corrmatrix <- matrix("",nrow=nrchars, ncol=nrchars)
		chartable <- data.frame(CharPair=character(), AbsentStates=integer(), BestFitMLmodel=character(), AkaikeWeightD=double(), BestrjMCMCmodel=character(), BFdi=double(), stringsAsFactors=FALSE)
		charpairindex <- 0
		BTdatafilenamelist <- c() #(only used with option prepareBayesTraits)
		for (i in 1:(nrchars-1)) {
			for (j in (i+1):nrchars) {

				#Designate character to analyze and read character-specific data
					charindex1 <- charlist[i]
					charindex2 <- charlist[j]
					charpair <- c(charindex1,charindex2)
					char1 <- readCharacter(charindex1,usrdata$charmatrix,usrdata$charnames,usrdata$treeUP$tip.label,usrcolors$co)
					char2 <- readCharacter(charindex2,usrdata$charmatrix,usrdata$charnames,usrdata$treeUP$tip.label,usrcolors$co)
					corrchar <- readcorrCharacter(charpair,usrdata$charmatrix,usrdata$charnames,usrdata$treeUP$tip.label,usrcolors$co)
					charpairindex <- charpairindex + 1
					chartable[charpairindex,1] <- corrchar$charname

				if (prepareBayesTraits=="") { #analysis mode (default)
				
					#Initiate output files
						corrcharID <- paste(runIDprefix, "_chars", char1$charnr, "-", char2$charnr, sep="")
						if (resformat=="Excel") {
							# sheet <- createSheet(wb, sheetName=char$charnr) #method for creating one sheet per character (paused for now)
						} else {
							csv1filename <- paste(statsfolder, corrcharID, "_summary.csv", sep="")
							csv2filename <- paste(statsfolder, corrcharID, "_beststandard.csv", sep="")
						}
						detrestable <- data.frame(usrdata$mapclades)
						colnames(detrestable)[1] <- "Node"
						cat(sep="","\n\n\n",rep("*",100),"\n\n")
						cat(sep="","Analyzing characters: ",corrchar$charname,"\n\n")

					#Prepare the data for analysis by rayDISC
						#Subset matrix to only include the two characters to analyze
							charpairmatrix <- cbind(data.frame(usrdata$charmatrix[,1]),data.frame(usrdata$charmatrix[charlist[i]+1]),data.frame(usrdata$charmatrix[charlist[j]+1]))
						#Defactorize character columns (required for replacements)
							charpairmatrix[,2] <- as.character(charpairmatrix[,2])
							charpairmatrix[,3] <- as.character(charpairmatrix[,3])
						#First, replace all polymorphisms with missing data (this could be made once for matrix as a whole, but would remove the info about polymorphisms vs. missing data reported for single-trait analyses)
							charpairmatrix[charpairmatrix[,2]=="0&1",2] <- "?"
							charpairmatrix[charpairmatrix[,3]=="0&1",3] <- "?"
						#Compile combined column
							charpairmatrix[,4] <- paste(charpairmatrix[,2], charpairmatrix[,3], sep="")
						#Translate combined column into single multistate character with states 1 to 4 (to follow standard symbols of combined characters)
							translatetable <- matrix(c("00", "1", "01", "2", "10", "3", "11", "4", "?0", "1&3", "?1", "2&4", "0?", "1&2", "1?", "3&4", "??", "?"), ncol=2, byrow=TRUE)
							charpairmatrix[,5] <- translatetable[,2][match(charpairmatrix[,4],translatetable[,1])]
	
					#Summarize state frequency in combined character and report on absent states
						sumfreq <- count(as.factor(charpairmatrix[,5])) #note that rayDISC does this too and reports the same numbers
						colnames(sumfreq)[1] <- "state"
						cat(sep="","Summary of state frequencies in combined character:\n")
						print(sumfreq)
						missdata <- sumfreq[grep("\\?", sumfreq[,1]),2]
						poldata <- sum(sumfreq[grep("&", sumfreq[,1]),][,2])
						presdata <- sumfreq[grep("&", sumfreq[,1], invert=TRUE),]
						states <- c(1:4)
						translatetable <- matrix(c("00", "1", "01", "2", "10", "3", "11", "4", "?0", "1&3", "?1", "2&4", "0?", "1&2", "1?", "3&4", "??", "?"), ncol=2, byrow=TRUE)
						monostatefreqs <- sumfreq[,2][match(states,sumfreq[,1])]
						nrabsentstates <- sum(is.na(monostatefreqs))
						absentstatewarning <- ""
						if (nrabsentstates > 0) { #absent states detected
							absentstates_symbol <- paste(which(is.na(monostatefreqs)), collapse=", ")
							absentstates_verbal <- paste(paste("'", corrchar$states[which(is.na(monostatefreqs))], "'", sep=""), collapse=", ")
							if (nrabsentstates == 1) absentstatewarning <- "WARNING: 1 absent state detected" else absentstatewarning <- paste("WARNING: ", nrabsentstates, " absent states detected", sep="")
							absentstatewarning <- paste(absentstatewarning, " (", absentstates_symbol, " = ", absentstates_verbal, ")", sep="")
							cat(sep="","\n", absentstatewarning, "; absent states can produce strange results, especially with unconstrained models!\n\n")
						}
						chartable[charpairindex,2] <- nrabsentstates
						cat("\n")
	
					if (grepl("ML", methodoption)==TRUE) {

						#Loop through models and perform ASR for each one of them
							#Initialize
								# mkmodels <- c("ARD","ARDeq","ER")
								# nrmodels <- length(mkmodels)
								nrmodels <- 0
								# mkcorrmodels <- c("ARD","ARDeq","ARDnodual","ARDnodualeq","SYM","SYMeq","SYMnodual","SYMnodualeq")
								mkcorrmodels <- c("ARDnodual","ARDnodualeq","SYMnodual","SYMnodualeq")
								# mkcorrmodels <- c("SYMnodual")
								nrcorrmodels <- length(mkcorrmodels)
								# mkuncorrmodels <- c("ER","ERnodual","UNCORRnodual","UNCORRnodualeq")
								mkuncorrmodels <- c("ERnodual","UNCORRnodual","UNCORRnodualeq")
								# mkuncorrmodels <- c("ERnodual")
								nruncorrmodels <- length(mkuncorrmodels)
								nrtips <- length(usrdata$tree$tip.label)
								corrmodel <- c()
								corrmodeltype <- c()
								loglik <- c()
								npar <- c()
								bestcladestates <- NULL
							#Fit two independent models (uncorrelated evolution)
								# for (k in 1:nrmodels)
									# {
										# mkmodel <- mkmodels[k]
										# corrmodel <- c(corrmodel, paste(mkmodel, mkmodel, sep="x"))
										# corrmodeltype <- c(corrmodeltype,"independent")
										# ansAndmore1 <- asrAndmore(usrdata$tree,usrdata$treeUP,usrdata$charmatrix,charindex1,mkmodel,usrdata$mapclades,usrdata$mapnodesUP,char1,usrcolors$co,bct,method="ML")
										# ansAndmore2 <- asrAndmore(usrdata$tree,usrdata$treeUP,usrdata$charmatrix,charindex2,mkmodel,usrdata$mapclades,usrdata$mapnodesUP,char2,usrcolors$co,bct,method="ML")
										# loglik <- c(loglik,ansAndmore1$ans$loglik+ansAndmore2$ans$loglik)
										# npar <- c(npar,length(unique(as.vector(ansAndmore1$ans$index.mat)))-1+length(unique(as.vector(ansAndmore2$ans$index.mat)))-1) #fetch the number of parameters in the model from the rate matrix
										# bestcladestates[[k]] <- ansAndmore1$bestcladestates
									# }
							#Fit dependent models (correlated)
								for (k in 1:nrcorrmodels)
									{
										mkmodel <- mkcorrmodels[k]
										corrmodel <- c(corrmodel, mkmodel)
										corrmodeltype <- c(corrmodeltype,"D") #previously named 'dependent, correlated'
										corransAndmore <- corrasrAndmore(usrdata$tree,usrdata$treeUP,charpairmatrix,charpair,mkmodel,usrdata$mapclades,usrdata$mapnodesUP,corrchar,usrcolors$co,bct,method="ML")
										loglik <- c(loglik,corransAndmore$corrans$loglik)
										npar <- c(npar,length(unique(as.vector(corransAndmore$corrans$index.mat)))-1) #fetch the number of parameters in the model from the rate matrix
										bestcladestates[[nrmodels+k]] <- corransAndmore$bestcladestates
									}
							#Fit dependent models (uncorrelated)
								for (k in 1:nruncorrmodels)
									{
										mkmodel <- mkuncorrmodels[k]
										corrmodel <- c(corrmodel, mkmodel)
										corrmodeltype <- c(corrmodeltype,"I") #previously named 'dependent, uncorrelated'
										corransAndmore <- corrasrAndmore(usrdata$tree,usrdata$treeUP,charpairmatrix,charpair,mkmodel,usrdata$mapclades,usrdata$mapnodesUP,corrchar,usrcolors$co,bct,method="ML")
										loglik <- c(loglik,corransAndmore$corrans$loglik)
										npar <- c(npar,length(unique(as.vector(corransAndmore$corrans$index.mat)))-1) #fetch the number of parameters in the model from the rate matrix
										bestcladestates[[nrmodels+nrcorrmodels+k]] <- corransAndmore$bestcladestates
										# bestcladestates[[nrmodels+k]] <- corransAndmore$bestcladestates
									}

						#Summarize results in main log file and output summary tables to csv files
							#Prepare calculated data for the summary table
								AIC <- -2*loglik+2*npar #see Posada & Buckley 2004 (SystBiol)
								AICc <- AIC+(2*npar*(npar+1))/(nrtips-npar-1) #see Posada & Buckley 2004 (SystBiol)
								lowestAICc <- min(AICc)
								deltaAICc <- AICc - lowestAICc
								omegavar <- exp(-0.5*deltaAICc) #intermediate variable to calculate Akaike weights (omega)
								omega <- omegavar/sum(omegavar) #Akaike weight
								bestmodelindex <- which.min(AICc) #get best model index
							#Annotate best model with a significance level based on the Akaike weight (note: this depends entirely on the number of models tested!)
								siglevel <- "*" #initialize: default significance level (Akaike weight less than 0.5)
								if (omega[bestmodelindex] >= 0.5) siglevel <- "**" #moderately significant (cutoff arbitrary)
								if (omega[bestmodelindex] >= 0.95) siglevel <- "***" #very significant
								fmkmodels <- corrmodel
								fmkmodels[bestmodelindex] <- paste(fmkmodels[bestmodelindex],siglevel,sep="")
							#Format table columns one by one (because nr of digits to display depends on each column)
								fmkmodels <- format(fmkmodels,justify="left")
								fcorrmodeltype <- format(corrmodeltype,justify="left")
								floglik <- round(loglik,digits=2)
								fAIC <- round(AIC,digits=2)
								fAICc <- round(AICc,digits=2)
								fdeltaAICc <- round(deltaAICc,digits=2)
								fomega <- round(omega,digits=2)
							mlsumtable <- cbind(data.frame(fmkmodels), data.frame(fcorrmodeltype), data.frame(floglik), data.frame(npar), data.frame(fAIC), data.frame(fAICc), data.frame(fdeltaAICc), data.frame(fomega))
							colnames(mlsumtable) <- c("Model","Type","LogL","Npar","AIC","AICc","DeltaAICc","w")
							cat(sep="","Summary table of model comparison\n")
							print(mlsumtable, right=FALSE)
							cat(sep="","\n","The best standard model found by AICc is ",paste(corrmodel[bestmodelindex]," (", corrmodeltype[bestmodelindex], ")", sep=""),"\n")
							cat(sep="","\n","Most likely states at key nodes for best standard model\n")
							prettybestcladestates <- bestcladestates[[bestmodelindex]]
							prettybestcladestates[,1] <- format(prettybestcladestates[,1], justify="left")
							prettybestcladestates[,2:3] <- format(prettybestcladestates[,2:3], justify="centre")
							print(prettybestcladestates, right=FALSE)

						#Save summary results to output files
							detrestable = cbind(detrestable, bestcladestates[[bestmodelindex]][-1])

						#Plot ASR from best standard model on tree in PDF format
							# if (addfantree==TRUE) fabTree(shape="fan",usrdata$tree,ansAndmore=beststandardmodel,char,usrcolors$co,char$subco,usrcolors$circo,usrdata$mapclades,usrdata$mapnodes,usrdata$mapnodesUP,,,charID,treefolder,method="ML")
							# fabTree(shape="up",usrdata$treeUP,ansAndmore=beststandardmodel,char,usrcolors$co,char$subco,usrcolors$circo,usrdata$mapclades,usrdata$mapnodes,usrdata$mapnodesUP,simpmlsumtable,bestfit=TRUE,charID,treefolder,method="ML",graphpar)
						
					}

					if (grepl("rjMCMC", methodoption)==TRUE) {

						#Do reversible-jump MCMC ancestral state reconstruction (for now, read ad summarize results of an output file from a previously run BayesTraits analysis)
							corransAndmore <- corrasrAndmore(usrdata$tree,usrdata$treeUP,usrdata$charmatrix, charpair,"",usrdata$mapclades,usrdata$mapnodesUP,corrchar,usrcolors$co,bct,method="rjMCMC",BayesTraitsfolder=BayesTraitsfolder)

						#Save summary results to output files
							detrestable = cbind(detrestable, corransAndmore$bestcladestates[-1])
							detrestable <- detrestable[,1:(ncol(detrestable)-3)]

					}

					#Save summary results to output files and print to console

						newcol <- matrix(c(rep("",nrow(detrestable))),nrow=nrow(detrestable))
						cdetrestable <- cbind(newcol, detrestable)
						colnames(cdetrestable)[1] <- "Character"
						temprow <- matrix(c(rep("",length(cdetrestable))),nrow=1,ncol=length(cdetrestable))
						newemptyrow <- data.frame(temprow)
						colnames(newemptyrow) <- colnames(cdetrestable)
						temprow[1,1] <- corrchar$charname
						newcharrow <- data.frame(temprow)
						colnames(newcharrow) <- colnames(cdetrestable)
						cdetrestable <- rbind(newcharrow, cdetrestable)
						cdetrestable[] <- lapply(cdetrestable, as.character) #defactorize (required for next line)
						cdetrestable[2,1] <- absentstatewarning

						if (resformat=="Excel") {
							if ((i==1) & (j==2)) {
								addDataFrame(cdetrestable[0,], sheet, row.names=FALSE, startRow=1, startColumn=1)
							}
							lastrownr <- length(getRows(sheet))
							lastcolnr <- ncol(cdetrestable)
							addDataFrame(cdetrestable, sheet, col.names=FALSE, row.names=FALSE, startRow=lastrownr+1, startColumn=1)
							if (grepl("ML", methodoption)==TRUE) {
								freetitle <- "Summary of ML model comparisons"
								addDataFrame(freetitle, sheet, col.names=FALSE, row.names=FALSE, startRow=lastrownr+2, startColumn=lastcolnr+2)
								addDataFrame(mlsumtable, sheet, row.names=FALSE, startRow=lastrownr+4, startColumn=lastcolnr+2)
								lastcolnr <- lastcolnr + ncol(mlsumtable) + 1
							}
							if (grepl("rjMCMC", methodoption)==TRUE) {
								freetitle <- "Summary of rjMCMC results"
								addDataFrame(freetitle, sheet, col.names=FALSE, row.names=FALSE, startRow=lastrownr+2, startColumn=lastcolnr+2)
								addDataFrame(corransAndmore$rjsumtable$rjmodeltable, sheet, row.names=FALSE, startRow=lastrownr+4, startColumn=lastcolnr+2)
								addDataFrame(corransAndmore$rjsumtable$rjstatstable, sheet, row.names=TRUE, startRow=lastrownr+nrow(corransAndmore$rjsumtable$rjmodeltable)+6, startColumn=lastcolnr+2)
							}
						} else {
							write.csv(detrestable,file=csv2filename,row.names=FALSE)
							cat(sep="","\n","Summary results output to ",csv2filename,"\n\n")
							if (grepl("ML", methodoption)==TRUE) {
								write.csv(mlsumtable,file=csv1filename,row.names=FALSE)
								cat(sep="","\n","Summary table of model comparisons output to ",csv1filename,"\n\n")
							}
						}

						cat(sep="","\n","SUMMARY: reconstructed ancestral states for this character\n\n")
						prettycdetrestable <- cdetrestable
						prettycdetrestable[,2] <- format(prettycdetrestable[,2], justify="left")
						prettycdetrestable[,-(1:2)] <- format(prettycdetrestable[,-(1:2)], justify="centre")
						print(rbind(newemptyrow, prettycdetrestable), right=FALSE)

						if ((i==1) & (j==2)) {
							sumrestable <- rbind(newemptyrow, cdetrestable)
							prettysumrestable <- rbind(newemptyrow, prettycdetrestable)
						} else {
							sumrestable <- rbind(sumrestable, cdetrestable)
							prettysumrestable <- rbind(prettysumrestable, prettycdetrestable)
						}

					#Optional: save the data for re-use
						if (writedata==TRUE) save(list=ls(all=TRUE),file=paste(runfolder, charID, ".RData", sep=""))
					
					#Write summary of correlation test in the relevant cells of the correlation matrix and the character table
						if (grepl("ML", methodoption)==TRUE) {
							chartable[charpairindex,3] <- corrmodel[bestmodelindex]
							chartable[charpairindex,4] <- sum(mlsumtable[mlsumtable[,2]=="D",8])
							# if (grepl(" correlated",corrmodeltype[bestmodelindex])==TRUE) corrmatrix[i,j]="YES" else corrmatrix[i,j]="no"
							corrmatrix[i,j] <- sum(mlsumtable[mlsumtable[,2]=="D",8])
							# corrmatrix[j,i] <- ""
						}
						if (grepl("rjMCMC", methodoption)==TRUE) {
							chartable[charpairindex,5] <- as.character(corransAndmore$rjsumtable$rjmodeltable[1,1])
							chartable[charpairindex,6] <- round(corransAndmore$corrans$BFDI,digits=2)
							corrmatrix[j,i] <- round(corransAndmore$corrans$BFDI,digits=2)
						}

				} else { #BayesTraits data (paired traits) preparation mode

						#Subset matrix to only include the two characters to analyze
							charpairmatrix <- cbind(data.frame(usrdata$charmatrix[,1]),data.frame(usrdata$charmatrix[charpair[1]+1]),data.frame(usrdata$charmatrix[charpair[2]+1]))
						#Convert polymorphisms and missing data into BayesTraits format
							charpairmatrix[,2] <- gsub("&","",charpairmatrix[,2])
							charpairmatrix[,2] <- gsub("\\?","-",charpairmatrix[,2])
							charpairmatrix[,3] <- gsub("&","",charpairmatrix[,3])
							charpairmatrix[,3] <- gsub("\\?","-",charpairmatrix[,3])
						#Write or prepare output files
							filesuffix <- gsub(" / ","_",corrchar$charnr)
							BTdatafilename <- paste("bayestraitsData_",filesuffix,".txt",sep="")
							write.table(charpairmatrix,file=paste(runfolder,BTdatafilename,sep=""),quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)
							BTdatafilenamelist <- c(BTdatafilenamelist, BTdatafilename)

				}

			}
		}

		if (prepareBayesTraits=="") { #analysis mode (default)

			#Tidy up correlation matrix
				charnames <- substr(usrdata$charnames[,1],1,5)[charlist]
				diag(corrmatrix) <- c(rep(".", nrchars))
				corrmatrix <- data.frame(corrmatrix)
				rownames(corrmatrix) <- usrdata$charnames[,1][charlist]
				colnames(corrmatrix) <- charnames

			#Compute basic global stats from chartable
				nrpairs <- nrow(chartable)
				trimchartable <- chartable[chartable[,2] < 2,] #create a trimmed version of chartable excluding all pairs with two (or more!) absent states
				nrtrimpairs <- nrow(trimchartable)
				MLcorrsummary <- ""
				MLtrimcorrsummary <- ""
				rjMCMCcorrsummary <- ""
				rjMCMCtrimcorrsummary <- ""
				if (grepl("ML", methodoption)==TRUE) {
					nrpairsStrongD_ML <- nrow(chartable[chartable[,4] >= 0.95,])
					percpairsStrongD_ML <- paste(round(100*nrpairsStrongD_ML/nrpairs, 2), "%", sep="")
					nrtrimpairsStrongD_ML <- nrow(trimchartable[trimchartable[,4] >= 0.95,])
					perctrimpairsStrongD_ML <- paste(round(100*nrtrimpairsStrongD_ML/nrtrimpairs, 2), "%", sep="")
					MLcorrsummary <- paste("ML: Out of ", nrpairs, " character pairs tested, ", nrpairsStrongD_ML, " (", percpairsStrongD_ML, ") are strongly correlated (cumulative Akaike weight of dependent models >= 0.95)", sep="")
					MLtrimcorrsummary <- paste("ML (excl. pairs with 2+ absent states): Out of ", nrtrimpairs, " character pairs tested, ", nrtrimpairsStrongD_ML, " (", perctrimpairsStrongD_ML, ") are strongly correlated (cumulative Akaike weight of dependent models >= 0.95)", sep="")
					MLmodelfreq <- count(as.factor(chartable[,3]))
					colnames(MLmodelfreq)[1] <- "MLmodel"
				}
				if (grepl("rjMCMC", methodoption)==TRUE) {
					nrpairsStrongD_rjMCMC <- nrow(chartable[chartable[,6] >= 3,])
					percpairsStrongD_rjMCMC <- paste(round(100* nrpairsStrongD_rjMCMC/nrpairs, 2), "%", sep="")
					nrtrimpairsStrongD_rjMCMC <- nrow(trimchartable[trimchartable[,6] >= 3,])
					perctrimpairsStrongD_rjMCMC <- paste(round(100* nrtrimpairsStrongD_rjMCMC/nrtrimpairs, 2), "%", sep="")
					rjMCMCcorrsummary <- paste("rjMCMC: Out of ", nrpairs, " character pairs tested, ", nrpairsStrongD_rjMCMC, " (", percpairsStrongD_rjMCMC, ") are strongly correlated (BFdi >= 3)", sep="")
					rjMCMCtrimcorrsummary <- paste("rjMCMC (excl. pairs with 2+ absent states): Out of ", nrtrimpairs, " character pairs tested, ", nrtrimpairsStrongD_rjMCMC, " (", perctrimpairsStrongD_rjMCMC, ") are strongly correlated (BFdi >= 3)", sep="")
				}

			#Print final summaries to screen and output log file
				cat(sep="","\n\n\n",rep("#",100),"\n\n")
				cat(sep="","\n","FINAL SUMMARY (1): reconstructed ancestral states for all characters\n\n")
				print(prettysumrestable, right=FALSE)
				cat(sep="","\n","FINAL SUMMARY (2): matrix of all pairwise character correlations tested\n\n")
				print(corrmatrix)
				cat(sep="","\n","FINAL SUMMARY (3): table of all pairwise character correlations tested\n\n")
				print(chartable)
				cat(sep="","\n", MLcorrsummary, "\n", rjMCMCcorrsummary, "\n\n")
				cat(sep="","\n", MLtrimcorrsummary, "\n", rjMCMCtrimcorrsummary, "\n\n")
				if (grepl("ML", methodoption)==TRUE) {
					cat(sep="", "Best-fit ML model frequencies over all pairs tested:\n")
					print(MLmodelfreq)
				}

			#Print final summaries to Excel (or CSV) output file
				if (resformat=="Excel") {
					sheet <- createSheet(wb, sheetName="Correlation matrix")
					freetitle <- "MATRIX OF ALL PAIRWISE CHARACTER CORRELATIONS TESTED"
					addDataFrame(freetitle, sheet, col.names=FALSE, row.names=FALSE, startRow=1, startColumn=1)
					addDataFrame(corrmatrix, sheet, row.names=TRUE, startRow=3, startColumn=1)
					sheet <- createSheet(wb, sheetName="Character pair summary")
					freetitle <- "SUMMARY TABLE OF ALL PAIRWISE CHARACTER CORRELATIONS TESTED"
					addDataFrame(freetitle, sheet, col.names=FALSE, row.names=FALSE, startRow=1, startColumn=1)
					chartable[] <- lapply(chartable, as.character) #Excel does not like Inf values stored as double
					addDataFrame(chartable, sheet, col.names=TRUE, row.names=FALSE, startRow=3, startColumn=1)
					lastrownr <- length(getRows(sheet))
					addDataFrame(MLcorrsummary, sheet, col.names=FALSE, row.names=FALSE, startRow=lastrownr+3, startColumn=1)
					addDataFrame(rjMCMCcorrsummary, sheet, col.names=FALSE, row.names=FALSE, startRow=lastrownr+4, startColumn=1)
					addDataFrame(MLtrimcorrsummary, sheet, col.names=FALSE, row.names=FALSE, startRow=lastrownr+6, startColumn=1)
					addDataFrame(rjMCMCtrimcorrsummary, sheet, col.names=FALSE, row.names=FALSE, startRow=lastrownr+7, startColumn=1)
					saveWorkbook(wb, file=resfilename)
					cat(sep="","\n","All summary results output to ",resfilename,"\n\n")
				} else {
					if (testmode==FALSE) write.csv(sumrestable,file=sumresfilename,row.names=FALSE)
					cat(sep="","\n","All summary results output to ",sumresfilename,"\n\n")
				}

		} else { #BayesTraits data (paired traits) preparation mode
			runlines <- paste("./BayesTraitsV2_OpenMP eFLOWER_C_1706_ingroup.tre ", BTdatafilenamelist, " < rjMCMCcorr.txt", sep="")
			runlines <- c("#!/bin/sh","",runlines) #first two lines of BayesTraits run script (only used with option prepareBayesTraits)
			write.table(runlines,file=paste(runfolder,"rjCcorr.sh",sep=""),quote=FALSE,sep=" ",row.names=FALSE,col.names=FALSE)
		}

		sink()
}