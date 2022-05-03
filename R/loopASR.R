#Function by Hervé Sauquet (2016)

#Description of parameters:
	# bct: branch (edge) coloring threshold (minimum probability of most likely state required for coloring a branch with this state)
	# modeloption: sets which models to go through (depending on number of states, see function listMkmodels)
		# "comprehensive": ARD and ER + UNI01 and UNI10 (if binary) + SYM, ORD, and ORDSYM (if 3-state)
		# "simple": ARD and ER only
	# rootmode: sets which root prior option(s) to use with unequal-rate models (ARD, SYM, ORD, ORDSYM)
		# "both": systematically tests both the flat and equilibrium options, treating them effectively as two separate models in the AIC comparison
		# "flat": equal state probabilities (default for rayDISC and most R ASR functions)
		# "equilibrium": equilibrium frequencies implied by the model (default for Mesquite)
	# mapallmodels: whether to output PDF character reconstructions for all tested models or just the best-fit model
	# charlist: (optional) restricts list of characters to loop through; e.g., c(x,y) or c(x:y)
	# testmode: test mode will only fit the "ARD" and "ER" models and will write in the R console rather than log files

loopASR <- function(usrdata, methodoption="ML", modeloption="comprehensive", rootmode="both", bct=0.5, colorscheme="bright", mapallmodels=FALSE, addfantree=FALSE, writedata=FALSE, charlist="", testmode=TRUE, runIDprefix="", resformat="Excel", BayesTraitsfolder="./", graphpar)
{

	#Set character state colors
		usrcolors <- setColors(colorscheme)
	
	#Set up log and output files and objects
		uniquedatetime <- format(Sys.time(), "%Y-%m-%d_%H-%M")
		runID <- paste(runIDprefix,"_", uniquedatetime, sep="")
		runfolder <- paste("./", runID, "/", sep="")
		treefolder <- paste(runfolder, "TreesByCharacter/", sep="")
		dir.create(runfolder)
		dir.create(treefolder)
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
	
	#Loop through characters and perform a number of actions for each of them
		if (charlist=="") charlist <- c(1:nrow(usrdata$charnames))
		for (charindex in charlist)
			{

				#Designate character to analyze and read character-specific data
					char <- readCharacter(charindex,usrdata$charmatrix,usrdata$charnames,usrdata$treeUP$tip.label,usrcolors$co)

				if (char$nrstates > 1) {
				
					#Initiate output files
						charID <- paste(runIDprefix, "_char", char$charnr, sep="")
						char$charID <- charID
						if (resformat=="Excel") {
							#sheet <- createSheet(wb, sheetName=char$charnr) #method for creating one sheet per character (paused for now)
						} else {
							csv1filename <- paste(statsfolder, charID, "_summary.csv", sep="")
							csv2filename <- paste(statsfolder, charID, "_beststandard.csv", sep="")
						}
						detrestable <- data.frame(usrdata$mapclades)
						colnames(detrestable)[1] <- "Node"
						cat(sep="","\n\n\n",rep("*",100),"\n\n")
						cat(sep="","Analyzing character: ",char$charname,"\n\n")
						cat(sep="","This character has ",char$nrstates," states: ",toString(char$states),"\n\n")

					if (grepl("MP", methodoption)==TRUE) {

						ansAndmore <- asrAndmore(usrdata$tree,usrdata$treeUP,usrdata$charmatrix,charindex,"",usrdata$mapclades,usrdata$mapnodesUP,char,usrcolors$co,bct,method="MP")

						#Save summary results to output files
							detrestable = cbind(detrestable, ansAndmore$bestcladestates[-1])

						#Plot ASR from best standard model on tree in PDF format
							if (addfantree==TRUE) fabTree(shape="fan",usrdata$tree,ansAndmore=ansAndmore,char,usrcolors$co,char$subco,usrcolors$circo,usrdata$mapclades,usrdata$mapnodes,usrdata$mapnodesUP,,,charID,treefolder,method="MP")
							fabTree(shape="up",usrdata$treeUP,ansAndmore=ansAndmore,char,usrcolors$co,char$subco,usrcolors$circo,usrdata$mapclades,usrdata$mapnodes,usrdata$mapnodesUP,"",bestfit=TRUE,charID,treefolder,method="MP",graphpar)
					}
					
					if (grepl("ML", methodoption)==TRUE) {

						#Loop through models and perform ASR for each one of them
							mkmodels <- listMkmodels(char$nrstates, char$charname, modeloption, rootmode)
							if (testmode==TRUE) mkmodels <- c("ARD","ER")
							nrmodels <- length(mkmodels)
							nrstmodels <- nrmodels
							loglik <- c()
							AIC <- c()
							AICc <- c()
							npar <- c()
							lowestAICc <- 100000
							ratelist <- matrix(,nrow=0,ncol=char$nrstates*char$nrstates) #empty matrix to store transition rates
							for (i in 1:nrmodels)
								{
									mkmodel <- mkmodels[i]
									ansAndmore <- asrAndmore(usrdata$tree,usrdata$treeUP,usrdata$charmatrix,charindex,mkmodel,usrdata$mapclades,usrdata$mapnodesUP,char,usrcolors$co,bct,method="ML")
									loglik <- c(loglik,ansAndmore$ans$loglik)
									AIC <- c(AIC,ansAndmore$ans$AIC) #fetched from rayDISC output, but could easily be re-calculated (=2*loglik+2*npar)
									AICc <- c(AICc,ansAndmore$ans$AICc) #fetched from rayDISC output, but could easily be re-calculated (function of loglik, npar, and sample size = total number of taxa)
									npar <- c(npar,length(unique(as.vector(ansAndmore$ans$index.mat)))-1) #fetch the number of parameters in the model from the rate matrix
									if  ((ansAndmore$ans$AICc < lowestAICc) & (i <= nrstmodels)) { #find the best model by comparing AICc (exclude the last model, identified with rjMCMC, which will always be the best)
										beststandardmodel <- ansAndmore
										lowestAICc <- ansAndmore$ans$AICc
									}
									rates <- as.vector(t(ansAndmore$ans$solution)) #decomposed rate matrix
									ratelist <- rbind(ratelist,rates) #append matrix of rates
									if (mapallmodels==TRUE) fabTree(shape="up",usrdata$treeUP,ansAndmore,char,usrcolors$co,char$subco,usrcolors$circo,usrdata$mapclades,usrdata$mapnodes,usrdata$mapnodesUP,mlsumtable,bestfit=FALSE,charID,treefolder,method="ML",graphpar)
								}

						#Summarize results in main log file and output summary tables to csv files
							ratelist <- ratelist[,-char$dropdiag] #drop diagonal (NA) rates
							#Prepare calculated data for the summary table
								deltaAICc <- AICc - lowestAICc
								omegavar <- exp(-0.5*deltaAICc) #intermediate variable to calculate Akaike weights (omega)
								omega <- omegavar/sum(omegavar) #Akaike weight
								bestmodelindex <- which(mkmodels==beststandardmodel$mkmodel) #get best model index (could also be stored in the loop above)
							#Annotate best model with a significance level based on the Akaike weight (note: this depends entirely on the number of models tested!)
								siglevel <- "*" #initialize: default significance level (Akaike weight less than 0.5)
								if (omega[bestmodelindex] >= 0.5) siglevel <- "**" #moderately significant (cutoff arbitrary)
								if (omega[bestmodelindex] >= 0.95) siglevel <- "***" #very significant
								fmkmodels <- mkmodels
								fmkmodels[bestmodelindex] <- paste(fmkmodels[bestmodelindex],siglevel,sep="")
							#Format table columns one by one (because nr of digits to display depends on each column)
								fmkmodels <- format(fmkmodels,justify="left")
								floglik <- round(loglik,digits=2)
								fAIC <- round(AIC,digits=2)
								fAICc <- round(AICc,digits=2)
								fdeltaAICc <- round(deltaAICc,digits=2)
								fomega <- round(omega,digits=2)
								fratelist <- round(ratelist,digits=4)
							mlsumtable <- cbind(data.frame(fmkmodels), data.frame(floglik), data.frame(npar), data.frame(fAIC), data.frame(fAICc), data.frame(fdeltaAICc), data.frame(fomega), fratelist)
							colnames(mlsumtable) <- c("Model","LogL","Npar","AIC","AICc","DeltaAICc","w",char$ratecolnames)
							cat(sep="","Summary table of model comparison\n")
							print(mlsumtable, right=FALSE)
							if (char$nrstates>2) { #if more than two states, the rate columns take took much space to be printed on figure so we need to truncate the table
								simpmlsumtable <- mlsumtable[,1:9]
								simpmlsumtable[,9] <- "..."
								colnames(simpmlsumtable)[9] <- "..."
							} else {
								simpmlsumtable <- mlsumtable
							}
							cat(sep="","\n","The best standard model found by AICc is ",beststandardmodel$mkmodel,"\n")
							cat(sep="","\n","Most likely states at key nodes for best standard model (",beststandardmodel$mkmodel,")\n")
							prettybestcladestates <- beststandardmodel$bestcladestates
							prettybestcladestates[,1] <- format(prettybestcladestates[,1], justify="left")
							prettybestcladestates[,2:3] <- format(prettybestcladestates[,2:3], justify="centre")
							print(prettybestcladestates, right=FALSE)

						#Save summary results to output files
							detrestable = cbind(detrestable, beststandardmodel$bestcladestates[-1])

						#Plot ASR from best standard model on tree in PDF format
							if (addfantree==TRUE) fabTree(shape="fan",usrdata$tree,ansAndmore=beststandardmodel,char,usrcolors$co,char$subco,usrcolors$circo,usrdata$mapclades,usrdata$mapnodes,usrdata$mapnodesUP,,,charID,treefolder,method="ML")
							fabTree(shape="up",usrdata$treeUP,ansAndmore=beststandardmodel,char,usrcolors$co,char$subco,usrcolors$circo,usrdata$mapclades,usrdata$mapnodes,usrdata$mapnodesUP,simpmlsumtable,bestfit=TRUE,charID,treefolder,method="ML",graphpar)
						
					}

					if (grepl("rjMCMC", methodoption)==TRUE) {

						#Do reversible-jump MCMC ancestral state reconstruction (for now, read ad summarize results of an output file from a previously run BayesTraits analysis)
							ansAndmore <- asrAndmore(usrdata$tree,usrdata$treeUP,usrdata$charmatrix,charindex,"",usrdata$mapclades,usrdata$mapnodesUP,char,usrcolors$co,bct,method="rjMCMC",BayesTraitsfolder=BayesTraitsfolder)

						#Save summary results to output files
							detrestable = cbind(detrestable, ansAndmore$bestcladestates[-1])

					}

					#Calculate a confidence score based on cross-method comparison:
						#*** if and only if the three methods are congruent (same unique best state) and the lower bound of the 95% HPD of the rjMCMC probability of the best state is more than (or equal to) 0.95
						#** if ML and rjMCMC results are congruent and rjMCMC lower bound > 0.5 and either MP, ML, and rjMCMC are congruent or rjMCMC lower bound >= 0.95
						#* otherwise
						
						if (methodoption=="MPMLrjMCMC") {
							conf <- matrix(c(rep("*",nrow(detrestable))),ncol=1)
							#De-factorize relevant columns (for some reason, probably the way it was built, detrestable is entirely factorized)
								MPstate <- as.character(detrestable[,2])
								MLstate <- as.character(detrestable[,3])
								rjMCMCstate <- as.character(detrestable[,5])
								rjMCMC_hpdlb <- as.numeric(as.character(detrestable[,8]))
							for (i in 1:nrow(detrestable)) {
								if (rjMCMC_hpdlb[i] >= 0.95) {
									if (MLstate[i]==rjMCMCstate[i]) {
											if (MPstate[i]==MLstate[i]) {
												conf[i] <- "***"
											} else {
												conf[i] <- "**"
											}
									}
								} else {
									if ((rjMCMC_hpdlb[i] > 0.5) & (MLstate[i]==rjMCMCstate[i]) & (MPstate[i]==MLstate[i])) {
										conf[i] <- "**"
									}
								}
							}
							detrestable <- detrestable[,1:(ncol(detrestable)-3)]
							detrestable <- cbind(detrestable, conf)
							colnames(detrestable)[ncol(detrestable)] <- "Confidence"
						} else {
							if (grepl("rjMCMC", methodoption)==TRUE) { detrestable <- detrestable[,1:(ncol(detrestable)-3)] }
						}

					#Save summary results to output files and print to console

						newcol <- matrix(c(rep("",nrow(detrestable))),nrow=nrow(detrestable))
						cdetrestable <- cbind(newcol, detrestable)
						colnames(cdetrestable)[1] <- "Character"
						temprow <- matrix(c(rep("",length(cdetrestable))),nrow=1,ncol=length(cdetrestable))
						newemptyrow <- data.frame(temprow)
						colnames(newemptyrow) <- colnames(cdetrestable)
						temprow[1,1] <- char$charname
						newcharrow <- data.frame(temprow)
						colnames(newcharrow) <- colnames(cdetrestable)
						cdetrestable <- rbind(newcharrow, cdetrestable)

						if (resformat=="Excel") {
							if (charindex==charlist[1]) {
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
								addDataFrame(ansAndmore$rjsumtable$rjmodeltable, sheet, row.names=FALSE, startRow=lastrownr+4, startColumn=lastcolnr+2)
								addDataFrame(ansAndmore$rjsumtable$rjstatstable, sheet, row.names=TRUE, startRow=lastrownr+nrow(ansAndmore$rjsumtable$rjmodeltable)+6, startColumn=lastcolnr+2)
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
						print(rbind(newemptyrow, newcharrow, prettycdetrestable), right=FALSE)

						if (charindex==charlist[1]) {
							sumrestable <- rbind(newemptyrow, cdetrestable)
							prettysumrestable <- rbind(newemptyrow, prettycdetrestable)
						} else {
							sumrestable <- rbind(sumrestable, cdetrestable)
							prettysumrestable <- rbind(prettysumrestable, prettycdetrestable)
						}

					#Optional: save the data for re-use
						if (writedata==TRUE) save(list=ls(all=TRUE),file=paste(runfolder, charID, ".RData", sep=""))

				} else {
					cat(sep="","\n","WARNING: character is constant or pseudoconstant (only a single monomorphic state sampled in the data). ASR skipped.\n")
				}

			}

		#Close and write global output files
			cat(sep="","\n\n\n",rep("#",100),"\n\n")
			cat(sep="","\n","FINAL SUMMARY: reconstructed ancestral states for all characters\n\n")
			print(prettysumrestable, right=FALSE)
			if (resformat=="Excel") {
				if (resformat=="Excel") saveWorkbook(wb, file=resfilename)
				cat(sep="","\n","All summary results output to ",resfilename,"\n\n")
			} else {
				if (testmode==FALSE) write.csv(sumrestable,file=sumresfilename,row.names=FALSE)
				cat(sep="","\n","All summary results output to ",sumresfilename,"\n\n")
			}
			sink()
}