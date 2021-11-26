#Function by Hervé Sauquet (2016-2017)

readBayesTraits <- function(datafolder, char, burnin_frac=0.1, corr=FALSE)
{

	#List the files in the designated folder that match the expected pattern:
		charnr <- gsub(" / ","_",char$charnr)
		charfiles <- list.files(path=datafolder, pattern=charnr, ignore.case=TRUE)
		charfiles <- charfiles[grepl("out", charfiles)==FALSE]
		charfiles <- charfiles[grepl("par", charfiles)==FALSE]

		if (length(charfiles)>0) {

			#Read the BayesTraits raw output file and store its contents in a character vector (rawoutput)
				if (length(charfiles)==1) {
					charfilename <- charfiles
				} else {
					charfilename <- charfiles[1]
					cat("! Found more than one file matching the search pattern. Will use the first in the list.\n")
				}
				cat("Reading ", charfilename, "...\n\n", sep="")
				charfile <- file(paste(datafolder, charfilename, sep=""))
				rawoutput <- readLines(charfile)
				close(charfile)
			
			#Extract the actual MCMC log part of the file (i.e., discard the header)
				rawchainstart <- max(grep("Iteration", rawoutput))
				rawchain <- rawoutput[rawchainstart:length(rawoutput)]
			
			#Convert the tab-delimited rawchain into matrix
				rawchain <- strsplit(rawchain, "\t")
				rawchain <- matrix(unlist(rawchain), nrow=length(rawchain), byrow=TRUE)
			
			#Rename the matrix columns (using the first row), then delete the first row
				colnames(rawchain) <- rawchain[1,]
				rawchain <- rawchain[-1,]
				rawchain <- rawchain[-1,] #delete the first iteration (BayesTraitsV3 onwards: first = 0)
			
			#Extract and report on basic information about the raw chain
				sampfreq <- as.numeric(rawchain[1,1])
				nrsamples <- nrow(rawchain)
				nrpar <- ncol(rawchain)
				cat("Successfully read the MCMC sample, with the following detected features:\n")
				cat("\tsampling frequency = ", sampfreq, "\n", sep="")
				cat("\tnr of sampled iterations = ", format(nrsamples, big.mark=","), "\n", sep="")
				cat("\ttotal nr of iterations = ", format(nrsamples*sampfreq, big.mark=",", scientific=FALSE), "\n", sep="")
				cat("\tnr of reported statistics = ", nrpar, "\n\n", sep="")
			
			#Discard the burnin (using the value input for this function)
				burnin_nrsamp <- nrsamples * burnin_frac
				cleanchain <- rawchain[(burnin_nrsamp+1):nrsamples,]
				nrcleansamples <- nrow(cleanchain)
				cat("Discarded the first ", format(burnin_nrsamp, big.mark=",", scientific=FALSE), " samples (or ", format(burnin_nrsamp * sampfreq, big.mark=",", scientific=FALSE), " iterations, or ", paste(round(100*burnin_nrsamp/nrsamples, 2), "%", sep=""), ") of the chain\n", sep="")
				cat("\tnr of sampled iterations after burnin removal = ", format(nrcleansamples, big.mark=","), "\n\n", sep="")
			
			#Extract basic information about model visitation (assuming this is a rjMCMC analysis)
				#modelstrcol <- normModelstr(cleanchain[,7]) #BayesTraitsV2
				modelstrcol <- normModelstr(cleanchain[,6]) #BayesTraitsV3
				modelstring <- levels(modelstrcol)
				nrmodels <- length(modelstring)
				transmodels <- matrix(c("0 0 ","ER","0 1 ","ARD","0 Z ","UNI01","Z 0 ","UNI10"), ncol=2, byrow=TRUE)
				nrsampmodel <- NULL
				model <- NULL
				for (i in 1:nrmodels) {
					nrsampmodel <- c(nrsampmodel, length(modelstrcol[modelstrcol==modelstring[i]]))
					modelindex <- which(transmodels[,1]==modelstring[i])
					if (length(modelindex>0)) {
						model <- c(model, transmodels[modelindex,2])
					} else {
						model <- c(model, "")
					}
				}
				postprob <- round(nrsampmodel/nrcleansamples, 2)
				rjmodeltable <- data.frame(modelstring, model, nrsampmodel, postprob)
				names(rjmodeltable) <- c("Model_str","Model","NrSamples","PostProb")
				rjmodeltable <- rjmodeltable[order(-nrsampmodel),]
				rownames(rjmodeltable) <- NULL
				cat("Total number of unique models visited by the rjMCMC chain = ", nrmodels, "\n", sep="")
				if (corr==FALSE) {
					cat("(Total number of possible models for a ", char$nrstates, "-state character = ", bell(char$nrstates*(char$nrstates-1)+1)-1, ")\n\n", sep="")
				} else {
					cat("(Total number of possible models for a combined (4x4) matrix of two binary characters, without dual transitions = 21146)\n\n", sep="")
				}
				if (nrow(rjmodeltable) > 5) {
					rjmodeltable <- head(rjmodeltable, n=5)
					cat("Top five models visited by the rjMCMC chain:\n\n")
				} else {
					cat("Summary of models visited by the rjMCMC chain:\n\n")
				}
				print(rjmodeltable)

			#Extract information about dependent vs. independent model visit frequency (applicable only to combined characters)
				if (corr==TRUE) {
					#depcol <- cleanchain[,8] #BayesTraitsV2
					depcol <- cleanchain[,7] #BayesTraitsV3
					D <- length(depcol[depcol=="D"])
					I <- length(depcol[depcol=="I"])
					BFDI <- D / I / (21146-51) * log(51)
					#Translate BF into qualitative support (see Posada and Buckley 2004 and Pagel and Meade 2006)
						if (BFDI<3) {
							corrsupport <- "virtually no"
						} else {
							if (BFDI<12) {
								corrsupport <- "positive"
							} else {
								if (BFDI<150) {
									corrsupport <- "strong"
								} else {
									corrsupport <- "very strong"
								}
							}
						}
					cat("Out of ", nrcleansamples, " sampled models, ", D, " represent dependent evolution between the two characters and ", I, " represent independent evolution.\n\n", sep="")
					cat("The Bayes Factor comparing dependent to independent models is BFdi = ", BFDI, ", indicating ", corrsupport, " support for correlated evolution between the two characters (", char$charnr, ").\n\n", sep="")
					#cleanchain <- cleanchain[,-8] #BayesTraitsV2
					cleanchain <- cleanchain[,-7] #BayesTraitsV3
				}

			#Discard model string column and convert cleanchain to numeric
				#cleanchain <- cleanchain[,-7] #BayesTraitsV2
				cleanchain <- cleanchain[,-6] #BayesTraitsV3
				class(cleanchain) <- "numeric"
				if (corr==TRUE) {
					list(cleanchain=cleanchain, rjmodeltable=rjmodeltable, BFDI=BFDI)
				} else {
					list(cleanchain=cleanchain, rjmodeltable=rjmodeltable)
				}

		} else {
			cat("!! No files could be found in folder ", datafolder, " that match the pattern of ", char$charnr, " in their name. BayesTraits extraction skipped.\n", sep="")
		}

}