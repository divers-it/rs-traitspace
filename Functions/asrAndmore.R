#Function by Hervé Sauquet (2016)
#Reconstructs ancestral states and summarizes the results in several useful ways for future extraction and graphical mapping
#Ancestral state reconstruction (ASR) is performed here with ML using the rayDISC function of corHMM for a discrete character, given:
#	- a tree (treeUP)
#	- a character matrix (charmatrix) and character number (charindex)
#	- a base model and root option (mkmodel)
#The results are then summarized as follows:
#	- ancestral states for key focal nodes (defined by mapclades and mapnodesUP) are extracted
#	- edge (branch) colors are set according to the most likely state or the tip state (sourced from char$charData) they subtend, given a minimum threshold of probability (bct) [note that this method if for intuitive graphical display only; coloring branches with a single color is in fact misleading because the models used allow multiple changes along a branch even if the start and end node have the same likeliest state]
#Here I use the "UP" version of the tree for ASR optimization, but prepare edge colors for both versions

asrAndmore <- function(tree, treeUP, charmatrix, charindex, mkmodel, mapclades, mapnodesUP, char, co, bct, method="ML", BayesTraitsfolder="./")
{

	nrtips <- length(tree$tip.label)
	rjsumtable <- ""

	if (method=="MP") {
		
		cat(sep="","Performing maximum parsimony ASR using ancestral.pars (character ", char$charnr, ")...\n\n")
		flush.console() #print immediately to console to follow progress (otherwise delayed until function or loop finished)
		
		#Prepare single-character matrix in the strange format required by phangorn
			#! some or all of this code could be moved to readCharacter.R
			#! here I subset the total matrix to the single character in focus, but ancestral.pars works just as fast on the total matrix (though I have not tried building a global contrast matrix yet)

				pmatrix <- charmatrix[charindex+1]
				pmatrix <- matrix(apply(pmatrix, 1, as.character))
				rownames(pmatrix) <- charmatrix[,1]
				
				#Build a custom so-called 'contrast' matrix to define polymorphic states
					#! this may see, very complicated, but there is no other way: phangorn does not seem to have a native system to read polymorphic data for the custom "USER" type (e.g., morphological), unlike corHMM
					#! the code below assumed that polymorphic data have been provided as requested by corHMM (e.g., 0/1 = "0&1")
					#! the code below will only work for characters with symbols in "0123456789" (more work needed for allowing non-numeric symbols)
					#! the alternative is to transform all polymorphic data into missing data, but this is not the exact same inference; using contrasts is much better and replicates exactly what Mesquite does

						contrast <- diag(char$nrstates)
						levels <- levels(factor(pmatrix))
						symbols <- levels(factor(as.numeric(levels)))
						dimnames(contrast) <- list(symbols, symbols)
						for (i in 1:length(levels)) {
							if (is.na(as.numeric(levels)[i])==TRUE) { #NA values introduced by this conversion correspond to polymorphic and missing data codes
								if (levels[i]=="?") {
									newcontrast <- rep(1,char$nrstates)
									contrast <- rbind(contrast, newcontrast)
								} else {
									atomicsymbols <- unlist(strsplit(levels[i], "&"))
									newcontrast <- rep(0,char$nrstates)
									contrast <- rbind(contrast, newcontrast)
									contrast[nrow(contrast), atomicsymbols] <- rep(1, length(atomicsymbols))
								}
								rownames(contrast)[nrow(contrast)]=levels[i]
							}
						}

				pdata <- phyDat(pmatrix, type="USER", contrast=contrast)

		#Reconstruct ancestral states with Fitch (unordered) parsimony using ancestral.pars (phangorn package)
			mpstates <- ancestral.pars(treeUP, pdata, type="MPR")
			nrsteps <- parsimony(treeUP, pdata, method="fitch")
			ans <- list(mpstates=mpstates, nrsteps=nrsteps)

		cat(sep="","Parsimony score: ", nrsteps, " steps\n\n")
		flush.console() #print immediately to console to follow progress (otherwise delayed until function or loop finished)
		ancstates <- matrix(unlist(ans$mpstates), nrow=length(ans$mpstates), byrow=T) #(as for nredges below, the total number of internal nodes should be nrtips*2-1, unless there are polytomies; using length(ans) is one way to get this number right)
		ancstates <- ancstates[-(1:nrtips),]

	}

	if (method=="ML") {

		cat(sep="","Performing maximum likelihood ASR (using rayDISC) with model ",mkmodel," (character ", char$charnr, ")...\n\n")
		flush.console() #print immediately to console to follow progress (otherwise delayed until function or loop finished)
		basemkmodel <- gsub("eq","",mkmodel)
	
		#Prepare rate matrix (not needed for ER and ARD models)
			if (! ((basemkmodel=="ER") | (basemkmodel=="ARD")) ) {
				rate.mat <- rate.mat.maker(rate.cat=1,hrm=FALSE,ntraits=1,nstates=char$nrstates,model="ARD") #prepares most general matrix for a single binary trait
				#IMPORTANT: models and rate matrices are defined below for a 2-state character
					if (basemkmodel=="UNI01") {
						rate.mat <- rate.par.drop(rate.mat,drop.par=c(1))
					}
					if (basemkmodel=="UNI10") {
						rate.mat <- rate.par.drop(rate.mat,drop.par=c(2))
					}
				#IMPORTANT: models and rate matrices are defined below for any multistate character (i.e., 3 or more states)
					if (grepl("ORD",basemkmodel)==TRUE) {
						for (i in 1:(char$nrstates-1)) {
							for (j in 1:(char$nrstates-2)) {
								rate.mat <- rate.par.drop(rate.mat,drop.par=c(i*2))
							}
						}
					}
					if (grepl("SYM",basemkmodel)==TRUE) { #this actually already exists natively in rayDISC (but looks like it does something different)
						for (i in 1:(char$nrstates*(char$nrstates-1)/2)) {
							rowi <- which(rate.mat==i,arr.ind=TRUE)[,1]
							coli <- which(rate.mat==i,arr.ind=TRUE)[,2]
							rate.mat <- rate.par.eq(rate.mat,eq.par=c(i,rate.mat[coli,rowi]))
						}
					}
					if (grepl("ER",basemkmodel)==TRUE) {
						rate.mat <- rate.mat/rate.mat
					}
				cat(sep="","Rate matrix\n")
				print(rate.mat)
				cat(sep="","\n")
				flush.console() #print immediately to console to follow progress (otherwise delayed until function or loop finished)
			}
	
		#Estimate transition rates and reconstruct ancestral states with rayDISC (corHMM package)
			if (grepl("eq",mkmodel)==TRUE) rootoption="maddfitz" else rootoption=NULL #specify the root prior option: by default, rayDISC uses a flat prior; when the suffix "eq" is added to the model name, a prior with equilibrium frequencies will be used instead (as is the default option in Mesquite); this concerns only models with unequal transition rates (e.g., ARD, SYM, ORD) and cannot be applied to the unidirectional models (for which the root state is instead determined by the model itself)
			if ((basemkmodel=="ER") | (basemkmodel=="ARD")) {
				ans <- rayDISC(treeUP,charmatrix,ntraits=1,charnum=charindex,model=basemkmodel,node.states="marginal",root.p=rootoption)
			} else {
				ans <- rayDISC(treeUP,charmatrix,ntraits=1,charnum=charindex,rate.mat=rate.mat,node.states="marginal",root.p=rootoption)
			}
	
		print(ans)
		ancstates <- ans$states
		flush.console() #print immediately to console to follow progress (otherwise delayed until function or loop finished)
		
	}

	if (method=="rjMCMC") {
	
		cat(sep="","Summarizing results from reversible-jump MCMC BayesTraits log file (character ", char$charnr, ")...\n\n")
		flush.console() #print immediately to console to follow progress (otherwise delayed until function or loop finished)
		
		#Read results from a previously run BayesTraits analysis (and summarize models visited)
			ans <- readBayesTraits(datafolder=BayesTraitsfolder, char=char)
			ancstates <- ""
		
		#Summarize basic statistics on number of parameters and transition rates (cleanchain column indices adjusted for BayesTraitsV3 outputs)
		
			#First, summarize number of free and zero parameters (for each model visited by the chain)
				meannrpar <- format(round(mean(ans$cleanchain[,4]), digits=1), scientific=FALSE) #(note that the first column in ans$cleanchain is teh Iteration column)
				hpdnrpar <- format(round(boa.hpd(ans$cleanchain[,4], 0.05), digits=1), scientific=FALSE)
				meannrzero <- format(round(mean(ans$cleanchain[,5]), digits=1), scientific=FALSE)
				hpdnrzero <- format(round(boa.hpd(ans$cleanchain[,5], 0.05), digits=1), scientific=FALSE)
				rjstatstable <- data.frame(rbind(cbind(meannrpar, meannrzero), cbind(hpdnrpar, hpdnrzero)))
			
			#Then, summarize transition rates
				for (i in 1:length(char$ratecolnames)) { #(note that length(char$ratecolnames) should be equal to char$nrstates * (char$nrstates - 1))
					meanq <- format(round(mean(ans$cleanchain[,5+i]), digits=4), scientific=FALSE)
					hpdq <- format(round(boa.hpd(ans$cleanchain[,5+i], 0.05), digits=4), scientific=FALSE)
					rjstatstable <- cbind(rjstatstable, matrix(c(meanq, hpdq), ncol=1))
					# mcmcPlot(char$charID, , iter=ans$cleanchain[,1], ans$cleanchain[,5+i], char$ratecolnames[i])
					# mcmcPlot(char$charID, , , ans$cleanchain[,5+i], char$ratecolnames[i])
				}
			
			#Tidy up and print
				colnames(rjstatstable) <- c("Npar","Nzero", char$ratecolnames)
				rownames(rjstatstable) <- c("mean","95% lb","95% ub")
				cat(sep="","\n","Summary statistics for number of parameters and transition rates\n")
				prettyrjstatstable <- format(rjstatstable, justify="centre")
				print(prettyrjstatstable, quote=FALSE, right=FALSE)
				rjsumtable <- list(rjmodeltable=ans$rjmodeltable, rjstatstable=rjstatstable)

	}
	
	#Map and summarize focal node states

		nrmapnodes <- length(mapclades)

		if ((method=="MP") | (method=="ML")) {
			mapstates <- ancstates #suboptimal (but works), I would do this more elegantly now by creating an empty matrix as below
			for (i in 1:nrmapnodes)
				{
					mapstates[i,] <- ancstates[mapnodesUP[i]-nrtips,]
				}
			allcladestates <- cbind(data.frame(mapclades), round(mapstates[1:nrmapnodes,],digits=4))

			colnames(allcladestates) <- c("Node",char$states)
			if (method=="MP") {
				bestcladestates <- matrix(,nrow=nrmapnodes,ncol=1) #empty matrix to store most parsimonious clade states
				for (i in 1:nrmapnodes)
					{
						favstate <- which(allcladestates[i,-1]>0) #find most parsimonious state indices
						bestcladestates[i,1] <- paste(char$states[favstate], sep=" / ", collapse=" / ")
					}
				bestcladestates <- cbind(data.frame(mapclades), data.frame(bestcladestates))
				colnames(bestcladestates) <- c("Node","MP state(s)")
				cat(sep="","Ancestral states for key nodes\n")
				prettyallcladestates <- allcladestates
				prettyallcladestates[,1] <- format(prettyallcladestates[,1], justify="left")
				prettyallcladestates[,2:char$nrstates+1] <- format(prettyallcladestates[,2:char$nrstates+1], justify="centre")
				print(prettyallcladestates, right=FALSE)
				cat(sep="","\n","Summary: most parsimonious state(s) at key nodes\n")
				prettybestcladestates <- bestcladestates
				prettybestcladestates[,1] <- format(prettybestcladestates[,1], justify="left")
				prettybestcladestates[,2] <- format(prettybestcladestates[,2], justify="centre")
				print(prettybestcladestates, right=FALSE)
			}
			if (method=="ML") {
				bestcladestates <- matrix(,nrow=nrmapnodes,ncol=2) #empty matrix to store most likely clade states
				for (i in 1:nrmapnodes)
					{
						favstate <- which.max(allcladestates[i,-1])[[1]] #find most probable state index
						bestcladestates[i,1] <- char$states[favstate]
						bestcladestates[i,2] <- allcladestates[i,favstate+1]
					}
				bestcladestates <- cbind(data.frame(mapclades), data.frame(bestcladestates))
				colnames(bestcladestates) <- c("Node","ML state","Prob")
				cat(sep="","\n","Ancestral states for key nodes\n")
				allcladestates[,2:char$nrstates+1] <- format(allcladestates[,2:char$nrstates+1], scientific=FALSE)
				prettyallcladestates <- allcladestates
				prettyallcladestates[,1] <- format(prettyallcladestates[,1], justify="left")
				prettyallcladestates[,2:char$nrstates+1] <- format(prettyallcladestates[,2:char$nrstates+1], justify="centre")
				print(prettyallcladestates, right=FALSE)
				cat(sep="","\n","Summary: most likely state at key nodes\n")
				prettybestcladestates <- bestcladestates
				prettybestcladestates[,1] <- format(prettybestcladestates[,1], justify="left")
				prettybestcladestates[,2:3] <- format(prettybestcladestates[,2:3], justify="centre")
				print(prettybestcladestates, right=FALSE)
				cat(sep="","\n\n")
			}
		}

		if (method=="rjMCMC") {
			mapstates <- matrix(,nrow=nrmapnodes,ncol=char$nrstates)
			mapmean <- matrix(,nrow=nrmapnodes,ncol=char$nrstates)
			maphpdlb <- matrix(,nrow=nrmapnodes,ncol=char$nrstates)
			maphpdub <- matrix(,nrow=nrmapnodes,ncol=char$nrstates)
			bestcladestates <- matrix(,nrow=nrmapnodes,ncol=5) #empty matrix to store most probable clade states
			for (i in 1:nrmapnodes) {
				for (j in 1:char$nrstates) {
					btcolname <- paste(mapclades[i], " P(", j-1, ")", sep="")
					colnr <- which(colnames(ans$cleanchain)==btcolname)
					if (length(colnr)==0) {
						if (mapclades[i]=="Superasteridae") {btcolname <- paste("superasterids", " P(", j-1, ")", sep="")}
						if (mapclades[i]=="Superrosidae") {btcolname <- paste("superrosids", " P(", j-1, ")", sep="")}
						if (mapclades[i]=="Monocotyledoneae") {btcolname <- paste("Monocotyledonae", " P(", j-1, ")", sep="")}
						if (mapclades[i]=="Eudicotyledoneae") {btcolname <- paste("Eudicotyledonae", " P(", j-1, ")", sep="")}
						colnr <- which(colnames(ans$cleanchain)==btcolname)
					}
					mapmean[i,j] <- format(round(mean(ans$cleanchain[,colnr]), digits=4), scientific=FALSE)
					hpd <- format(round(boa.hpd(ans$cleanchain[,colnr], 0.05), digits=4), scientific=FALSE)
					maphpdlb[i,j] <- hpd[1]
					maphpdub[i,j] <- hpd[2]
					mapstates[i,j] <- paste(mapmean[i,j], " (", maphpdlb[i,j], "-", maphpdub[i,j], ")", sep="")
				}
				favstate <- which.max(mapmean[i,])[[1]] #find most probable state index
				bestcladestates[i,1] <- char$states[favstate]
				bestcladestates[i,2] <- mapstates[i,favstate]
				bestcladestates[i,3] <- mapmean[i,favstate]
				bestcladestates[i,4] <- maphpdlb[i,favstate]
				bestcladestates[i,5] <- maphpdub[i,favstate]
			}
			allcladestates <- cbind(data.frame(mapclades), mapstates)
			colnames(allcladestates) <- c("Node",char$states)
			cat(sep="","\n","Ancestral states for key nodes (mean and 95% HPD)\n")
			prettyallcladestates <- allcladestates
			prettyallcladestates[,1] <- format(prettyallcladestates[,1], justify="left")
			prettyallcladestates[,2:(char$nrstates+1)] <- format(prettyallcladestates[,2:(char$nrstates+1)], justify="centre")
			print(prettyallcladestates, right=FALSE)
			bestcladestates <- cbind(data.frame(mapclades), data.frame(bestcladestates))
			colnames(bestcladestates) <- c("Node","rjMCMC state","Prob (mean and 95% HPD)","P (mean)","P (95% lb)","P (95% ub)")
			cat(sep="","\n","Summary: most probable state at key nodes\n")
			prettybestcladestates <- bestcladestates[,1:3]
			prettybestcladestates[,1] <- format(prettybestcladestates[,1], justify="left")
			prettybestcladestates[,2:3] <- format(prettybestcladestates[,2:3], justify="centre")
			print(prettybestcladestates, right=FALSE)
			cat(sep="","\n\n")
		}
		
		flush.console() #print immediately to console to follow progress (otherwise delayed until function or loop finished)

	#Prepare edge colors

		tips <- tree$tip.label #Prepare a vector with the complete list of tips in "tree" ("tips")
		#Calculate number of branches (required for mapping colors in asr function): this would normally be equal to 2*nrtips-2 unless there are polytomies
			nredges <- nrtips + treeUP$Nnode - 1 #total number of branches (edges) in tree: number of terminal branches (tips) + number of internal branches (= number of internal nodes - 1)
		edgeco <- rep("grey",nredges)
		edgecoUP <- rep("grey",nredges) #separate set of edge colors for upright (inverse-ladderized) tree

		if ((method=="MP") | (method=="ML")) {
			
			#Prepare edge colors for 'tree'
				for (i in 1:nredges)
					{
						edgeendnode <- tree$edge[[i,2]]
						if (edgeendnode<=nrtips) { #edge is terminal -> color with tip state
							tipmatch <-match(node.leaves(tree, edgeendnode),names(char$charData))
							if (is.na(tipmatch)==FALSE) { #tip is matched in character matrix
								tipstate <- char$charData[[tipmatch]]
								if (is.na(tipstate)==FALSE) { #tip state is not missing data (NA)
									edgeco[i] <- char$subco[tipstate+1]
								}
							}
						} else { #edge is internal -> color with most probable state of the node it leads to
							asrnode <- getMRCA(treeUP,node.leaves(tree,edgeendnode)) #match node of 'tree' (being mapped) with node of 'treeUP' (which was used for ASR); although the two nodes are phylogenetically identical, they have different numbers (because ladderization affects node number in ape)
							if (method=="ML") {
								favstate <- which.max(ancstates[asrnode-nrtips,])[[1]] #find most probable state index
								if (ancstates[[asrnode-nrtips,favstate]]>bct) { #color only if probability of most probable state is > bct (branch coloring threshold set in parameters)
									edgeco[i] <- char$subco[favstate]
								}
							}
							if (method=="MP") {
								favstate <- which(ancstates[asrnode-nrtips,]==1) #find unique most parsimonious state index (return no value in case of multiple equally most parsimonious states)
								if (length(favstate)>0) {
									edgeco[i] <- char$subco[favstate]
								}
							}
						}
					}

			#Prepare edge colors for 'treeUP'
				for (i in 1:nredges)
					{
						edgeendnode <- treeUP$edge[[i,2]]
						if (edgeendnode<=nrtips) { #edge is terminal -> color with tip state
							tipmatch <-match(node.leaves(treeUP, edgeendnode),names(char$charData))
							if (is.na(tipmatch)==FALSE) { #tip is matched in character matrix
								tipstate <- char$charData[[tipmatch]]
								if (is.na(tipstate)==FALSE) { #tip state is not missing data (NA)
									edgecoUP[i] <- char$subco[tipstate+1]
								}
							}
						} else { #edge is internal -> color with most probable state of the node it leads to
							if (method=="ML") {
								favstate <- which.max(ancstates[edgeendnode-nrtips,])[[1]] #find most probable state index
								if (ancstates[[edgeendnode-nrtips,favstate]]>bct) { #color only if probability of most probable state is > bct (branch coloring threshold set in parameters)
									edgecoUP[i] <- char$subco[favstate]
								}
							}
							if (method=="MP") {
								favstate <- which(ancstates[edgeendnode-nrtips,]==1) #find unique most parsimonious state index (return no value in case of multiple equally most parsimonious states)
								if (length(favstate)>0) {
									edgecoUP[i] <- char$subco[favstate]
								}
							}
						}
					}
		}

	list(mkmodel=mkmodel,ans=ans,ancstates=ancstates,mapstates=mapstates,bestcladestates=bestcladestates,edgeco=edgeco,edgecoUP=edgecoUP,rjsumtable=rjsumtable)
}