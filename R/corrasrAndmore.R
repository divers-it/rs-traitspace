#Function by Hervé Sauquet (2017)

#Primer (set-up mode)
	# tree = usrdata$tree
	# treeUP = usrdata$treeUP
# #	charmatrix = usrdata$charmatrix
	# charindexlist=c(2,21)
	# mkmodel = "ARDeqnodual"
	# mapclades = usrdata$mapclades
	# mapnodesUP = usrdata$mapnodesUP
	# co = usrcolors$co
	# bct = bct
	# # method="rjMCMC"

corrasrAndmore <- function(tree, treeUP, charpairmatrix, charindexlist, mkmodel, mapclades, mapnodesUP, corrchar, co, bct, method="ML", BayesTraitsfolder="./")
{

	nrtips <- length(tree$tip.label)
	rjsumtable <- ""

	if (method=="ML") {

		cat(sep="","Performing maximum likelihood ASR (using rayDISC) with model ",mkmodel," (character ", corrchar$charnr, ")...\n\n")
		flush.console() #print immediately to console to follow progress (otherwise delayed until function or loop finished)
		basemkmodel <- gsub("eq","",mkmodel)
	
		#Prepare rate matrix (not needed for ER and ARD models)
			if (! ((basemkmodel=="ER") | (basemkmodel=="ARD")) ) {
				rate.mat <- rate.mat.maker(rate.cat=1,hrm=FALSE,ntraits=1,4,model="ARD") #prepares most general matrix for a pair of binary traits
				#IMPORTANT: models and rate matrices are defined below for a pair of 2-state characters
					if (grepl("nodual",basemkmodel)==TRUE) { #drop diagonal of dual transitions (replicates corDISC)
						rate.mat <- rate.par.drop(rate.mat,drop.par=c(10))
						rate.mat <- rate.par.drop(rate.mat,drop.par=c(8))
						rate.mat <- rate.par.drop(rate.mat,drop.par=c(5))
						rate.mat <- rate.par.drop(rate.mat,drop.par=c(3))
					}
					if (grepl("SYM",basemkmodel)==TRUE) { #this actually already exists natively in rayDISC (but looks like it does something different)
						for (i in 1:(corrchar$nrstates*(corrchar$nrstates-1)/2)) {
							rowi <- which(rate.mat==i,arr.ind=TRUE)[,1]
							coli <- which(rate.mat==i,arr.ind=TRUE)[,2]
							rate.mat <- rate.par.eq(rate.mat,eq.par=c(i,rate.mat[coli,rowi]))
						}
					}
					if (grepl("UNCORR",basemkmodel)==TRUE) { #use the symmetry code above, then rotate (this works only for the nodual version, which has identical NA diagonals; not sure the code would be general enough for larger matrices, need to check)
						sym.mat <- rate.mat
						for (i in 1:(corrchar$nrstates*(corrchar$nrstates-1)/2)) {
							rowi <- which(sym.mat==i,arr.ind=TRUE)[,1]
							coli <- which(sym.mat==i,arr.ind=TRUE)[,2]
							sym.mat <- rate.par.eq(sym.mat,eq.par=c(i,sym.mat[coli,rowi]))
						}
						rate.mat[,] <- t(apply(sym.mat,2,rev))[,]
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
				corrans <- rayDISC(treeUP,charpairmatrix,ntraits=1,charnum=4,model=basemkmodel,node.states="marginal",root.p=rootoption)
				# corrans <- corDISC(treeUP,charpairmatrix,ntraits=2,model=basemkmodel,node.states="marginal",root.p=rootoption)
			} else {
				corrans <- rayDISC(treeUP,charpairmatrix,ntraits=1,charnum=4,rate.mat=rate.mat,node.states="marginal",root.p=rootoption)
				# corrans <- corDISC(treeUP,charpairmatrix,ntraits=2,rate.mat=rate.mat,node.states="marginal",root.p=rootoption)
			}
	
		print(corrans)
		ancstates <- corrans$states
		flush.console() #print immediately to console to follow progress (otherwise delayed until function or loop finished)
		
	}

	if (method=="rjMCMC") {
	
		cat(sep="","Summarizing results from reversible-jump MCMC BayesTraits log file (character ", corrchar$charnr, ")...\n\n")
		flush.console() #print immediately to console to follow progress (otherwise delayed until function or loop finished)
		
		#Read results from a previously run BayesTraits analysis (and summarize models visited)
			corrans <- readBayesTraits(datafolder=BayesTraitsfolder, char=corrchar, , corr=TRUE)
			ancstates <- ""
		
		#Summarize basic statistics on number of parameters and transition rates
		
			#First, summarize number of free and zero parameters (for each model visited by the chain)
				meannrpar <- format(round(mean(corrans$cleanchain[,5]), digits=1), scientific=FALSE) #(note that the first column in corrans$cleanchain is teh Iteration column)
				hpdnrpar <- format(round(boa.hpd(corrans$cleanchain[,5], 0.05), digits=1), scientific=FALSE)
				meannrzero <- format(round(mean(corrans$cleanchain[,6]), digits=1), scientific=FALSE)
				hpdnrzero <- format(round(boa.hpd(corrans$cleanchain[,6], 0.05), digits=1), scientific=FALSE)
				rjstatstable <- data.frame(rbind(cbind(meannrpar, meannrzero), cbind(hpdnrpar, hpdnrzero)))
			
			#Then, summarize transition rates
				# for (i in 1:length(char$ratecolnames)) { #(note that length(char$ratecolnames) should be equal to corrchar$nrstates * (corrchar$nrstates - 1))
					# meanq <- format(round(mean(corrans$cleanchain[,6+i]), digits=4), scientific=FALSE)
					# hpdq <- format(round(boa.hpd(corrans$cleanchain[,6+i], 0.05), digits=4), scientific=FALSE)
					# rjstatstable <- cbind(rjstatstable, matrix(c(meanq, hpdq), ncol=1))
					# mcmcPlot(char$charID, , iter=corrans$cleanchain[,1], corrans$cleanchain[,6+i], char$ratecolnames[i])
					# mcmcPlot(char$charID, , , corrans$cleanchain[,6+i], char$ratecolnames[i])
				# }
			
			#Tidy up and print
				# colnames(rjstatstable) <- c("Npar","Nzero", char$ratecolnames)
				colnames(rjstatstable) <- c("Npar","Nzero")
				rownames(rjstatstable) <- c("mean","95% lb","95% ub")
				# cat(sep="","\n","Summary statistics for number of parameters and transition rates\n")
				cat(sep="","\n","Summary statistics for number of parameters\n")
				prettyrjstatstable <- format(rjstatstable, justify="centre")
				print(prettyrjstatstable, quote=FALSE, right=FALSE)
				rjsumtable <- list(rjmodeltable=corrans$rjmodeltable, rjstatstable=rjstatstable)

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

			colnames(allcladestates) <- c("Node",corrchar$states)
			if (method=="MP") {
				bestcladestates <- matrix(,nrow=nrmapnodes,ncol=1) #empty matrix to store most parsimonious clade states
				for (i in 1:nrmapnodes)
					{
						favstate <- which(allcladestates[i,-1]>0) #find most parsimonious state indices
						bestcladestates[i,1] <- paste(corrchar$states[favstate], sep=" / ", collapse=" / ")
					}
				bestcladestates <- cbind(data.frame(mapclades), data.frame(bestcladestates))
				colnames(bestcladestates) <- c("Node","MP state(s)")
				cat(sep="","Ancestral states for key nodes\n")
				prettyallcladestates <- allcladestates
				prettyallcladestates[,1] <- format(prettyallcladestates[,1], justify="left")
				prettyallcladestates[,2:corrchar$nrstates+1] <- format(prettyallcladestates[,2:corrchar$nrstates+1], justify="centre")
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
						bestcladestates[i,1] <- corrchar$states[favstate]
						bestcladestates[i,2] <- allcladestates[i,favstate+1]
					}
				bestcladestates <- cbind(data.frame(mapclades), data.frame(bestcladestates))
				colnames(bestcladestates) <- c("Node","ML state","Prob")
				cat(sep="","\n","Ancestral states for key nodes\n")
				allcladestates[,2:corrchar$nrstates+1] <- format(allcladestates[,2:corrchar$nrstates+1], scientific=FALSE)
				prettyallcladestates <- allcladestates
				prettyallcladestates[,1] <- format(prettyallcladestates[,1], justify="left")
				prettyallcladestates[,2:corrchar$nrstates+1] <- format(prettyallcladestates[,2:corrchar$nrstates+1], justify="centre")
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
			mapstates <- matrix(,nrow=nrmapnodes,ncol=corrchar$nrstates)
			mapmean <- matrix(,nrow=nrmapnodes,ncol=corrchar$nrstates)
			maphpdlb <- matrix(,nrow=nrmapnodes,ncol=corrchar$nrstates)
			maphpdub <- matrix(,nrow=nrmapnodes,ncol=corrchar$nrstates)
			bestcladestates <- matrix(,nrow=nrmapnodes,ncol=5) #empty matrix to store most probable clade states
			for (i in 1:nrmapnodes) {
				for (j in 1:corrchar$nrstates) {
					btcolname <- paste(mapclades[i], " - P(", corrchar$corrsymbols[j], ")", sep="")
					colnr <- which(colnames(corrans$cleanchain)==btcolname)
					if (length(colnr)==0) {
						if (mapclades[i]=="Superasteridae") {btcolname <- paste("superasterids", " P(", j-1, ")", sep="")}
						if (mapclades[i]=="Superrosidae") {btcolname <- paste("superrosids", " P(", j-1, ")", sep="")}
						if (mapclades[i]=="Monocotyledoneae") {btcolname <- paste("Monocotyledonae", " P(", j-1, ")", sep="")}
						if (mapclades[i]=="Eudicotyledoneae") {btcolname <- paste("Eudicotyledonae", " P(", j-1, ")", sep="")}
						colnr <- which(colnames(corrans$cleanchain)==btcolname)
					}
					mapmean[i,j] <- format(round(mean(corrans$cleanchain[,colnr]), digits=4), scientific=FALSE)
					hpd <- format(round(boa.hpd(corrans$cleanchain[,colnr], 0.05), digits=4), scientific=FALSE)
					maphpdlb[i,j] <- hpd[1]
					maphpdub[i,j] <- hpd[2]
					mapstates[i,j] <- paste(mapmean[i,j], " (", maphpdlb[i,j], "-", maphpdub[i,j], ")", sep="")
				}
				favstate <- which.max(mapmean[i,])[[1]] #find most probable state index
				bestcladestates[i,1] <- corrchar$states[favstate]
				bestcladestates[i,2] <- mapstates[i,favstate]
				bestcladestates[i,3] <- mapmean[i,favstate]
				bestcladestates[i,4] <- maphpdlb[i,favstate]
				bestcladestates[i,5] <- maphpdub[i,favstate]
			}
			allcladestates <- cbind(data.frame(mapclades), mapstates)
			colnames(allcladestates) <- c("Node",corrchar$states)
			cat(sep="","\n","Ancestral states for key nodes (mean and 95% HPD)\n")
			prettyallcladestates <- allcladestates
			prettyallcladestates[,1] <- format(prettyallcladestates[,1], justify="left")
			prettyallcladestates[,2:(corrchar$nrstates+1)] <- format(prettyallcladestates[,2:(corrchar$nrstates+1)], justify="centre")
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

		# tips <- tree$tip.label #Prepare a vector with the complete list of tips in "tree" ("tips")
		# #Calculate number of branches (required for mapping colors in asr function): this would normally be equal to 2*nrtips-2 unless there are polytomies
			# nredges <- nrtips + treeUP$Nnode - 1 #total number of branches (edges) in tree: number of terminal branches (tips) + number of internal branches (= number of internal nodes - 1)
		# edgeco <- rep("grey",nredges)
		# edgecoUP <- rep("grey",nredges) #separate set of edge colors for upright (inverse-ladderized) tree

		# if ((method=="MP") | (method=="ML")) {
			
			# #Prepare edge colors for 'tree'
				# for (i in 1:nredges)
					# {
						# edgeendnode <- tree$edge[[i,2]]
						# if (edgeendnode<=nrtips) { #edge is terminal -> color with tip state
							# tipmatch <-match(node.leaves(tree, edgeendnode),names(char$charData))
							# if (is.na(tipmatch)==FALSE) { #tip is matched in character matrix
								# tipstate <- char$charData[[tipmatch]]
								# if (is.na(tipstate)==FALSE) { #tip state is not missing data (NA)
									# edgeco[i] <- char$subco[tipstate+1]
								# }
							# }
						# } else { #edge is internal -> color with most probable state of the node it leads to
							# asrnode <- getMRCA(treeUP,node.leaves(tree,edgeendnode)) #match node of 'tree' (being mapped) with node of 'treeUP' (which was used for ASR); although the two nodes are phylogenetically identical, they have different numbers (because ladderization affects node number in ape)
							# if (method=="ML") {
								# favstate <- which.max(ancstates[asrnode-nrtips,])[[1]] #find most probable state index
								# if (ancstates[[asrnode-nrtips,favstate]]>bct) { #color only if probability of most probable state is > bct (branch coloring threshold set in parameters)
									# edgeco[i] <- char$subco[favstate]
								# }
							# }
							# if (method=="MP") {
								# favstate <- which(ancstates[asrnode-nrtips,]==1) #find unique most parsimonious state index (return no value in case of multiple equally most parsimonious states)
								# if (length(favstate)>0) {
									# edgeco[i] <- char$subco[favstate]
								# }
							# }
						# }
					# }

			# #Prepare edge colors for 'treeUP'
				# for (i in 1:nredges)
					# {
						# edgeendnode <- treeUP$edge[[i,2]]
						# if (edgeendnode<=nrtips) { #edge is terminal -> color with tip state
							# tipmatch <-match(node.leaves(treeUP, edgeendnode),names(char$charData))
							# if (is.na(tipmatch)==FALSE) { #tip is matched in character matrix
								# tipstate <- char$charData[[tipmatch]]
								# if (is.na(tipstate)==FALSE) { #tip state is not missing data (NA)
									# edgecoUP[i] <- char$subco[tipstate+1]
								# }
							# }
						# } else { #edge is internal -> color with most probable state of the node it leads to
							# if (method=="ML") {
								# favstate <- which.max(ancstates[edgeendnode-nrtips,])[[1]] #find most probable state index
								# if (ancstates[[edgeendnode-nrtips,favstate]]>bct) { #color only if probability of most probable state is > bct (branch coloring threshold set in parameters)
									# edgecoUP[i] <- char$subco[favstate]
								# }
							# }
							# if (method=="MP") {
								# favstate <- which(ancstates[edgeendnode-nrtips,]==1) #find unique most parsimonious state index (return no value in case of multiple equally most parsimonious states)
								# if (length(favstate)>0) {
									# edgecoUP[i] <- char$subco[favstate]
								# }
							# }
						# }
					# }
		# }

	list(mkmodel=mkmodel,corrans=corrans,ancstates=ancstates,mapstates=mapstates,bestcladestates=bestcladestates,rjsumtable=rjsumtable)
}