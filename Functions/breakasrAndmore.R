#Function by Hervé Sauquet (2017)

# Temp set-up
	# subtree=subtreeB
	# tree=usrdata$tree
	# treeUP=usrdata$treeUP
	# charmatrix=usrdata$charmatrix
	# charindex=charindex
	# mkmodel="ER"
	# mapclades=mapcladesB
	# mapnodesUP=mapnodesB
	# char=char
	# co=usrcolors$co
	# bct=bct
	# method="ML"
	# BayesTraitsfolder="./"

breakasrAndmore <- function(subtree, tree, treeUP, charmatrix, charindex, mkmodel, mapclades, mapnodesUP, char, co, bct, method="ML", BayesTraitsfolder="./")
{

	nrtips <- length(tree$tip.label)
	nrsubtips <- length(subtree$tip.label)
	rjsumtable <- ""

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
				ans <- rayDISC(subtree,charmatrix,ntraits=1,charnum=charindex,model=basemkmodel,node.states="marginal",root.p=rootoption)
			} else {
				ans <- rayDISC(subtree,charmatrix,ntraits=1,charnum=charindex,rate.mat=rate.mat,node.states="marginal",root.p=rootoption)
			}
	
		print(ans)
		ancstates <- ans$states
		flush.console() #print immediately to console to follow progress (otherwise delayed until function or loop finished)
		
	}

	#Map and summarize focal node states

		mapstates <- ancstates #needed in fabTree only (to plot pie charts for key clades only in fan mode)
		nrmapnodes <- length(mapclades)
		
		if (nrmapnodes > 0) { #the partial rewrite below was required for specific cases of subtrees with 0 or 1 mapclade only
			if ((method=="MP") | (method=="ML")) {
				allcladestates <- matrix(nrow=nrmapnodes, ncol=length(char$states)+1)
				allcladestates <- data.frame(allcladestates)
				colnames(allcladestates) <- c("Node",char$states)
				for (i in 1:nrmapnodes)
					{
						allcladestates[i,1] <- mapclades[i]
						allcladestates[i,2:(char$nrstates+1)] <- round(ancstates[mapnodesUP[i]-nrsubtips,], digits=4)
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
		} else {
			bestcladestates=""
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
							asrnode <- getMRCA(subtree,tips(tree,edgeendnode)) #match node of 'tree' (being mapped) with node of 'treeUP' (which was used for ASR); although the two nodes are phylogenetically identical, they have different numbers (because ladderization affects node number in ape)
							if (is.null(asrnode)==TRUE) {
								edgeco[i] <- "red"
							} else {
								favstate <- which.max(ancstates[asrnode-nrsubtips,])[[1]] #find most probable state index
								if (ancstates[[asrnode-nrsubtips,favstate]]>bct) { #color only if probability of most probable state is > bct (branch coloring threshold set in parameters)
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
							asrnode <- getMRCA(subtree,node.leaves(treeUP,edgeendnode)) #match node of 'tree' (being mapped) with node of 'treeUP' (which was used for ASR); although the two nodes are phylogenetically identical, they have different numbers (because ladderization affects node number in ape)
							if (is.null(asrnode)== FALSE) {
								favstate <- which.max(ancstates[asrnode-nrsubtips,])[[1]] #find most probable state index
								if (ancstates[[asrnode-nrsubtips,favstate]]>bct) { #color only if probability of most probable state is > bct (branch coloring threshold set in parameters)
									edgecoUP[i] <- char$subco[favstate]
								}
							}
						}
					}
		}

	list(mkmodel=mkmodel,ans=ans,ancstates=ancstates,mapstates=mapstates,bestcladestates=bestcladestates,edgeco=edgeco,edgecoUP=edgecoUP,rjsumtable=rjsumtable)
}