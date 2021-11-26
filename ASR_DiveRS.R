#Ancestral state reconstruction for DiveRS (Hervé Sauquet, November 2021)
#(adapted from ASR_eFLOWER-1209.R)

#This script reconstructs ancestral character states on a given phylogeny using maximum likelihood
	# Inputs: see readData()
	# Outputs: PDF(s) of character optimization(s), log files with summary results of ASR, table with summary data for model comparison

#Clear buffer of all objects
	rm(list=ls())

#Set working directory
	#UNIX: no path, comment out the three lines above; to run the script, type "R CMD BATCH ASR_eFLOWER-vX.R"
	#PC from Notepad++ with correct plugin: no path required; to run the script, select all and press F8
		#setwd("C:/Users/Herv?/Dropbox (MAGNIPHY)/Documents/WorkDocs/eFLOWER/eFLOWER Analyses/R") #PC path (only if executing code directly from the R GUI)
		#setwd("/Users/hervesauquet/Dropbox (MAGNIPHY)/Documents/WorkDocs/eFLOWER/eFLOWER Analyses/R") #Mac path
		#Sys.setlocale('LC_ALL','C') #apparently now required on my Mac

#Load packages
	library(ape) #framework for handling phylogenetic trees
	library(corHMM) #used for ML ancestral state reconstruction (rayDISC function)
	library(phangorn) #used for MP ancestral state reconstruction (ancestral.pars function)
	library(geiger) #required for function tips
	library(plotrix) #required for plotting tables
	library(xlsx) #required for writing summary results in multi-sheet Excel workbooks
	library(numbers) #required for calculating Bell numbers
	library(boa) #required for calculating 95% HPDs

#Load custom functions
	sapply(list.files(pattern="[.]R$", path="./Functions/", full.names=TRUE), source);

#Optional debug mode
#	options(error = recover)

#Read data
	dsversion <- ""
	series <- ""
	suffix <- "traits"
	treefilename <- "tree_298.tre"
	usrdata <- readData(
		datafolder="./Data/", 
		treefile=treefilename, 
		maptreefile="", 
		datafile="cleandata.csv", 
		charnamefile="datanames.csv", 
		ingroupdef="", 
		cladedefsfile="CladeDefinitions_DiveRS.csv",
		trimtotree=TRUE,
		stretchtotree=TRUE
	)

#Reconstruct ancestral states for all (or a subset of) characters, fitting multiple models for each:
	runIDprefix <- paste("ASR_DiveRS_", dsversion, series, suffix, sep="")
	loopASR(
		usrdata, 
#		methodoption="MPMLrjMCMC", 
		methodoption="MPML", 
		modeloption="comprehensive", 
#    modeloption="simple", 
		rootmode="both", 
#    rootmode="equilibrium", 
		bct=0.5, 
		colorscheme="pastel4", 
		mapallmodels=FALSE, 
		addfantree=TRUE, 
		writedata=FALSE, 
		charlist="", 
		testmode=FALSE, 
		runIDprefix=runIDprefix, 
		resformat="Excel", 
		BayesTraitsfolder=paste("./BayesTraits_", dsversion, "_", series, suffix, "/", sep=""), 
		graphpar=c(45,8,20,0.12,0.05,4.5,3.5,-40)
	)
