#Function by Hervé Sauquet (2016)
#Sets the list of discrete (Markov) models to test, according to the number of character states (nrstates), type of character (determined from charname), model option (modeloption), and root mode (rootmode)

listMkmodels <- function(nrstates, charname, modeloption, rootmode)
{
	mkmodels <- c("ARD","ARDeq","ER") #default (if unchanged, corresponding to modeloption="simple" and rootmode="both")
	if (modeloption=="comprehensive") {
		if (nrstates==2) { mkmodels <- c("ARD","ARDeq","ER","UNI01","UNI10") }
		if (nrstates>=3) {
			mkmodels <- c("ARD","ARDeq","ER","SYM","SYMeq")
			if ((grepl("number",charname,ignore.case=TRUE)==TRUE) | (grepl("merism",charname,ignore.case=TRUE)==TRUE)) {
				mkmodels <- c(mkmodels,"ORD","ORDeq","ORDSYM","ORDSYMeq","ORDER")
			}
		}
	}
	if (modeloption=="simple") {
		if (nrstates>=3) {
			mkmodels <- c("ER","SYM","SYMeq")
#		  mkmodels <- c("ER")
		}
	}
	if (rootmode=="flat") { mkmodels <- mkmodels[grepl("eq",mkmodels)==FALSE] } #trim (subset) models containing "eq"
	if (rootmode=="equilibrium") { mkmodels <- mkmodels[(mkmodels %in% c("ARD","SYM","ORD","ORDSYM")==FALSE)] } #trim (subset) models containing "eq"
	return(mkmodels)
}