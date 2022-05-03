#Function by Hervé Sauquet (2016, 2019)

#Store node state colors in a vector (colors must be provided in the same exact order as above)
#Important: make sure there are at least as many elements in the lists below (currently 12) as in the character with the highest number of character states in the matrix

setColors <- function(scheme="bright")
{
	co <- c("yellow","blue","green","red","purple","orange","brown","black","black","black","black","black") #plain colors
	#pastel colors optimized for figure production
		if (scheme=="pastel1") co <- c("springgreen3","gold1","dodgerblue3","firebrick1","darkorchid","tan4") #Fumarioideae scheme used in Sauquet et al (2015)
		if (scheme=="pastel2") co <- c("gold1","dodgerblue3","springgreen3","firebrick1","darkorchid","tan4") #Proteaceae scheme used in Citerne et al (2017)
		#if (scheme=="pastel3") co <- c("dodgerblue3","springgreen3","gold1","firebrick1","darkorchid","tan4") #superbiome scheme used in Ramírez-Barahona et al In prep (AngioTimeTree2)
		if (scheme=="pastel3") co <- c("dodgerblue3","gold1","firebrick1","springgreen3","darkorchid","tan4") #superbiome scheme used in Ramírez-Barahona et al In prep (AngioTimeTree2)
		if (scheme=="pastel4") co <- c("gold1","dodgerblue3","springgreen3","firebrick1","darkorchid","tan4","yellow","blue","green","red","purple","orange","brown","black","black","black","black","black")
		if (scheme=="pastel5") co <- c("firebrick1","gold1","springgreen3","dodgerblue3","darkorchid","tan4","yellow","blue","green","red","purple","orange","brown","black","black","black","black","black")
		circo <- c("black","black","black","black","black","black","black","black","black","black","black","black") #color for drawing circle border of tip states
	list(co=co,circo=circo)
}