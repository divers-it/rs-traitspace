#A chi2 test is not valid if: (1) at least one of the expected values in the frequency table is equal to zero, (2) more than a fifth of the expected values in the frequency table are inferior to 5.
makechi2=function(categ){ #function performing a chi2 test between each pair of columns given in "categ"
  for(i in 1:ncol(categ)){categ[,i]=as.factor(categ[,i])}
  ncomp=sum(1:(ncol(categ)-1)) #nombre de comparaisons
  RESchi2=vector("list",3*ncomp)
  N=1
  oldnames=NULL
  nc=ncol(categ)
  
  for(a in 1:(nc-1)){
    for (b in (a+1):nc){
      subcateg=as.data.frame(categ[,c(a,b)]) #sampling of two categories, that will be tested for correlations
      
      i=1
      while(i < nrow(subcateg)){
        while((is.na(subcateg[i,1])==T|is.na(subcateg[i,2])==T)&(i<=nrow(subcateg))){    #to remove lines containing NAs
          subcateg=subcateg[-i,]
        }
        
        if((nchar(as.character(subcateg[i,1]))>1)&(i<=nrow(subcateg))){ #is there polymorphism in the first column
          states=unlist(strsplit(as.character(subcateg[i,1]),"&"))
          subcateg[i,1]=states[1]
          for(k in 2:length(states)){
            newrow=as.data.frame(t(c(states[k],as.character(subcateg[i,2]))))
            colnames(newrow)=colnames(subcateg)
            subcateg=rbind(subcateg,newrow)
          }
        }
        if((nchar(as.character(subcateg[i,2]))>1)&(i<=nrow(subcateg))){ #is there polymorphism in the second column
          states=unlist(strsplit(as.character(subcateg[i,2]),"&"))
          subcateg[i,2]=states[1]
          for(k in 2:length(states)){
            newrow=as.data.frame(t(c(as.character(subcateg[i,1]),states[k])))
            colnames(newrow)=colnames(subcateg)
            subcateg=rbind(subcateg,newrow)
          }
        }
        i=i+1
      }
      #construction of a frequency table to perform the chi2 test
      lignes=levels(droplevels(subcateg[,1])) #levels have to be updated
      colones=levels(droplevels(subcateg[,2]))
      FreqTabl=matrix(0,ncol=length(colones), nrow=length(lignes))
      rownames(FreqTabl)=lignes
      colnames(FreqTabl)=colones
      for(i in 1:length(lignes)){
        for(j in 1:length(colones)){
          FreqTabl[i,j]=length(which((subcateg[,1]==lignes[i])&(subcateg[,2]==colones[j])))
        }  
      }
      
      Chitest=chisq.test(FreqTabl)
      
      RESchi2[[N]]=FreqTabl
      N=N+1
      RESchi2[[N]]=Chitest$expected
      N=N+1
      RESchi2[[N]]=Chitest
      N=N+1
      RESchi2[[N]]=Chitest$residuals
      N=N+1
      
      nomTests=paste(colnames(subcateg[1]),colnames(subcateg[2]),sep=" VS ")
      newnames=c(paste(nomTests," FREQ", sep=""), paste(nomTests," EXP", sep=""), paste(nomTests," CHI2", sep=""), paste(nomTests," PEARS", sep=""))
      oldnames=append(oldnames,newnames)
    }  
  }
  
  names(RESchi2)=oldnames
  return(RESchi2)
}
