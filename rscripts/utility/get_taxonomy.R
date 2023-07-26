rm(list=ls())

#load libraries
library(taxize)

#load data set
df<-readRDS(file = here::here("outputs/df_filt_trans.rds"))

#get list of species in dataset
specieslist <- rownames(df)

#loop through species and get family and order
#try in case there are errors
for(i in 1:length(specieslist)){
  if(i == 1){
    try(res<-tax_name(query = specieslist[i], get = c("family","order"), db = "ncbi"))
  } else {
    try(res<-rbind(res,tax_name(query = specieslist[i], get = c("family","order"), db = "ncbi")))
  }
}
head(res)

#make new results dataframe
res2<-res

#get species that were not retrieved first time
sp_with_err <- setdiff(rownames(df),res$query)

#get info for species that didn't work first time
#NOTE: this could serve as first loop but probably won't need to run again
for(i in 1:length(sp_with_err)){
  
  #keeps j null each iteration
  j<-NULL 
  
  #keeps trying to get info while j is null
  while(is.null(j)) {
    try(
      j <- tax_name(query = sp_with_err[i], get = c("family","order"), db = "ncbi")
    )
  }
  
  res2<-rbind(res2,j)
  
}

#check all species are present
setdiff(rownames(df),res2$query)
setdiff(res2$query,rownames(df))

#reorder taxonomy results
res2<-res2[order(res2$query),]

#check order
rownames(df)==res2$query

#reformat df
res2_out<-res2[,c(3:4)]
rownames(res2_out)<-res2$query

#NOTE: Temporary fix for error because of synonym
#Pontederia crassipes (formerly Eichhornia crassipes)

#family
res2_out["Eichhornia crassipes",][1]<-"Pontederiaceae"
#order
res2_out["Eichhornia crassipes",][2]<-"Commelinales"

#NOTE: Temporary fix for unknown error in Quintinia hyehenensis

#family
res2_out["Quintinia hyehenensis",][1]<-"Paracryphiaceae"
#order
res2_out["Quintinia hyehenensis",][2]<-"Paracryphiales"

#check to see if results make sense
res2_out[res2_out$order=="Poales",]

#output taxonomy
saveRDS(res2_out, file = here::here("outputs/taxonomy.rds"))
