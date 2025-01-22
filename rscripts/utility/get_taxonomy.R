rm(list=ls())

# load libraries
library(taxize)

# load data set
df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

# get list of species in dataset
specieslist <- rownames(df)

# 

####
## Get number of families with TNRS ----
####

# check names with TNRS
tnrs_check<-TNRS::TNRS(specieslist, source = "wcvp")

# number of families         
length(table(tnrs_check$Name_matched_accepted_family))

####
## Get number of orders (and families) with taxize ----
####

# Alternate function that does the same as loop below 
# taxize::classification(specieslist, db = 'ncbi')

# loop through species and get family and order
# try() in case there are errors
for(i in 1:length(specieslist)){
  if(i == 1){
    try(res<-tax_name(query = specieslist[i], get = c("family","order"), db = "ncbi"))
  } else {
    try(res<-rbind(res,tax_name(query = specieslist[i], get = c("family","order"), db = "ncbi")))
  }
}
head(res)

# make new results dataframe
res2<-res

# get species that were not retrieved first time
sp_with_err <- setdiff(rownames(df),res$query)

# get info for species that didn't work first time
# NOTE: this could serve as first loop but probably won't need to run again
for(i in 1:length(sp_with_err)){
  
  # keeps j null each iteration
  j<-NULL 
  
  # keeps trying to get info while j is null
  while(is.null(j)) {
    try(
      j <- tax_name(query = sp_with_err[i], get = c("family","order"), db = "ncbi")
    )
  }
  
  res2<-rbind(res2,j)
  
}

# check all species are present
setdiff(rownames(df),res2$query)
setdiff(res2$query,rownames(df))

# reorder taxonomy results
res2<-res2[order(res2$query),]

# check order
table(rownames(df)==res2$query)

# reformat df
res2_out<-res2[,c(3:4)]
rownames(res2_out)<-res2$query

# NOTE: Temporary fix for error because of synonymy
# Mappia nimmoniana (formerly Nothapodytes nimmoniana)

# family
res2_out["Mappia nimmoniana",][1]<-"Icacinaceae"
# order
res2_out["Mappia nimmoniana",][2]<-"Icacinales"

# NOTE: Temporary fix for unknown error in Quintinia hyehenensis

# family
res2_out["Quintinia hyehenensis",][1]<-"Paracryphiaceae"
# order
res2_out["Quintinia hyehenensis",][2]<-"Paracryphiales"

# NOTE: temporary fix for Pandanus

# family
res2_out[grep("Pandanus",rownames(res2_out)),]$family <- "Pandanaceae"
# order
res2_out[grep("Pandanus",rownames(res2_out)),]$order <- "Pandanales"

# NOTE: temporary fix for Ternstroemia
# from POWO

# family
res2_out[grep("Ternstroemia",rownames(res2_out)),]$family <- "Pentaphylacaceae"
# order
res2_out[grep("Ternstroemia",rownames(res2_out)),]$order <- "Ericales"


# check to see if results make sense
res2_out[res2_out$order=="Poales",]

# NOT RUN: output taxonomy
saveRDS(res2_out, file = here::here("outputs/taxonomy.rds"))

# load data
res2_out<-readRDS(file = here::here("outputs/taxonomy.rds"))

# get frequencies of families/orders in our data set
length(sort(table(res2_out$family),decreasing = TRUE))
length(sort(table(res2_out$order),decreasing = TRUE))

# make data frame
ord_df <- data.frame(sort(table(res2_out$order),decreasing = TRUE))

colnames(ord_df) <- c("Order","Frequency")

write.csv(ord_df,
          row.names = FALSE,
          file = here::here("outputs/orders.csv"))

             