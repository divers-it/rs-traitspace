library(rphylopic)

#load data set
df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

phylopic_uuid_list<-list()

for(i in 1:length(rownames(df))){
  
  print(paste("Getting: ",rownames(df)[i]))
  
  try(
    # Get a single image uuid
    phylopic_uuid_list[[i]] <- get_uuid(name = rownames(df)[i], n = 5), silent=TRUE
    )
  
}

names(phylopic_uuid_list)<-rownames(df)

phylopic_uuid_list[lengths(phylopic_uuid_list) != 0]

saveRDS(phylopic_uuid_list[lengths(phylopic_uuid_list) != 0],"outputs/phylopic_uuids.rds")
