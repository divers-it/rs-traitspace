rm(list=ls())

#load libraries
library(visdat)

# Read Dataset ----

df <- read.csv(
  file      = here::here("outputs/", "4_proteus_combined.csv"),
  header    = TRUE,
  row.names = 1
)

# Traits Preparation and Categorization ----

#dataset structure
str(df)

#character to factor
df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)],
                                       as.factor)

#integer to numeric
df[sapply(df, is.integer)] <- lapply(df[sapply(df, is.integer)],
                                     as.numeric)
#check structure
str(df)

#visualise missing data
vis_miss(df)
ggsave("figures/5_missing_data_pre_clean.png",
       height = 20,
       width = 40,
       bg="white",
       units = 'cm')

# Remove traits with too much missing data ----
# limit currently at 50%

#names of columns to be removed
colnames(df[ , (colSums(is.na(df)) > length(df[,1])*0.5)])

#remove columns
df <- df[ , (colSums(is.na(df)) <= length(df[,1])*0.5)]
str(df)

# Remove line with too much missing data ----
# limit currently at 50%

#names of rows to be removed
rownames(df[(rowSums(is.na(df)) > length(df[1,])*0.5), ])

#remove rows
df <- df[(rowSums(is.na(df)) <= length(df[1,])*0.5), ]
str(df)

## ----- Outlier removal -----

# Note: Removal thresholds are subjective

#remove outlier in no. structural carpels
df$Numberofstructuralcarpels[df$Numberofstructuralcarpels>998]<-NA

#remove outlier in no. ovules per carpel
df$Numberofovulesperfunctionalcarpel[df$Numberofovulesperfunctionalcarpel>998]<-NA

#remove outlier in no. ovules per carpel
df$Numberoffertilestamens[df$Numberoffertilestamens==0]<-0.0001

#remove outlier in no. ovules per carpel
df$Numberofstructuralcarpels[df$Numberofstructuralcarpels==0]<-0.0001

#visualise missing data after cleaning
vis_miss(df)
ggsave("figures/5_missing_data_post_clean.png",
       height = 20,
       width = 40,
       bg="white",
       units = 'cm')

####
# ---- Correct species names
####

#get species list
spec_list <- rownames(df)

#read in database
#The Plant List (TPL) taken from:
#https://github.com/nameMatch/Database/
# db1<-read.csv("data/Plants_TPL_database_part1.csv")
# db2<-read.csv("data/Plants_TPL_database_part2.csv")
# db3<-read.csv("data/Plants_TPL_database_part3.csv")
# db<-rbind(db1,db2,db3)

db1<-read.csv("data/Plants_WCVP_database_part1.csv")
db2<-read.csv("data/Plants_WCVP_database_part2.csv")
db3<-read.csv("data/Plants_WCVP_database_part3.csv")
db<-rbind(db1,db2,db3)


#get standardized list of taxon names
corr_list <- nameMatch(spec_list,spSource=db)

#NOTE: There are some species with issues
corr_list[corr_list$Fuzzy==1,]
corr_list[corr_list$name.dist>1,]

#reformat standardized list of names
spec_df <- corr_list[,c("Accepted_SPNAME","Genus_in_database","Family")]

#check differences
diffs<-cbind(setdiff(rownames(df),spec_df$Accepted_SPNAME),
      setdiff(spec_df$Accepted_SPNAME,rownames(df)))

colnames(diffs)<-c("old","new")

diffs

#set correct row names
rownames(df)<-spec_df$Accepted_SPNAME

#order rows based on new names
df<-df[order(rownames(df)),]

# Save filtered dataset
saveRDS(df, file = here::here("outputs/5_df_filt.rds"))


