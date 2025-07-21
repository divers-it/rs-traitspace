rm(list=ls())

# load libraries
library(visdat)

# Read Dataset ----

df <- read.csv(
  file      = here::here("outputs/", "one_hot_4_proteus_combined.csv"),
  header    = TRUE,
  row.names = 1
)

# Read non-one-hot dataset
df2 <- read.csv(
  file      = here::here("outputs/", "4_proteus_combined.csv"),
  header    = TRUE,
  row.names = 1
)

# check order of data sets match
rownames(df)==rownames(df2)

# dataset structure
str(df)

# integer to factor
df[sapply(df, is.integer)] <- lapply(df[sapply(df, is.integer)],
                                     as.factor)
# check structure
str(df)

# Remove traits with too much missing data ----

# names of columns to be removed
colnames(df[ , (colSums(is.na(df)) > length(df[,1])*0.5)])

# remove columns
df <- df[ , (colSums(is.na(df)) <= length(df[,1])*0.5)]
str(df)

# Remove line with too much missing data IN NORMAL (NON-ONE-HOT DATA) ----
# This is done so comparisons can be made across data sets (requires same set of species)

# names of rows to be removed
rownames(df2[ , (colSums(is.na(df2)) <= length(df2[,1])*0.5)])

# remove rows
df2 <- df2[ , (colSums(is.na(df2)) <= length(df2[,1])*0.5)]
df <- df[(rowSums(is.na(df2)) <= length(df2[1,])*0.5), ]
str(df)

#NOT RUN:
# But this might inlude some species that have more missing data than 50% after transformation
# check for species removed in one-hot but not in original 
# df_temp<-df[(rowSums(is.na(df)) <= length(df[1,])*0.5), ]
# setdiff(rownames(df),rownames(df_temp))

## Outlier removal ----

# Note: Removal thresholds are subjective

# remove outlier in no. structural carpels
df$Number.of.carpels[df$Number.of.carpels>998]<-NA

# remove outlier in no. ovules per carpel
df$Number.of.ovules.per.carpel[df$Number.of.ovules.per.carpel>998]<-NA

# remove outlier in no. ovules per carpel
df$Number.of.fertile.stamens[df$Number.of.fertile.stamens==0]<-0.0001

# remove outlier in no. ovules per carpel
df$Number.of.carpels[df$Number.of.carpels==0]<-0.0001

####
## Correct species names ----
####

# get species list
spec_list <- rownames(df)

# run TNRS to check species (best result only)
check_species <- TNRS(spec_list, matches="best", sources="wcvp")

# how many name issues?
table(check_species$Name_submitted == check_species$Accepted_species)

# species with issues
issue_species <- check_species[check_species$Name_submitted != check_species$Accepted_species,]

# examine species
issue_species[,c("Name_submitted",
                 "Accepted_name")]
# NOTE: verified on POWO: https://powo.science.kew.org/

# set correct row names
table(rownames(df) == check_species$Name_submitted)
rownames(df) <- check_species$Accepted_species

# order rows based on new names
df<-df[order(rownames(df)),]

# Save filtered dataset
saveRDS(df, file = here::here("outputs/one_hot_5_df_filt.rds"))


