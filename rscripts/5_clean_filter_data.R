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

# Save filtered dataset
saveRDS(df, file = here::here("outputs/5_df_filt.rds"))


