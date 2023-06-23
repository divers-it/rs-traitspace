#' @header *********************************************************************
#' @dataset (01) PROTEUS 2022
#' @header *********************************************************************

# Read Dataset ----

df <- read.csv(
  file      = here::here("outputs/", "proteus_combined.csv"),
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

#limit missing data
library(visdat)

png("figures/missing_data.png",height = 500,width = 1000)
vis_miss(df)
dev.off()

# Remove traits with too much missing data ----
df <- df[ , (colSums(is.na(df)) < length(df[,1])*0.5)]
str(df)

# Remove line with too much missing data ----
df <- df[(rowSums(is.na(df)) < length(df[1,])*0.5), ]

 ## ----- Outlier removal -----

## Think about removal thresholds

#remove outlier in no. structural carpels
df$Numberofstructuralcarpels[df$Numberofstructuralcarpels>998]<-NA

#remove outlier in no. ovules per carpel
df$Numberofovulesperfunctionalcarpel[df$Numberofovulesperfunctionalcarpel>998]<-NA

#remove outlier in no. ovules per carpel
df$Numberoffertilestamens[df$Numberoffertilestamens==0]<-0.0001

#remove outlier in no. ovules per carpel
df$Numberofstructuralcarpels[df$Numberofstructuralcarpels==0]<-0.0001

# Save filtered dataset
saveRDS(df, file = here::here("outputs/df_filt.rds"))


