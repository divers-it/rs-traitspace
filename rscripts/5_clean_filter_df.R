#' @header *********************************************************************
#' @dataset (01) PROTEUS 2022
#' @header *********************************************************************

# Read Dataset ----

df <- read.csv(
  file      = here::here("data", "proteus_combined.csv"),
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

#identify and remove outliers
boxplot(df$Numberoffertilestamens, plot=FALSE)$out


# Remove traits with too much missing data ----
df <- df[ , (colSums(is.na(df)) < length(df[,1])*0.6)]
str(df)

# Remove line with too much missing data ----
df <- df[(rowSums(is.na(df)) < length(df[1,])*0.5), ]

# Save filtered dataset
saveRDS(df, file = here::here("outputs/df_filt.rds"))


