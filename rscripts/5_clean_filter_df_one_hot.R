#' @header *********************************************************************
#' @dataset (01) PROTEUS 2022
#' @header *********************************************************************

# Read Dataset ----

df <- read.csv(
  file      = here::here("outputs/", "proteus_combined_one_hot.csv"),
  header    = TRUE,
  row.names = 1
)

# Read non-onehot dataset
df2 <- read.csv(
  file      = here::here("outputs/", "proteus_combined.csv"),
  header    = TRUE,
  row.names = 1
)

rownames(df)==rownames(df2)

# Traits Preparation and Categorization ----

#dataset structure
str(df)

#character to factor
#df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)], as.factor)
#for(i in seq(1,58)){
#df[[i]]=as.factor(df[[i]])
#}

#integer to factor
df[sapply(df, is.integer)] <- lapply(df[sapply(df, is.integer)],
                                     as.factor)
#check structure
str(df)

#limit missing data
library(visdat)

pdf("figures/missing_data_one_hot.pdf")
vis_miss(df)
dev.off()

# Remove traits with too much missing data ----
df <- df[ , (colSums(is.na(df)) < length(df[,1])*0.5)]
str(df)

# Remove line with too much missing data IN NORMAL (NON-ONEHOT DATA)----
df2 <- df2[ , (colSums(is.na(df2)) < length(df2[,1])*0.5)]
df <- df[(rowSums(is.na(df2)) < length(df2[1,])*0.5), ]

## ----- Outlier removal -----

## Think about removal thresholds

#remove outlier in no. structural carpels
df$Number.of.carpels[df$Number.of.carpels>998]<-NA

#remove outlier in no. ovules per carpel
df$Number.of.ovules.per.carpel[df$Number.of.ovules.per.carpel>998]<-NA

#remove outlier in no. ovules per carpel
df$Number.of.fertile.stamens[df$Number.of.fertile.stamens==0]<-0.0001

#remove outlier in no. ovules per carpel
df$Number.of.carpels[df$Number.of.carpels==0]<-0.0001

# Save filtered dataset
saveRDS(df, file = here::here("outputs/df_filt_one_hot.rds"))


