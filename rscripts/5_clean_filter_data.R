rm(list=ls())

# load libraries
library(visdat)
library(TNRS)

# read data

df <- read.csv(
  file      = here::here("outputs/", "4_proteus_combined.csv"),
  header    = TRUE,
  row.names = 1
)

# Trait preparation and categorization ----

# dataset structure
str(df)

# character to factor
df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)],
                                       as.factor)

# integer to numeric
df[sapply(df, is.integer)] <- lapply(df[sapply(df, is.integer)],
                                     as.numeric)
# check structure
str(df)

par(mar=c(5,5,5,5))

# visualise missing data
vis_miss(df) + ggplot2::theme(plot.margin = ggplot2::margin(5,40,5,5))

# ggsave("figures/5_missing_data_pre_clean.png",
#        height = 20,
#        width = 40,
#        bg="white",
#        units = 'cm')

# Remove traits with too much missing data ----
# limit currently at 50%

# names of columns to be removed
colnames(df[ , (colSums(is.na(df)) > length(df[,1])*0.5)])

# remove columns
df <- df[ , (colSums(is.na(df)) <= length(df[,1])*0.5)]
str(df)

# Remove line with too much missing data ----
# limit currently at 50%

#names of rows to be removed
rownames(df[(rowSums(is.na(df)) > length(df[1,])*0.5), ])

#remove rows
df <- df[(rowSums(is.na(df)) <= length(df[1,])*0.5), ]
str(df)

####
## Outlier removal -----
####
# NOTE: Removal thresholds are subjective based on distribution of values for each trait

#remove outlier in no. structural carpels
hist(df$Numberofstructuralcarpels)
head(sort(df$Numberofstructuralcarpels,decreasing = T))
sort(df$Numberofstructuralcarpels)
# No longer needed
# df$Numberofstructuralcarpels[df$Numberofstructuralcarpels>998]<-NA
# df$Numberofstructuralcarpels[df$Numberofstructuralcarpels==0]<-0.0001

# remove outlier in no. ovules per carpel
hist(log(df$Numberofovulesperfunctionalcarpel))
head(sort(df$Numberofovulesperfunctionalcarpel,decreasing = T),20)
df$Numberofovulesperfunctionalcarpel[df$Numberofovulesperfunctionalcarpel>998]<-NA

# remove outlier in no. fertile stamens
hist(df$Numberoffertilestamens)
head(sort(df$Numberoffertilestamens,decreasing = T))
sort(df$Numberoffertilestamens)
# No longer needed
# df$Numberoffertilestamens[df$Numberoffertilestamens==0]<-0.0001

# remove outlier in fusion of ovaries
hist(df$Fusionofovaries)
hist(asin(sqrt(df$Fusionofovaries))) #doesn't do much
head(sort(df$Fusionofovaries,decreasing = T))
sort(df$Fusionofovaries)

# remove outlier in flower size
hist(df$flowerSize)
head(sort(df$flowerSize,decreasing = T))
head(sort(df$flowerSize))

# remove outlier in seed mass
hist(df$seedMass)
head(sort(df$seedMass,decreasing = T))
head(sort(df$seedMass))

# remove outlier in height
hist(df$Maximumverticalheight)
head(sort(df$Maximumverticalheight,decreasing = T))
head(sort(df$Maximumverticalheight))

# Figure S1 ----

# visualise missing data after cleaning
vis_miss(df) + ggplot2::theme(plot.margin = ggplot2::margin(5,40,5,5))
ggsave("figures/figure_S1_missing_data_post_clean.png",
       height = 20,
       width = 40,
       bg="white",
       units = 'cm')

# Save filtered dataset
saveRDS(df, file = here::here("outputs/5_df_filt.rds"))


