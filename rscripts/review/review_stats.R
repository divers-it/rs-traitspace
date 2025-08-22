rm(list=ls())

# Lastly, dioecy was earlier found to be overrepresented among lianas/climbers (Renner & Ricklefs 1995). 
# You list ‘climbing’ among your traits – was there any signal?

# load data set
df <- readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

# how many climbers
table(df$Climbing)

# get climbers
df_climb <- df[df$Climbing=="climbing",]

# table of sexual system
table(df_climb$SexualSystem)

# read reproductive systems
ors <- read.csv("outputs/original_reproductive_systems.csv")

# get rs of climbers
ors_climb <- ors[ors$species %in% rownames(df_climb),]

# table of rs
table(ors_climb$RS)
