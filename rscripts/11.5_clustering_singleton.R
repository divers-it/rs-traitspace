rm(list = ls())
library(cluster)
library(ggplot)

df <- readRDS(file = here::here("outputs/df_filt_trans.rds"))
gower_df <- daisy(df,
                  metric = "gower" )
dataset_dist <- stats::as.dist(gower_df)
class(dataset_dist)="dist"

source("R/singleton.R")

dataset_single <- DPC(dataset_dist,
  metric         = "predefined",
  radius         = "automatic",
  density_filter = "global_threshold",
  rho_threshold  = 1
)
for(i in seq(1,12)){
print(dataset_single$rho[dataset_single$center[i],])
}

dclust=data.frame(cluster=dataset_single$cluster,core=dataset_single$cluster_core)
rownames(dclust)=rownames(dataset_single$rho)
names(dclust)="cluster"
df_dclust=merge(df,dclust,by=0)
names(df_dclust)[21]="cluster"

###
# Plot qualitative stats of robust groups
###
library(data.table)

df_labelled=df_dclust

df_labelled$cluster<-as.character(df_labelled$cluster)

#add group size to robust group label
for (i in 1:length(unique(df_labelled$cluster))) {
  df_labelled$cluster[df_labelled$cluster %in% sort(unique(df_labelled$cluster))[i]] <-
    paste(
      sort(unique(df_labelled$cluster))[i],
      " (n = ",
      table(df_labelled$cluster)[i],
      ")",
      sep = ""
    )
  
}

#make as factor for grouping
df_labelled$cluster<-as.factor(df_labelled$cluster)

#qualitative only
facts <- unlist(lapply(df_labelled, is.factor))
df_temp<-df_labelled[ , facts]

#change table to long form and count combinations
df_temp_melt<-data.table::melt(df_temp,id.vars="cluster")
df_temp_melt_counts <- df_temp_melt %>% group_by(cluster,variable,value) %>% summarise(count=n())

#add new column to remove text labels if counts are <5
df_temp_melt_counts$label<-df_temp_melt_counts$value
df_temp_melt_counts$label[df_temp_melt_counts$count<3]<-NA

#make new column for text size
#df_temp_melt_counts$text_size<-df_temp_melt_counts$count^(1/2)

#plot stacked barplots per robust group for each qualitative trait, with labels
ggplot(df_temp_melt_counts, aes(variable, count, fill = value)) +
  geom_col(position = 'stack') + facet_wrap(. ~ cluster, scales = "free")  + theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.title.x = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) + geom_text(aes(size = count,label = label),
                angle = 90,
                position = position_stack(vjust = .5))

ggsave("figures/stacked_barplots_singl.pdf",width=15,height=15)


