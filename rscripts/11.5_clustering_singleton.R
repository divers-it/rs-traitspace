rm(list = ls())
library(cluster)
library(ggplot2)

df <- readRDS(file = here::here("outputs/df_filt_trans_one_hot.rds"))
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
#names(dclust)="cluster"
df_dclust=merge(df,dclust,by=0)
#names(df_dclust)[21]="cluster"
#names(df_dclust)[22]="core"

###
# Plot qualitative stats of robust groups 
###
library(data.table)

#read human-readable traits
df_orig=readRDS(file = here::here("outputs/df_filt_trans.rds"))
df_orig=merge(df_orig,dclust,by=0)

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
df_labelled$cluster <- reorder(df_labelled$cluster,dclust$cluster)

#qualitative only
facts <- unlist(lapply(df_labelled, is.factor))
df_temp<-df_labelled[ , facts]

#change table to long form and count combinations
df_temp_melt<-data.table::melt(df_temp,id.vars="cluster")
df_temp_melt_counts <- df_temp_melt %>% group_by(cluster,variable,value) %>% summarise(count=n())

#add new column to remove text labels if counts are <5
df_temp_melt_counts$label<-df_temp_melt_counts$value
df_temp_melt_counts$fraction=NA
for(i in levels(df_temp_melt_counts$cluster)){
	for(j in levels(df_temp_melt_counts[df_temp_melt_counts$cluster==i,]$variable)){
		df_temp_melt_counts[df_temp_melt_counts$cluster==i & df_temp_melt_counts$variable==j,]$fraction=df_temp_melt_counts[df_temp_melt_counts$cluster==i & df_temp_melt_counts$variable==j,]$count/sum(df_temp_melt_counts[df_temp_melt_counts$cluster==i & df_temp_melt_counts$variable==j,]$count)
	}
}
df_temp_melt_counts$label[df_temp_melt_counts$fraction<0.1]<-NA
#df_temp_melt_counts$label[df_temp_melt_counts$count<3]<-NA

#make new column for text size
#df_temp_melt_counts$text_size<-df_temp_melt_counts$count^(1/2)
c52 <- sample(rainbow(52))

#plot stacked barplots per robust group for each qualitative trait, with labels
ggplot(df_temp_melt_counts, aes(variable, count, fill = value)) +
  geom_col(position = 'stack') + facet_wrap(. ~ cluster, scales = "free")  + theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    axis.title.x = element_blank(),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) + scale_fill_manual(values=c52) + geom_text(aes(label = label),
                angle = 90,
                position = position_stack(vjust = .5))

ggsave("figures/stacked_barplots_singl.pdf",width=15,height=15)


dataset_pcoa <- ape::pcoa(dataset_dist)
library(ggrepel)

#get subset of species names to highlight
sp_names=vector(length=nrow(df_dclust),mode="character")
sp_names[sp_names==0]=NA
#prune down species name to make readable
for(i in seq(1,length(dataset_single$centers))){
sp_names[dataset_single$center[i]]=df_dclust$Row.names[dataset_single$center[i]]
}

c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)
options(ggrepel.max.overlaps = Inf)
ggplot(data.frame(dataset_pcoa$vectors),
       aes(
         x = Axis.1,
         y = Axis.2,
         shape = as.factor(df_dclust$core>0),
         colour = as.factor(df_dclust$cluster)
       )) +
  geom_point(
#    color = "black",
#    shape = 21,
    alpha = 0.8,
    size = 3,
    stroke = 0.5
   ) + scale_colour_manual(values=c25) +
   geom_text_repel(aes(label = sp_names, colour = as.factor(df_dclust$cluster)), size = 3.5) + 
#							 stat_ellipse(geom = "polygon",
#                                                   aes(fill = as.factor(clust.num)),
#                                                   alpha = 0.25) +
  xlab(paste(
    "Axis 1: relative eigenvalue =",
    round(dataset_pcoa$values$Relative_eig[1], 2)
  )) +
  ylab(paste(
    "Axis 2: relative eigenvalue =",
    round(dataset_pcoa$values$Relative_eig[2], 2)
  ))

ggsave("figures/pcoa_singletons_core_halo.png",
       width = 12,
       height = 10)






