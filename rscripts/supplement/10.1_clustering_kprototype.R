rm(list=ls())

# load packages
library(dplyr)
library(clustMixType)
library(wesanderson)
library(ggplot2)
library(cluster)
library(networkD3)
library(gridExtra)
library(patchwork)

# load formatted data
df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))


####
## K-prototypes clustering ----
####

# NOTE: Can comment out if clustering already run
# set up empty vectors
ss <- vector()
clust_memb <- vector()
kproto_list <- list()

# run clustering with different values of K up to 10
for (i in 2:10) {
  kproto_out <-
    kproto(
      df,
      k = i,
      lambda = NULL,
      iter.max = 1000,
      nstart = 10,
      na.rm = F
    )
  
  kproto_list[[i]] <- kproto_out
  
  ss[i] <- kproto_out$tot.withinss
  
  if (i == 2) {
    clust_memb <- kproto_out$cluster
  } else {
    clust_memb <- cbind(clust_memb, kproto_out$cluster)
  }
  
}

# save image as takes long to run
save.image(file = "outputs/10.1_kpro.RData")

# load previous Kproto run
load("outputs/10.1_kpro.RData")

### Select number of clusters (k) ----

# look at clustering output
head(clust_memb)

# check alignment of names
names(kproto_list[[2]]$cluster)==rownames(clust_memb)

# plot total ss to choose 
plot(ss,type='b')

# chosen value of k
kproto_out<-kproto_list[[5]]

####
## PCOA scatterplot with cluster annotation ----
####

# dissimilarity matrix calculation
gower_df <- daisy(df,
                  metric = "gower" )

# convert to distance matrix
dataset_dist <- stats::as.dist(gower_df)

# run PCOA
dataset_pcoa <- ape::pcoa(dataset_dist)

# Recalculate relative eigenvalues by removing negative eigenvalues as in Mouillot et al.  
ev_pcoa <- dataset_pcoa$values$Eigenvalues
ev_pcoa_g0 <- ev_pcoa[ev_pcoa>0]
rel_ev_pcoa_g0 <- ev_pcoa_g0/sum(ev_pcoa_g0)

# plot points on first two axes, coloured by cluster
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2, fill = as.factor(kproto_out$cluster))) +
  geom_point(
    color="black",
    shape=21,
    alpha=0.5,
    size=3,
    stroke = 0.5
  ) + 
  stat_ellipse(geom = "polygon",
               aes(fill = as.factor(kproto_out$cluster)), 
               alpha = 0.25) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2))) + 
  theme(legend.position = "none")

# save plot
# ggsave("figures/10.1_scatterplot_pcoa_kpro_k5_coloured_by_cluster.png",width = 10,height=10)

### 
# ---- Sankey plot ----
### 

# from: https://r-graph-gallery.com/321-introduction-to-interactive-sankey-diagram-2.html
# table of different k values (2-7)

for (i in 1:6) {
  if (i == 1) {
    clust.num.k.2.7 <- paste("k",i+1,"cluster",as.character(clust_memb[,i]),sep="_")
  } else {
    clust.num.k.2.7 <-
      cbind(clust.num.k.2.7, paste("k",i+1,"cluster",as.character(clust_memb[,i]),sep="_"))
  }
}

colnames(clust.num.k.2.7)<-c("2clusters",
                             "3clusters",
                             "4clusters",
                             "5clusters",
                             "6clusters",
                             "7clusters")
clust.num.k.2.7.df <-as.data.frame(clust.num.k.2.7)

rownames(clust.num.k.2.7.df)<-rownames(clust_memb)

saveRDS(clust.num.k.2.7.df, file = here::here("outputs/10.1_clust_num_k_2_7_kpro.rds"))

# A connection data frame is a list of flows with intensity for each flow

for(i in 1:(length(colnames(clust.num.k.2.7.df))-1)){
  if(i == 1){
    links<-as.data.frame(table(clust.num.k.2.7.df[,c(i,(i+1))]))
    colnames(links)<-c("source","target","value")
  } else {
    
    tmp<-as.data.frame(table(clust.num.k.2.7.df[,c(i,(i+1))]))
    colnames(tmp)<-c("source","target","value")
    links <-
      rbind(links, tmp)
  }
  
}


# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(
  name=c(as.character(links$source), 
         as.character(links$target)) %>% unique()
)

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name)-1 
links$IDtarget <- match(links$target, nodes$name)-1
links

# remove rows where values are 0
links<-links[links$value>0,]


# Make and plot the Network
p <- sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "IDsource", Target = "IDtarget",
                   Value = "value", NodeID = "name", fontSize = 14,
                   sinksRight=FALSE)
p

# save as html
saveNetwork(p, "figures/10.1_sankey_kpro.html")

####
## Identify robust groups ----
####

# make data frame of combo frequencies
combos <- as.data.frame(table(clust.num.k.2.7.df))

# remove no existant combos
combos <- combos[combos$Freq > 0, ]

# order
combos <- combos[order(combos$Freq, decreasing = T), ]

# change to strings
combos <-data.frame(lapply(combos, as.character), stringsAsFactors = FALSE)
plot(combos$Freq)

# empty list
robust<-list()

# empty vector
robust_vect_kpro<-rep(NA,length(rownames(df)))
names(robust_vect_kpro)<-rownames(df)

# loop through ordered table to extract robust groups
for(i in 1:length(combos$Freq[as.numeric(combos$Freq)>0])){
  foo<-rownames(clust.num.k.2.7.df[clust.num.k.2.7.df[, 1] == combos[i, 1] & 
                                     clust.num.k.2.7.df[, 2] == combos[i, 2] &
                                     clust.num.k.2.7.df[, 3] == combos[i, 3] &
                                     clust.num.k.2.7.df[, 4] == combos[i, 4] &
                                     clust.num.k.2.7.df[, 5] == combos[i, 5] &
                                     clust.num.k.2.7.df[, 6] == combos[i, 6],])
  
  robust[[i]]<-foo
  
  robust_vect_kpro[foo]<-i
  
}

# robust groups
robust

# keep 80% of the species in robust clusters ; others are NA
sum = 0

for (i in 1:max(robust_vect_kpro)) {
  
  if (sum > 0.62 * length(robust_vect_kpro)) {
    robust_vect_kpro[robust_vect_kpro == i] <- NA
    
  }
  
  sum = sum + length(robust[[i]])
  
}

table(robust_vect_kpro)

# complete vector of robust groups and non-robust 
robust_vect_kpro_full<-robust_vect_kpro

saveRDS(robust_vect_kpro_full, file = here::here("outputs/10.1_robust_vect_kpro_full.rds"))

#remove species not in robust groups
robust_vect_kpro<-na.omit(robust_vect_kpro)

# check order
rownames(dataset_pcoa$vectors)==rownames(clust.num.k.2.7.df)

# Plot robust groups on PCoA scatterplot
# plot points on first two axes, coloured by robust group, shaped by cluster
ggplot(
  data.frame(dataset_pcoa$vectors),
  aes(
    x = Axis.1,
    y = Axis.2,
    col = as.factor(robust_vect_kpro_full)
  )
) +
  geom_point(
    aes(shape = as.factor(clust.num.k.2.7.df$`5clusters`)),
    alpha = 0.5,
    size = 3,
    stroke = 0.5
  )

# save plot
# ggsave("figures/10.1_scatterplot_pcoa_kpro_k5_coloured_by_robust.png",width=12,height=10)

# Check amount of missing data in non-robust species
# species that dont belong to robust group
df_not_robust<-df[is.na(robust_vect_kpro_full),]
mean(is.na(df_not_robust))

# species that do belong to robust group
df_robust<-df[!is.na(robust_vect_kpro_full),]
mean(is.na(df_robust))

###
## Figure S15: Robust kpro ----
###

####
### Quantitative traits boxplots for robust clusters ----
####

# make label
robust_group<-paste("kpro_robust_",robust_vect_kpro_full,sep="")

# add label to group
df_labelled<-cbind(df,robust_group)

# palette from scatterplot
cols<-harrypotter::hp(8,option = "lunalovegood")

b1 <- ggplot(df_labelled, aes(x=robust_group, y=Maximumverticalheight, fill=robust_group)) + 
  geom_boxplot(alpha=0.7) + 
  geom_jitter(shape=21, position=position_jitter(0.1),alpha=0.7) + 
  scale_fill_manual(values = c(cols,"grey")) +
  scale_y_continuous(limits = quantile(df_labelled$Maximumverticalheight, c(0.05, 0.95),na.rm = TRUE)) +
  ylab("Maximum height") +
  theme(legend.position = "none",
        # add border 1)
        # panel.border = element_blank(),
        # color background 2)
        panel.background = element_rect(fill = "white"),
        # modify grid 3)
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color="grey70"),
        panel.grid.minor.y = element_blank(),
        # modify text, axis and colour 4) and 5)
        axis.line.y = element_line(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        axis.ticks.x = element_blank(),
        # axis.ticks.y = element_blank()
  )#   + 
# annotate("text", x = -.5, y = 1.5, label = "(b)", size = 8) +
# coord_cartesian(ylim = c(-1.5, 1.75), xlim = c(1, 8), clip = "off") +
# theme(plot.margin = unit(c(1,1,1,3), "lines"))

b1

b2 <- ggplot(df_labelled, aes(x=robust_group, y=flowerSize, fill=robust_group)) + 
  geom_boxplot(alpha=0.7) + 
  geom_jitter(shape=21, position=position_jitter(0.1),alpha=0.7) + 
  scale_fill_manual(values = c(cols,"grey")) +
  scale_y_continuous(limits = quantile(df_labelled$flowerSize, c(0.05, 0.95),na.rm = TRUE)) +
  ylab("Flower size") +
  theme(legend.position = "none",
        plot.margin = unit(c(1,1,1,1), "cm"),
        # add border 1)
        # panel.border = element_blank(),
        # color background 2)
        panel.background = element_rect(fill = "white"),
        # modify grid 3)
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color="grey70"),
        panel.grid.minor.y = element_blank(),
        # modify text, axis and colour 4) and 5)
        axis.line.y = element_line(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=16),
        axis.ticks.x = element_blank(),
        # axis.ticks.y = element_blank()
  )  

b2

b1 / b2

# ggsave("figures/10.1_robust_boxplots.png",width=15,height=10)

####
### Qualitative traits stacked barplots ----
####

# reset margins
par(mar=c(3,3,3,3))

# add group size to robust group label
for (i in 1:length(unique(df_labelled$robust_group))) {
  df_labelled$robust_group[df_labelled$robust_group %in% sort(unique(df_labelled$robust_group))[i]] <-
    paste(
      sort(unique(df_labelled$robust_group))[i],
      " (n = ",
      table(df_labelled$robust_group)[i],
      ")",
      sep = ""
    )
  
}

# make as factor for grouping
df_labelled$robust_group<-as.factor(df_labelled$robust_group)

# qualitative only
facts <- unlist(lapply(df_labelled, is.factor))
df_temp<-df_labelled[ , facts]

# change table to long form and count combinations
df_temp_melt<-reshape2::melt(df_temp,id.vars="robust_group")
df_temp_melt_counts <- df_temp_melt %>% group_by(robust_group,variable,value) %>% summarise(count=n())

# add new column to remove text labels if counts are <5
df_temp_melt_counts$label<-df_temp_melt_counts$value
df_temp_melt_counts$label[df_temp_melt_counts$count<3]<-NA

# theme
my_theme <- function() {
  theme(
    # add border 1)
    panel.border = element_blank(),
    # color background 2)
    panel.background = element_rect(fill = "white"),
    # modify grid 3)
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    # modify text, axis and colour 4) and 5)
    axis.text.x = element_text(size=16),
    axis.text.y = element_text(size=16),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_line(),
    axis.ticks.length=unit(.25, "cm"),
    axis.ticks.y = element_blank(),
    # legend
    legend.position = "none",
    # margin
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
  )
}

# get names of robust groups
robs<-unique(df_temp_melt_counts$robust_group)

# Robust group 1
cols<-harrypotter::hp(8,option = "lunalovegood")

# subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$robust_group==robs[1],]

# get palette based on max counts
pal1<-colorRampPalette(c("white",cols[1]))(max(rob_df$count))

# make new column for colours based on count
rob_df$pal <- pal1[rob_df$count]
rob_df$pal[is.na(rob_df$value)]<-"grey20"

p1 <- ggplot(rob_df, aes(x = variable, y = count, fill=pal, alpha=0.98)) +
  geom_bar(position="stack", stat="identity",col="black") +
  scale_fill_identity() + 
  geom_text(aes(size = count, label = label), position = position_stack(vjust = 0.5)) + 
  coord_flip() +
  my_theme()

p1

## Robust group 2

# subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$robust_group==robs[2],]

# get palette based on max counts
pal2<-colorRampPalette(c("white",cols[2]))(max(rob_df$count))

# make new column for colours based on count
rob_df$pal <- pal2[rob_df$count]
rob_df$pal[is.na(rob_df$value)]<-"grey20"

p2 <- ggplot(rob_df, aes(x = variable, y = count, fill=pal, alpha=0.98)) +
  geom_bar(position="stack", stat="identity",col="black") +
  scale_fill_identity() + 
  geom_text(aes(size = count, label = label), position = position_stack(vjust = 0.5)) + 
  coord_flip() +
  my_theme() +
  theme(
    axis.text.y = element_blank()
  )

p2


## Robust group 3

# subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$robust_group==robs[3],]

# get palette based on max counts
pal3<-colorRampPalette(c("white",cols[3]))(max(rob_df$count))

# make new column for colours based on count
rob_df$pal <- pal3[rob_df$count]
rob_df$pal[is.na(rob_df$value)]<-"grey20"

p3 <- ggplot(rob_df, aes(x = variable, y = count, fill=pal, alpha=0.98)) +
  geom_bar(position="stack", stat="identity",col="black") +
  scale_fill_identity() + 
  geom_text(aes(size = count, label = label), position = position_stack(vjust = 0.5)) + 
  coord_flip() +
  my_theme() +
  theme(
    axis.text.y = element_blank()
  )

p3

## Robust group 4

# subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$robust_group==robs[4],]

# get palette based on max counts
pal3<-colorRampPalette(c("white",cols[4]))(max(rob_df$count))

# make new column for colours based on count
rob_df$pal <- pal3[rob_df$count]
rob_df$pal[is.na(rob_df$value)]<-"grey20"

p4 <- ggplot(rob_df, aes(x = variable, y = count, fill=pal, alpha=0.98)) +
  geom_bar(position="stack", stat="identity",col="black") +
  scale_fill_identity() + 
  geom_text(aes(size = count, label = label), position = position_stack(vjust = 0.5)) + 
  coord_flip() +
  my_theme()

p4

## Robust group 5

# subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$robust_group==robs[5],]

# get palette based on max counts
pal3<-colorRampPalette(c("white",cols[5]))(max(rob_df$count))

# make new column for colours based on count
rob_df$pal <- pal3[rob_df$count]
rob_df$pal[is.na(rob_df$value)]<-"grey20"

p5 <- ggplot(rob_df, aes(x = variable, y = count, fill=pal, alpha=0.98)) +
  geom_bar(position="stack", stat="identity",col="black") +
  scale_fill_identity() + 
  geom_text(aes(size = count, label = label), position = position_stack(vjust = 0.5)) + 
  coord_flip() +
  my_theme() +
  theme(
    axis.text.y = element_blank()
  )

p5

## Robust group 6

# subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$robust_group==robs[6],]

# get palette based on max counts
pal3<-colorRampPalette(c("white",cols[6]))(max(rob_df$count))

# make new column for colours based on count
rob_df$pal <- pal3[rob_df$count]
rob_df$pal[is.na(rob_df$value)]<-"grey20"

p6 <- ggplot(rob_df, aes(x = variable, y = count, fill=pal, alpha=0.98)) +
  geom_bar(position="stack", stat="identity",col="black") +
  scale_fill_identity() + 
  geom_text(aes(size = count, label = label), position = position_stack(vjust = 0.5)) + 
  coord_flip() +
  my_theme() +
  theme(
    axis.text.y = element_blank()
  )

p6

## Robust group 7

# subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$robust_group==robs[7],]

# get palette based on max counts
pal3<-colorRampPalette(c("white",cols[7]))(max(rob_df$count))

# make new column for colours based on count
rob_df$pal <- pal3[rob_df$count]
rob_df$pal[is.na(rob_df$value)]<-"grey20"

p7 <- ggplot(rob_df, aes(x = variable, y = count, fill=pal, alpha=0.98)) +
  geom_bar(position="stack", stat="identity",col="black") +
  scale_fill_identity() + 
  geom_text(aes(size = count, label = label), position = position_stack(vjust = 0.5)) + 
  coord_flip() +
  my_theme()

p7

## Robust group 8

# subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$robust_group==robs[8],]

# get palette based on max counts
pal3<-colorRampPalette(c("white","grey"))(max(rob_df$count))

# make new column for colours based on count
rob_df$pal <- pal3[rob_df$count]
rob_df$pal[is.na(rob_df$value)]<-"grey20"

p8 <- ggplot(rob_df, aes(x = variable, y = count, fill=pal, alpha=0.98)) +
  geom_bar(position="stack", stat="identity",col="black") +
  scale_fill_identity() + 
  geom_text(aes(size = count, label = label), position = position_stack(vjust = 0.5)) + 
  coord_flip() +
  my_theme() +
  theme(
    axis.text.y = element_blank()
  )

p8

(p1 + p2 + p3) / (p4 + p5 + p6) / (p7 + p8 + plot_spacer() )

### Combined plot ----

patch <- ( p1 + p2 + p3 ) / (p4 + p5 + p6) / (p7 + p8 + plot_spacer()) / ( b1 + b2 ) + plot_layout(heights=c(1, 1 , 1, 1))
patch + plot_annotation(tag_levels = 'a',tag_prefix="(",tag_suffix=")") & theme(plot.tag = element_text(size = 14))

# save plot
ggsave("figures/figure_S15_robust_kpro.png",width=20,height=25)


