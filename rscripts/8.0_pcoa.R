rm(list=ls())

# load libraries
library(dplyr)
library(ggplot2)
library(cluster)
library(patchwork)
library(wesanderson)
library(ape)
library(vegan)
library(rphylopic)

# load data set
df <- readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

# dissimilarity matrix calculation
gower_df <- daisy(df,
                  metric = "gower" )
summary(gower_df)

## Compare original and imputed data sets ----

# read in imputed data set
df_no_miss <-read.csv("outputs/imputed_with_phylo.csv",row.names=1,stringsAsFactors = TRUE)

# dissimilarity matrix calculation
gower_df_no_miss <- daisy(df_no_miss,
                          metric = "gower" )
summary(gower_df_no_miss)

# check names
table(labels(gower_df)==gsub("_"," ",labels(gower_df_no_miss)))

# check correlation
mant <- vegan::mantel(as.matrix(gower_df),as.matrix(gower_df_no_miss))
mant

# compare pairwise distances of matrices with missing data and with imputed
# png("figures/8.0_scatterplot_dist_missing_vs_imputed.png",width = 500,height = 500)
plot(gower_df,gower_df_no_miss,xlim=c(0,1),ylim=c(0,1),xlab="Original",ylab="Imputed")
abline(0,1,lty=2,col="red",lwd=2)
text(x=0.15, y=0.9, labels=paste("Mantel statistic r =",round(mant$statistic,3)))
# dev.off()

## Run and plot PCOA ----

# make into distance object
dataset_dist <- stats::as.dist(gower_df)

# run PCoA on distance matrix
dataset_pcoa <- ape::pcoa(dataset_dist)

# Recalculate relative eigenvalues by removing negative eigenvalues as in Mouillot et al.  
ev_pcoa <- dataset_pcoa$values$Eigenvalues
ev_pcoa_g0 <- ev_pcoa[ev_pcoa>0]
rel_ev_pcoa_g0 <- ev_pcoa_g0/sum(ev_pcoa_g0)

# check names
table(rownames(dataset_pcoa$vectors)==rownames(df))

# missing data per row (species) for plot
missDat <- rowSums(apply(is.na(df),2,as.numeric))/ncol(df)

# plot PCOA points on first two axes coloured by missing data
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    aes(color=missDat),
    shape=16,
    alpha=0.75,
    size=2,
    stroke = 0.5
  ) +
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

# save plot
# ggsave("figures/8.0_scatterplot_pcoa_missing.png",
#        width = 20,
#        height = 15,
#        units = 'cm')

# Make data frame of first 9 relative eigenvalues
eig_df<-data.frame(c(1:9),rel_ev_pcoa_g0[1:9])
colnames(eig_df)<-c("pcoa_axis","relative_eigenvalue")
eig_df$pcoa_axis<-as.character(eig_df$pcoa_axis)

# Proportion of variation first two axes explain
sum(rel_ev_pcoa_g0[1:2])

# plot barplot of first 9 relative eigenvalues
ggplot(eig_df, aes(x=pcoa_axis, y=relative_eigenvalue)) + 
  geom_bar(stat = "identity")
# ggsave("figures/8.0_barplot_relative_eigenvalues_pcoa.png")

####
##  Phylogenetic distance vs. trait distance ----
####

# read in phylogenetic tree
phy<-read.tree("outputs/pruned_tree.tre")

# add underscores to names
df2 <- df
rownames(df2) <- gsub(" ", "_", x = rownames(df2))

# in dataset but not in phylo
setdiff(rownames(df2),phy$tip.label)

# in phylo but not in dataset
setdiff(phy$tip.label,rownames(df2))

# get pairwise phylogenetic distances
phy_dm <- cophenetic.phylo(phy)

#dissimilarity matrix calculation
gower_df2 <- daisy(df2, metric = "gower" )
gower_df2 <- as.matrix(gower_df2)

#order rows and columns and check
phy_dm <- phy_dm[order(rownames(phy_dm)), order(colnames(phy_dm))]
gower_df2 <- gower_df2[order(rownames(gower_df2)), order(colnames(gower_df2))]
table(rownames(gower_df2)==rownames(phy_dm))
table(colnames(gower_df2)==colnames(phy_dm))

#compare distance matrices
plot(gower_df2,phy_dm,xlab="Trait distance",ylab="Phylogenetic distance")
vegan::mantel(gower_df2,phy_dm,method = "spear")

####
## Get phylopics ----
####

#choose from the different ones for a grass species
img <- pick_phylopic(name = "Oryza sativa", n = 5)
# Get a single image uuid
uuid <- get_uuid(name = "Oryza sativa", n = 1)
# Get the image for that uuid
oryza_pp <- get_phylopic(uuid = uuid)

uuid <- get_uuid(name = "Zea mays", n = 5)
zea_pp <- get_phylopic(uuid = uuid[4])

uuid <- get_uuid(name = "Phoenix dactylifera", n = 1)
phoenix_pp <- get_phylopic(uuid = uuid)

uuid <- get_uuid(name = "Solanum dulcamara", n = 2)
solanum_pp <- get_phylopic(uuid = uuid[2])

uuid <- get_uuid(name = "Plantago lanceolata", n = 2)
plantago_pp <- get_phylopic(uuid = uuid)

uuid <- get_uuid(name = "Trithuria submersa", n = 2)
trithuria_pp <- get_phylopic(uuid = uuid)

uuid <- get_uuid(name = "Commelina communis", n = 2)
commelina_pp <- get_phylopic(uuid = uuid[1])

uuid <- get_uuid(name = "Coffea arabica", n = 3)
coffea_pp <- get_phylopic(uuid = uuid[2])

uuid <- get_uuid(name = "Cornus florida", n = 1)
cornus_pp <- get_phylopic(uuid = uuid)

####
## Figure 2: Plot reproductive systems (back-engineer them first) ----
####

# load package
library(ggnewscale)

# make copy of df
df_ors <- df

### Coding classical reproductive system ----

# read in sexual system classification
ss <- read.csv("data/sexual_systems_to_verify.csv")

# simplify df and remove duplicates
ss <- data.frame(ss$NTaxDat,ss$Simple)
ss <- unique(ss)

# loop through table and recode based on sexual system and df_ors$Mating
for(i in 1:length(ss[,1])){
  
  if(is.na(ss$ss.Simple[i])){
    
    ss$ss.Simple[i] <- "unknown"
    
  }
  
  if(ss$ss.Simple[i] == "dioecious"){
    
    ss$ss.Simple[i] <- "Dioecy"
    
  }
  
  if(ss$ss.Simple[i] == "monoecious"){
    
    ss$ss.Simple[i] <- "Monoecy"
    
  }
  
  if(ss$ss.Simple[i] == "bisexual"){
    
    if(is.na(df_ors$Mating[grep(ss$ss.NTaxDat[i],rownames(df_ors))])){
      
      ss$ss.Simple[i] <- "unknown"
      
    }
    
    else if(df_ors$Mating[grep(ss$ss.NTaxDat[i],rownames(df_ors))] == "selfing"){
      
      ss$ss.Simple[i] <- "Monocliny: selfing"
      
    }
    
    else if(df_ors$Mating[grep(ss$ss.NTaxDat[i],rownames(df_ors))] == "mixed"){
      
      ss$ss.Simple[i] <- "Monocliny: mixed"
      
    }
    
    else if(df_ors$Mating[grep(ss$ss.NTaxDat[i],rownames(df_ors))] == "outcrossing"){
      
      ss$ss.Simple[i] <- "Monocliny: outcrossing"
      
    }
    
  }
  
}

#check recoding
ss

# recode flower sex, sexual system and mating system
df_ors$RS=rep("unknown",length(df_ors[,1]))
df_ors[df_ors$SexualSystem %in% "dimorphic" & df_ors$FlowerSex %in% "unisexual",]$RS="Dioecy"
df_ors[df_ors$SexualSystem %in% "monomorphic" & df_ors$FlowerSex %in% "unisexual",]$RS="Monoecy"
df_ors[df_ors$SexualSystem %in% "monomorphic" & df_ors$FlowerSex %in% "bisexual" & df_ors$Mating %in% "selfing",]$RS="Monocliny: selfing"
df_ors[df_ors$SexualSystem %in% "monomorphic" & df_ors$FlowerSex %in% "bisexual" & df_ors$Mating %in% "mixed",]$RS="Monocliny: mixed"
df_ors[df_ors$SexualSystem %in% "monomorphic" & df_ors$FlowerSex %in% "bisexual" & df_ors$Mating %in% "outcrossing",]$RS="Monocliny: outcrossing"

# reorder is close to order, but is made to change the order of the factor levels.
df_ors$RS <- factor(df_ors$RS, levels = c("Dioecy", "Monoecy", "Monocliny: selfing", "Monocliny: mixed", "Monocliny: outcrossing", "unknown"))

# combine with checked species
for(j in 1:length(ss[,1])){
  
  df_ors$RS[grep(ss$ss.NTaxDat[j],rownames(df_ors))] <- ss$ss.Simple[j]
  
}

# check RS
data.frame(rownames(df_ors),df_ors$RS)
table(df_ors$RS)

# write original reproductive systems
ors_for_csv <- data.frame(rownames(df_ors),df_ors$RS)
colnames(ors_for_csv) <- c("species", "RS")
write.csv(ors_for_csv, "outputs/original_reproductive_systems.csv")

# color palette
# colour-blind friendly
# pal 1
# one colour for dioecy, one for monoecy, then three on a gradient for bisexual selfing->outcrossing, and one for NA
# mypal <- c("#009E73","#CC79A7", colorRampPalette(c( "#E69F00","#D55E00"))(3), "black")

# pal 2
# one colour for dioecy, one for monoecy, then three on a gradient for bisexual selfing->outcrossing, and one for NA
mypal <- c("#882255","#CC6677", colorRampPalette(c( "lightblue","#332288"))(3), "black")


a1 <- ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2))+
  ggforce::geom_mark_hull(concavity = 5,expand=0,radius=0, aes(colour=as.factor(df_ors$RS), fill=as.factor(df_ors$RS),filter = df_ors$RS != 'unknown')) +
  ggforce::geom_mark_hull(concavity = 5,expand=0,radius=0, colour = "darkgrey") +
  scale_fill_manual(values=mypal) +
  scale_colour_manual(values=mypal, name = "Reproductive system") +
  geom_point(aes(colour=as.factor(df_ors$RS), shape=as.factor(df_ors$RS))) +
  scale_shape_manual(values=c(16,16,16,16,16,21)) +
  guides(fill ="none", shape = "none", color = guide_legend(override.aes = list(shape = c(16,16,16,16,16,21) ) ) ) +
  xlab(paste("PCoA Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("PCoA Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2))) +
  xlim(-0.6,0.6) + 
  ylim(-0.48,0.6) + 
  theme_bw() + theme(
    panel.border = element_blank(),
    #panel.grid.major = element_line(colour = "darkgrey"),
    #panel.grid.minor = element_line(colour = "grey"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.85, 0.85),
    legend.text = element_text(size=12),
    legend.title = element_text(size=14),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size=12),
    axis.title = element_text(size=14)) + add_phylopic(img=solanum_pp,x = 0.28, y=-0.225, ysize = 0.15,col = "grey30") +
  add_phylopic(img=commelina_pp,x = 0.44, y=0, ysize = 0.15,col = "grey30") +
  add_phylopic(img=phoenix_pp,x = -0.52, y=0, ysize = 0.15,col = "grey30") +
  add_phylopic(img=coffea_pp,x = -0.07, y=-0.38, ysize = 0.125,col = "grey30") +
  add_phylopic(img=trithuria_pp,x = 0.07, y=0.55, ysize = 0.125,col = "grey30") +
  add_phylopic(img=cornus_pp,x = -0.35, y=-0.26, ysize = 0.125,col = "grey30") +
  add_phylopic(img=zea_pp,x = -0.39, y=0.38, ysize = 0.125,col = "grey30") +
  annotate("text", x=-0.07, y=-0.46, label= "paste(italic(Coffea), ' ',italic(arabica))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=-0.52, y=-0.09, label= "paste(italic(Phoenix), ' ',italic(dactylifera))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=0.28, y=-0.32, label= "paste(italic(Solanum), ' ',italic(dulcamara))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=0.44, y=-0.1, label= "paste(italic(Commelina), ' ',italic(communis))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=0.07, y=0.47, label= "paste(italic(Trithuria), ' ',italic(submersa))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=-0.35, y=-0.34, label= "paste(italic(Cornus), ' ',italic(florida))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=-0.39, y=0.3, label= "paste(italic(Zea), ' ',italic(mays))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("segment", #phoenix
           #linetype=2,
           linewidth=0.75,
           x=-0.46,
           xend=dataset_pcoa$vectors["Phoenix dactylifera",][1], 
           y=0.03,
           yend=dataset_pcoa$vectors["Phoenix dactylifera",][2], 
           color = "black",
           alpha=0.5) + 
  annotate("segment", #solanum
           #linetype=2,
           linewidth=0.75,
           x=0.25,
           xend=dataset_pcoa$vectors["Solanum dulcamara",][1], 
           y=-0.2,
           yend=dataset_pcoa$vectors["Solanum dulcamara",][2], 
           color = "black",
           alpha=0.5) +
  annotate("segment", #coffea
           #linetype=2,
           linewidth=0.75,
           x=-0.045,
           xend=dataset_pcoa$vectors["Coffea arabica",][1], 
           y=-0.28,
           yend=dataset_pcoa$vectors["Coffea arabica",][2], 
           color = "black",
           alpha=0.5) + 
  annotate("segment", #trithuria
           #linetype=2,
           linewidth=0.75,
           x=0.07,
           xend=dataset_pcoa$vectors["Trithuria submersa",][1], 
           y=0.45,
           yend=dataset_pcoa$vectors["Trithuria submersa",][2], 
           color = "black",
           alpha=0.5) + 
  annotate("segment", #cornus
           #linetype=2,
           linewidth=0.75,
           x=-0.3,
           xend=dataset_pcoa$vectors["Cornus florida",][1],
           y=-0.21,
           yend=dataset_pcoa$vectors["Cornus florida",][2], 
           color = "black",
           alpha=0.5) +
  annotate("segment", #zea
           #linetype=2,
           linewidth=0.75,
           x=-0.35,
           xend=dataset_pcoa$vectors["Zea mays",][1],
           y=0.35,
           yend=dataset_pcoa$vectors["Zea mays",][2], 
           color = "black",
           alpha=0.5) +
  annotate("segment", #Commelina
           #linetype=2,
           linewidth=0.75,
           x=0.38,
           xend=dataset_pcoa$vectors["Commelina communis",][1],
           y=0,
           yend=dataset_pcoa$vectors["Commelina communis",][2], 
           color = "black",
           alpha=0.5)

a1

####
### Quantitative traits boxplots ----
####

# palette from scatterplot
cols<-mypal[1:5]

b1 <- ggplot(df_ors[df_ors$RS!="unknown",], aes(x=RS, y=Maximumverticalheight, fill=RS)) + 
  geom_boxplot(alpha=0.7, outlier.color=NA) + 
  geom_jitter(shape=21, position=position_jitter(0.1),alpha=0.7) + 
  scale_fill_manual(values = c(cols)) +
  # scale_y_continuous(limits = quantile(df_ors$Maximumverticalheight, c(0.025, 0.975),na.rm = TRUE)) +
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
        axis.title.y = element_text(size=16,vjust=0), #change vjust in case of patchwork issue with axis title positioning
        axis.ticks.x = element_blank(),
        # axis.ticks.y = element_blank()
  )#   + 
# annotate("text", x = -.5, y = 1.5, label = "(b)", size = 8) +
# coord_cartesian(ylim = c(-1.5, 1.75), xlim = c(1, 8), clip = "off") +
# theme(plot.margin = unit(c(1,1,1,3), "lines"))

b1

b2 <- ggplot(df_ors[df_ors$RS!="unknown",], aes(x=RS, y=flowerSize, fill=RS)) + 
  geom_boxplot(alpha=0.7, outlier.color=NA) + 
  geom_jitter(shape=21, position=position_jitter(0.1),alpha=0.7) + 
  scale_fill_manual(values = c(cols)) +
  # scale_y_continuous(limits = quantile(df_ors$flowerSize, c(0.025, 0.975),na.rm = TRUE)) +
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

# qualitative only
facts <- unlist(lapply(df_ors, is.factor))
df_temp<-df_ors[ , facts]

# change table to long form and count combinations
df_temp_melt<-reshape2::melt(df_temp,id.vars="RS")
df_temp_melt_counts <- df_temp_melt %>% group_by(RS,variable,value) %>% summarise(count=n())

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
robs <- unique(df_temp_melt_counts$RS)

## Dioecious

# subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$RS==robs[1],]

# get palette based on max counts
pal1 <- colorRampPalette(c("white",cols[1]))(max(rob_df$count))

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

## monoecious

# subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$RS==robs[2],]

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


## Monocliny: selfiing

# subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$RS==robs[3],]

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
  my_theme()

p3

## Monocliny: mixed

# subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$RS==robs[4],]

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
  my_theme() +
  theme(
    axis.text.y = element_blank()
  )

p4

## Monocliny: outcrossing

# subset df
rob_df <- df_temp_melt_counts[df_temp_melt_counts$RS==robs[5],]

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
  my_theme()

p5

### Combined plot ----

a1_fix <- a1 + coord_fixed()

patch <- ( (a1_fix / (b1 + b2 )) + plot_layout(heights=c(2, 1)) | (p1 + p2) / (p3 + p4) / (p5 + plot_spacer()) ) + plot_layout(widths=c(2, 1))
patch + plot_annotation(tag_levels = 'a',tag_prefix="(",tag_suffix=")") & theme(plot.tag = element_text(size = 14))

# save plot
ggsave("figures/figure_2_pcoa.png",
       width = 25,
       height = 20)

####
## Table S3: Correlations of traits with PCOA axes ----
####

# make empty matrix
corr_mat <- matrix(nrow=4,ncol=length(colnames(df)))

# loop through traits and calculate correlations with first four PCOA axes
for(i in 1:length(colnames(df))){
  
  if(is.numeric(df[,i])){
    
    corr_mat[1,i] <- cor(df[,i], dataset_pcoa$vectors[,1], method="pearson", use="complete.obs")
    corr_mat[2,i] <- cor(df[,i], dataset_pcoa$vectors[,2], method="pearson", use="complete.obs")
    corr_mat[3,i] <- cor(df[,i], dataset_pcoa$vectors[,3], method="pearson", use="complete.obs")
    corr_mat[4,i] <- cor(df[,i], dataset_pcoa$vectors[,4], method="pearson", use="complete.obs")
    
  } else {
    
    #https://stats.stackexchange.com/questions/119835/correlation-between-a-nominal-iv-and-a-continuous-dv-variable/124618#124618
    corr_mat[1,i] <- summary(lm(dataset_pcoa$vectors[,1] ~ df[,i]))$r.squared
    corr_mat[2,i] <- summary(lm(dataset_pcoa$vectors[,2] ~ df[,i]))$r.squared
    corr_mat[3,i] <- summary(lm(dataset_pcoa$vectors[,3] ~ df[,i]))$r.squared
    corr_mat[4,i] <- summary(lm(dataset_pcoa$vectors[,4] ~ df[,i]))$r.squared
    
  }
  
}

# set row/colnames
rownames(corr_mat)<-c("Axis.1","Axis.2","Axis.3","Axis.4")
colnames(corr_mat)<-colnames(df)

# transpose and order
corr_mat<-t(corr_mat)
corr_mat<-corr_mat[order(rownames(corr_mat)),]

# write output correlation table
write.csv(round(corr_mat,3),"outputs/8.0_correlation_pcoa_axes_traits.csv")

####
## Figure S13: Boxplots for all quantitative traits vs RS ----
####

# Maximum height
b1 <- ggplot(df_ors[df_ors$RS!="unknown",], aes(x=RS, y=Maximumverticalheight, fill=RS)) + 
  geom_boxplot(alpha=0.7, outlier.color=NA) + 
  geom_jitter(shape=21, position=position_jitter(0.1),alpha=0.7) + 
  scale_fill_manual(values = c(cols)) +
  # scale_y_continuous(limits = quantile(df_ors$Maximumverticalheight, c(0.05, 0.95),na.rm = TRUE)) +
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
        axis.text.y = element_text(size=12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14,vjust=0), #change vjust in case of patchwork issue with axis title positioning
        axis.ticks.x = element_blank(),
        # axis.ticks.y = element_blank()
  )

# Flower size
b2 <- ggplot(df_ors[df_ors$RS!="unknown",], aes(x=RS, y=flowerSize, fill=RS)) + 
  geom_boxplot(alpha=0.7, outlier.color=NA) + 
  geom_jitter(shape=21, position=position_jitter(0.1),alpha=0.7) + 
  scale_fill_manual(values = c(cols)) +
  # scale_y_continuous(limits = quantile(df_ors$flowerSize, c(0.05, 0.95),na.rm = TRUE)) +
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
        axis.title.y = element_text(size=14),
        axis.ticks.x = element_blank(),
        # axis.ticks.y = element_blank()
  )  

# Number of fertile stamens
b3 <- ggplot(df_ors[df_ors$RS!="unknown",], aes(x=RS, y=Numberoffertilestamens, fill=RS)) + 
  geom_boxplot(alpha=0.7, outlier.color=NA) + 
  geom_jitter(shape=21, position=position_jitter(0.1),alpha=0.7) + 
  scale_fill_manual(values = c(cols)) +
  # scale_y_continuous(limits = quantile(df_ors$Numberoffertilestamens, c(0.05, 0.95),na.rm = TRUE)) +
  ylab("Number of fertile stamens") +
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
        axis.title.y = element_text(size=14),
        axis.ticks.x = element_blank(),
        # axis.ticks.y = element_blank()
  )  

# Number of ovules per functional carpel
b4 <- ggplot(df_ors[df_ors$RS!="unknown",], aes(x=RS, y=Numberofovulesperfunctionalcarpel, fill=RS)) + 
  geom_boxplot(alpha=0.7, outlier.color=NA) + 
  geom_jitter(shape=21, position=position_jitter(0.1),alpha=0.7) + 
  scale_fill_manual(values = c(cols)) +
  # scale_y_continuous(limits = quantile(df_ors$Numberofovulesperfunctionalcarpel, c(0.05, 0.95),na.rm = TRUE)) +
  ylab("Number of ovules per functional carpel") +
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
        axis.title.y = element_text(size=12),
        axis.ticks.x = element_blank(),
        # axis.ticks.y = element_blank()
  )  

# Number of structural carpels
b5 <- ggplot(df_ors[df_ors$RS!="unknown",], aes(x=RS, y=Numberofstructuralcarpels, fill=RS)) + 
  geom_boxplot(alpha=0.7, outlier.color=NA) + 
  geom_jitter(shape=21, position=position_jitter(0.1),alpha=0.7) + 
  scale_fill_manual(values = c(cols)) +
  # scale_y_continuous(limits = quantile(df_ors$Numberofstructuralcarpels, c(0.05, 0.95),na.rm = TRUE)) +
  ylab("Number of structural carpels") +
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
        axis.title.y = element_text(size=14),
        axis.ticks.x = element_blank(),
        # axis.ticks.y = element_blank()
  )  

# Fusion of ovaries
b6 <- ggplot(df_ors[df_ors$RS!="unknown",], aes(x=RS, y=Fusionofovaries, fill=RS)) + 
  # geom_boxplot(alpha=0.7, outlier.color=NA) + 
  geom_jitter(shape=21, position=position_jitter(0.1),alpha=0.7) + 
  scale_fill_manual(values = c(cols)) +
  # scale_y_continuous(limits = quantile(df_ors$Fusionofovaries, c(0.05, 0.95),na.rm = TRUE)) +
  ylab("Fusion of ovaries") +
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
        axis.title.y = element_text(size=14),
        axis.ticks.x = element_blank(),
        # axis.ticks.y = element_blank()
  )  


# Seed mass
b7 <- ggplot(df_ors[df_ors$RS!="unknown",], aes(x=RS, y=seedMass, fill=RS)) + 
  geom_boxplot(alpha=0.7, outlier.color=NA) + 
  geom_jitter(shape=21, position=position_jitter(0.1),alpha=0.7) + 
  scale_fill_manual(values = c(cols)) +
  # scale_y_continuous(limits = quantile(df_ors$seedMass, c(0.05, 0.95),na.rm = TRUE)) +
  ylab("Seed mass") +
  theme(plot.margin = unit(c(1,1,1,1), "cm"),
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
        axis.title.y = element_text(size=14),
        axis.ticks.x = element_blank(),
        legend.title = element_blank()
        # axis.ticks.y = element_blank()
  )  


# combined plot
patch <- (b1 +  b2) / ( b3 + b4) / (b5 + b6 ) + ( b7 + plot_spacer())  + plot_layout(heights=c(1, 1, 1, 1),widths=c(1,1)) 
patch + plot_annotation(tag_levels = 'a',tag_prefix="(",tag_suffix=")") & theme(plot.tag = element_text(size = 14))

# save plot
ggsave("figures/figure_S13_boxplots_all_quant.png",
       width = 15,
       height = 15)


####
## Not used: PCoA with traits and density ----
####

# Edit factors for plotting to control colors
# add new factor level "none"
df$Woodiness = factor(df$Woodiness, levels=c(levels(df$Woodiness), "None"))
#df$SexualSystem = factor(df$SexualSystem, levels=c(levels(df$SexualSystem), "None"))

# convert all NA's to None
df$Woodiness[is.na(df$Woodiness)] = "None"
#df$SexualSystem = factor(df$SexualSystem, levels=c(levels(df$SexualSystem), "None"))

# y-axis scaling
# yax <- 11/19

# PCoA scatterplot with density polygons
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2)) +
  stat_density_2d(aes(fill = after_stat(level)), geom = "polygon", col=NA , n=100, bins=20) +
  scale_fill_distiller(palette = "Greys", direction = 1, guide = "none") +
  geom_point(
    aes(
      color = as.factor(df$FlowerSex),
      shape = as.factor(df$Woodiness)),
    # shape=21,
    alpha = 0.7,
    size = 2.5,
    stroke = 0.5) + 
  scale_color_manual(values=wes_palette("FantasticFox1", 3),
                     labels=c('Bisexual', 'Bisexual & Unisexual', 'Unisexual','No Data')) +
  scale_shape_discrete(labels=c('Herbaceous', 'Herbaceous & Woody', 'Woody','No Data')) +
  xlab(paste("PCoA Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("PCoA Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2))) +
  xlim(-0.6,0.6) + 
  ylim(-0.44,0.6) + 
  theme_bw() + theme(
    panel.border = element_blank(),
    #panel.grid.major = element_line(colour = "darkgrey"),
    #panel.grid.minor = element_line(colour = "grey"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.85, 0.8),
    legend.text = element_text(size=12),
    legend.title = element_text(size=14),
    axis.line = element_line(colour = "black"),
    axis.text = element_text(size=12),
    axis.title = element_text(size=14)) + 
  labs(
    colour = "Flower sex",
    shape = "Woodiness"
  ) +   
  add_phylopic(img=solanum_pp,x = 0.28, y=-0.225, ysize = 0.15,col = "grey30") +
  add_phylopic(img=commelina_pp,x = 0.44, y=0, ysize = 0.15,col = "grey30") +
  add_phylopic(img=phoenix_pp,x = -0.52, y=0, ysize = 0.15,col = "grey30") +
  add_phylopic(img=coffea_pp,x = -0.07, y=-0.35, ysize = 0.125,col = "grey30") +
  add_phylopic(img=trithuria_pp,x = 0.07, y=0.55, ysize = 0.125,col = "grey30") +
  add_phylopic(img=cornus_pp,x = -0.35, y=-0.26, ysize = 0.125,col = "grey30") +
  add_phylopic(img=zea_pp,x = -0.39, y=0.38, ysize = 0.125,col = "grey30") +
  annotate("text", x=-0.07, y=-0.43, label= "paste(italic(Coffea), ' ',italic(arabica))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=-0.52, y=-0.09, label= "paste(italic(Phoenix), ' ',italic(dactylifera))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=0.28, y=-0.32, label= "paste(italic(Solanum), ' ',italic(dulcamara))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=0.44, y=-0.1, label= "paste(italic(Commelina), ' ',italic(communis))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=0.07, y=0.47, label= "paste(italic(Trithuria), ' ',italic(submersa))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=-0.35, y=-0.34, label= "paste(italic(Cornus), ' ',italic(florida))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("text", x=-0.39, y=0.3, label= "paste(italic(Zea), ' ',italic(mays))",
           col="black", size=10 / .pt, family = "Helvetica", parse=TRUE) +
  annotate("segment", #phoenix
           #linetype=2,
           linewidth=0.75,
           x=-0.46,
           xend=dataset_pcoa$vectors["Phoenix dactylifera",][1], 
           y=0.03,
           yend=dataset_pcoa$vectors["Phoenix dactylifera",][2], 
           color = "black",
           alpha=0.5) + 
  annotate("segment", #solanum
           #linetype=2,
           linewidth=0.75,
           x=0.25,
           xend=dataset_pcoa$vectors["Solanum dulcamara",][1], 
           y=-0.2,
           yend=dataset_pcoa$vectors["Solanum dulcamara",][2], 
           color = "black",
           alpha=0.5) +
  annotate("segment", #coffea
           #linetype=2,
           linewidth=0.75,
           x=-0.045,
           xend=dataset_pcoa$vectors["Coffea arabica",][1], 
           y=-0.28,
           yend=dataset_pcoa$vectors["Coffea arabica",][2], 
           color = "black",
           alpha=0.5) + 
  annotate("segment", #trithuria
           #linetype=2,
           linewidth=0.75,
           x=0.07,
           xend=dataset_pcoa$vectors["Trithuria submersa",][1], 
           y=0.45,
           yend=dataset_pcoa$vectors["Trithuria submersa",][2], 
           color = "black",
           alpha=0.5) + 
  annotate("segment", #cornus
           #linetype=2,
           linewidth=0.75,
           x=-0.3,
           xend=dataset_pcoa$vectors["Cornus florida",][1],
           y=-0.21,
           yend=dataset_pcoa$vectors["Cornus florida",][2], 
           color = "black",
           alpha=0.5) +
  annotate("segment", #zea
           #linetype=2,
           linewidth=0.75,
           x=-0.35,
           xend=dataset_pcoa$vectors["Zea mays",][1],
           y=0.35,
           yend=dataset_pcoa$vectors["Zea mays",][2], 
           color = "black",
           alpha=0.5) +
  annotate("segment", #Commelina
           #linetype=2,
           linewidth=0.75,
           x=0.38,
           xend=dataset_pcoa$vectors["Commelina communis",][1],
           y=0,
           yend=dataset_pcoa$vectors["Commelina communis",][2], 
           color = "black",
           alpha=0.5)

# ggsave("figures/figure_2_pcoa.png",
#        width = 35,
#        #height = 35*yax,
#        height = 35,
#        units = 'cm')


####
## Not used: Plot taxonomy on PCOA ----
####

# load taxonomy
tax<-readRDS(file = here::here("outputs/taxonomy.rds"))

# check row order
table(rownames(df)==rownames(tax))

# empty column for common orders
tax$order_common<-NA

# top 10 most common orders in data
order_common<-names(sort(table(tax$order),decreasing = T)[1:8])
for(i in 1:length(order_common)){
  tax$order_common[grep(order_common[i],tax$order)]<-order_common[i]
}

# plot PCOA points on first two axes coloured by order
ggplot(data.frame(dataset_pcoa$vectors), aes(x = Axis.1, y = Axis.2)) +
  geom_point(
    aes(color=tax$order_common),
    shape=16,
    alpha=0.75,
    size=3,
    stroke = 0.5
  ) + 
  xlab(paste("Axis 1: relative eigenvalue =",round(rel_ev_pcoa_g0[1],2))) +
  ylab(paste("Axis 2: relative eigenvalue =",round(rel_ev_pcoa_g0[2],2)))

# ggsave("figures/8.0_scatterplot_pcoa_taxonomy.png",
#        width = 20,
#        height = 15,
#        units = 'cm')

