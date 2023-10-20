rm(list=ls())

#load libraries
library(corrr)
library(tidyverse)
library(rcompanion)
library(harrypotter)
library(vegan)

#load one-hot data set
df<-readRDS(file = here::here("outputs/one_hot_6_df_filt_trans.rds"))

#load original data set
#df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))

# Calculate a pairwise association between all variables in a data-frame. In particular nominal vs nominal with Chi-square, numeric vs numeric with Pearson correlation, and nominal vs numeric with ANOVA.

# Adopted from https://stackoverflow.com/a/52557631/590437
mixed_assoc = function(df, cor_method="spearman", adjust_cramersv_bias=TRUE){
  df_comb = expand.grid(names(df), names(df),  stringsAsFactors = F) %>% set_names("X1", "X2")
  
  is_nominal = function(x) class(x) %in% c("factor", "character")
  # https://community.rstudio.com/t/why-is-purr-is-numeric-deprecated/3559
  # https://github.com/r-lib/rlang/issues/781
  is_numeric <- function(x) { is.integer(x) || is_double(x)}
  
  f = function(xName,yName) {
    x =  pull(df, xName)
    y =  pull(df, yName)
    
    result = if(is_nominal(x) && is_nominal(y)){
      # use bias corrected cramersV as described in https://rdrr.io/cran/rcompanion/man/cramerV.html
      cv = cramerV(as.character(x), as.character(y), bias.correct = adjust_cramersv_bias)
      data.frame(xName, yName, assoc=cv, type="cramersV")
      
    }else if(is_numeric(x) && is_numeric(y)){
      correlation = cor(x, y, method=cor_method, use="complete.obs")
      data.frame(xName, yName, assoc=correlation, type="correlation")
      
    }else if(is_numeric(x) && is_nominal(y)){
      # from https://stats.stackexchange.com/questions/119835/correlation-between-a-nominal-iv-and-a-continuous-dv-variable/124618#124618
      r_squared = summary(lm(x ~ y))$r.squared
      data.frame(xName, yName, assoc=sqrt(r_squared), type="anova")
      
    }else if(is_nominal(x) && is_numeric(y)){
      r_squared = summary(lm(y ~x))$r.squared
      data.frame(xName, yName, assoc=sqrt(r_squared), type="anova")
      
    }else {
      warning(paste("unmatched column type combination: ", class(x), class(y)))
    }
    
    # finally add complete obs number and ratio to table
    result %>% mutate(complete_obs_pairs=sum(!is.na(x) & !is.na(y)), complete_obs_ratio=complete_obs_pairs/length(x)) %>% rename(x=xName, y=yName)
  }
  
  # apply function to each variable combination
  map2_df(df_comb$X1, df_comb$X2, f)
}

#correlation matrix
cor_mat<- df %>%
  mixed_assoc() %>%
  select(x, y, assoc) %>%
  spread(y, assoc) %>%
  column_to_rownames("x") %>%
  as.matrix %>%
  as_cordf

#look at correlations
cor_mat

###
# ---- Network plot adapted from code of corrr::network_plot ----
###

#minimum correlation allowed
min_cor<-0.3
curved=FALSE

#reformat matrix
rdf <- as_matrix(cor_mat, diagonal = 1)
distance <- 1 - abs(rdf)

#CHOOSE MDS APPROACH

# #do pcoa
# points<-stats::cmdscale(distance, k = 2)
# points <- data.frame(points)
# colnames(points) <- c("x", "y")
# points$id <- rownames(points)

 #do nmds
 nmds <-
   metaMDS(distance,
           #distance = "gower",
           k = 2,
           maxit = 999, 
           trymax = 500,
           wascores = TRUE)
 points <- data.frame(nmds$points)
 colnames(points) <- c("x", "y")
 points$id <- rownames(points)

# Create a proximity matrix of the paths to be plotted.
proximity <- abs(rdf)
proximity[upper.tri(proximity)] <- NA
diag(proximity) <- NA
proximity[proximity < min_cor] <- NA

# Produce a data frame of data needed for plotting the paths.
n_paths <- sum(!is.na(proximity))
paths <- data.frame(matrix(nrow = n_paths, ncol = 6))
colnames(paths) <- c("x", "y", "xend", "yend", "proximity", "sign")
path <- 1
for (row in 1:nrow(proximity)) {
  for (col in 1:ncol(proximity)) {
    path_proximity <- proximity[row, col]
    if (!is.na(path_proximity)) {
      path_sign <- sign(rdf[row, col])
      x <- points$x[row]
      y <- points$y[row]
      xend <- points$x[col]
      yend <- points$y[col]
      paths[path, ] <- c(x, y, xend, yend, path_proximity, path_sign)
      path <- path + 1
    }
  }
}

#look at output
paths

#add point size columns
points$freq<-rep(NA,length(points$id))

for(i in 1:length(rownames(points))){
  
  points$freq[i] <- sum(na.omit(df[,points$id[i]]) != 0)
  
}


#trait type column
points$type<-rep("reproductive",length(points$id))

points$type[grep("Aqua",points$id)] <- "vegetative"
points$type[grep("Climb",points$id)] <- "vegetative"
points$type[grep("Dispersal",points$id)] <- "vegetative"
points$type[grep("Lifespan",points$id)] <- "vegetative"
points$type[grep("height",points$id)] <- "vegetative"
points$type[grep("Wood",points$id)] <- "vegetative"

head(points)

#get colours
cols<-c(rev(harrypotter::hp(2,option="Ravenclaw")),harrypotter::hp(2,option="LunaLovegood"))

#make plot
ggplot() +
  geom_point(
    data = points,
    aes(x, y, size = freq),
    shape = 21,
    fill = "#EDEDED",
    colour = "#BDBDBD",
    alpha = 0.6,
  ) + 
  geom_curve(
    data = paths,
    alpha = 0.5,
    size = 1,
    curvature = 0.2,
    aes(
      x = x,
      y = y,
      xend = xend,
      yend = yend,
      colour = proximity * sign
    )
  ) +
  # geom_segment(
  #   data = paths,
  #   alpha = 0.5,
  #   size = 1,
  #   aes(
  #     x = x,
  #     y = y,
  #     xend = xend,
  #     yend = yend,
  #     colour = proximity * sign
  #   )
  # ) + 
  geom_point(
  data = points,
  aes(x, y, fill = type),
  size=3.5,
  shape = 21,
  alpha=0.7,
  colour = "#CDCDCD"
  ) +
  geom_point(
    data = points,
    aes(x, y),
    size=1.5,
    shape = 16,
    alpha=0.7,
    colour = "black"
  ) +
  scale_fill_manual(values = c(cols[4],cols[3]), name="Trait type") +
  scale_colour_gradientn(colours = c(cols[1],"white",cols[2]), name = "Correlation") +
  scale_size(range=c(1,15), name = "Frequency",breaks = c(25,50,100,200,300)) + 
  scale_alpha(range=c(0,1)) + 
  ggrepel::geom_text_repel(
    data = points,
    aes(x, y, label = id),
    size = 3,
    segment.size = 0.0,
    segment.color = "white"
  ) + 
  theme_void() +
  guides(alpha = "none",
         size = guide_legend(override.aes = list(shape = 21, alpha=0.5, fill = "#EDEDED", colour = "#CDCDCD")),
         #shape = guide_legend(override.aes = list(shape = 21, fill = "#EDEDED", colour = "#CDCDCD")),
         #fill = guide_legend(override.aes = list(shape = 21, fill = "#EDEDED", colour = "#CDCDCD"))
         ) +
  theme(legend.key = element_blank(),
        legend.background = element_blank(),
        legend.box = "horizontal",
        legend.position = c(0.85,0.1),
        panel.background = element_rect(fill='white'))


ggsave("figures/correlation_network.png",width=12.5,height=12.5)

# ###
# # ---- Original corrr function ----
# ###
# 
# network_plot(cor_mat,colours = c("indianred2", "white", "skyblue1"),min_cor=0.3)

# ###
# # ---- igraph network plot -----
# ###
# 
# library(GGally)
# library(network)
# library(sna)
# library(ggplot2)
# library(intergraph)
# library(tidyverse)
# library(corrr)
# library(igraph)
# library(ggraph)
# 
# net = rgraph(10, mode = "graph", tprob = 0.5)
# 
# #as matrix
# cor_mat <- as.matrix(cor_mat,rownames.force=TRUE)
# rownames(cor_mat)<-cor_mat[,1]
# cor_mat <- cor_mat[,-1]
# as.numeric(cor_mat)
# 
# # Keep only high correlations
# cor_mat[cor_mat<0.25] <- 0
# 
# # Make an Igraph object from this matrix:
# network <- graph_from_adjacency_matrix( cor_mat, weighted=T, mode="undirected", diag=F)
# 
# ggraph(network) +
#   geom_edge_link(aes(edge_alpha = abs(weight), edge_width = 3, color = weight)) +
#   guides(edge_alpha = "none", edge_width = "none") +
#   scale_edge_colour_gradientn(limits = c(min(cor_mat), max(cor_mat)), colors = c("white", "dodgerblue2")) +
#   geom_node_point(color = "grey", size = 5) +
#   geom_node_text(aes(label = name), repel = TRUE) +
#   theme_graph() +
#   labs(title = "Correlations between DiveRS traits")
# 
# ###
# # ---- Marion Chartier's association plot (with correlation) ----
# ###
# 
# #load data set
# df<-readRDS(file = here::here("outputs/6_df_filt_trans.rds"))
# 
# #correlation matrix
# cor_mat<- df %>%
#   mixed_assoc() %>%
#   select(x, y, assoc) %>%
#   spread(y, assoc) %>%
#   column_to_rownames("x") %>%
#   as.matrix %>%
#   as_cordf
# 
# #as matrix
# cor_mat <- as.matrix(cor_mat,rownames.force=TRUE)
# rownames(cor_mat)<-cor_mat[,1]
# cor_mat <- cor_mat[,-1]
# as.numeric(cor_mat)
# 
# #NMDS
# cor_dist <- as.dist(cor_mat)
# 
# abs(cor_dist)
# 
# library(vegan)
# 
# myMDS <- metaMDS(cor_dist,k=5,trymax=500, wascores=F)
# 
# #####> myMDS$stress
# #####[1] 0.09147451
# x=myMDS$points[,1]
# y=myMDS$points[,2]
# 
# ##     State Occurence Cex ordre    colo
# ## 1     Tre       209 6.0     1 #9fff9f
# ## 2     Shr       150 4.0     2 #9fff9f
# ## 3    Herb        72 3.0     4 #9fff9f
# ## 4    Trop       195 4.0     5 #ffe991
# ## 5     Ari        36 2.0     6 #ffe991
# ## 6     Tem       210 4.0     7 #ffe991
# ## 7  ColPol        56 3.0     8 #ffe991
# ## 8     NAm        69 3.0     9 #e3b699
# ## 9     Eur        71 3.0    10 #e3b699
# ## 10    SAm       106 4.0    11 #e3b699
# ## 11    Afr        61 3.0    12 #e3b699
# ## 12    Ind       106 4.0    13 #e3b699
# ## 13    Aus        12 1.3    14 #e3b699
# ## 14    Ope       101 3.0    15 #c7c7ff
# ## 15    For       247 4.0    16 #c7c7ff
# ## 16    Wet        27 1.3    17 #c7c7ff
# ## 17   Lian        16 1.3    18 #9fff9f
# 
# 
# taille<-matrix(nrow=length(colnames(df)),ncol=3)
# 
# for(i in 1:length(colnames(df))){
#   
#   taille[i,1] <- colnames(df)[i]
#   
#   tmp <- df[,i]
#   
#   taille[i,2] <- sum(na.omit(tmp) != 0)
#   
# }
# 
# taille <- data.frame(taille)
# colnames(taille)<-c("State","Occurence","Cex")
# 
# taille$Occurence<-as.numeric(taille$Occurence)
# 
# taille$Cex <- log(taille$Occurence)-2
# 
# taille$ordre<-seq(1:length(rownames(taille)))
# 
# taille$colo<-rep("#CDCDCD",length(rownames(taille)))
# 
# 
# par(mar=c(4,4,1,1))
# plot(x,y, pch=16, col=as.character(taille$colo), cex=0, xlim=c(-0.5,0.4), ylim=c(-0.4,0.45))
# 
# for(a in 1:(length(x)-1)){
#   for(b in (a+1):length(x)){
#     if(cor_mat[a,b]>(-0.4) & cor_mat[a,b]<=(-0.2)){
#       lines(c(x[a],x[b]),c(y[a],y[b]),col="lightskyblue", lwd=1.5)
#     }
#     if(cor_mat[a,b]>(-0.6) & cor_mat[a,b]<=(-0.4)){
#       lines(c(x[a],x[b]),c(y[a],y[b]),col="royalblue", lwd=1.5)
#     }
#     if(cor_mat[a,b]<=(-0.6)){
#       #lines(c(x[a],x[b]),c(y[a],y[b]),col="blue4", lwd=1.5)
#     }
#     
#   }
# }
# 
# points(x,y, pch=16, col=as.character(taille$colo), cex=(taille$Cex+1.5))
# 
# for(a in 1:(length(x)-1)){
#   for(b in (a+1):length(x)){
#     if(cor_mat[a,b]>=0.6){
#       lines(c(x[a],x[b]),c(y[a],y[b]),col="darkred", lwd=1.5)
#     }
#     if(cor_mat[a,b]>=0.4 & cor_mat[a,b]<0.6){
#       lines(c(x[a],x[b]),c(y[a],y[b]),col="firebrick1", lwd=1.5)
#     }
#     if(cor_mat[a,b]>=0.2 & cor_mat[a,b]<0.4){
#       #lines(c(x[a],x[b]),c(y[a],y[b]),col="#ff8f8f", lwd=1.5)
#     }
#   }
# }
# 
# points(x,y, pch=16)
# text(x,y,label=taille$State, cex=0.7, pos=1, col="darkgreen", offset=0.2)
# 
# 
# i=-0.45
# j=0.46
# text(i,j+0.005,label="No species:",cex=0.7)
# points(i,j-0.02,pch=16, cex=2.2, col="gray90")
# text(i,j-0.02,label="<20",cex=0.7)
# points(i,j-0.05,pch=16, cex=3, col="gray90")
# text(i,j-0.05,label="20-50",cex=0.7)
# points(i,j-0.09,pch=16, cex=4, col="gray90")
# text(i,j-0.09,label="50-100",cex=0.7)
# points(i,j-0.14,pch=16, cex=5, col="gray90")
# text(i,j-0.14,label="100-200",cex=0.7)
# points(i,j-0.21,pch=16, cex=7, col="gray90")
# text(i,j-0.21,label="<200",cex=0.7)
# rm(i,j)
# 
# 
# text(0.32,-0.35,label=paste("Stress = ",round(myMDS$stress, digit=3), sep=""), cex=0.7)
