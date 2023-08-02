####
# Sankey plot
####

#https://r-graph-gallery.com/321-introduction-to-interactive-sankey-diagram-2.html

clust_ward_one_hot=readRDS(file = here::here("outputs/one_hot_11_clust_num_k_2_7_ward.rds"))
clust_ward=readRDS(file = here::here("outputs/11_clust_num_k_2_7_ward.rds"))
clust_kpro_one_hot=readRDS(file = here::here("outputs/one_hot_11.1_clust_num_k_2_7_kpro.rds"))
clust_kpro=readRDS(file = here::here("outputs/11.1_clust_num_k_2_7_kpro.rds"))
clust_pam_one_hot=readRDS(file = here::here("outputs/one_hot_11.2_clust_num_k_2_7_pam.rds"))
clust_pam=readRDS(file = here::here("outputs/11.2_clust_num_k_2_7_pam.rds"))
clust_density_one_hot=readRDS(file = here::here("outputs/one_hot_11.4_clust_num_k_2_7_density.rds"))
clust_density=readRDS(file = here::here("outputs/11.4_clust_num_k_2_7_density.rds"))

#check orders
rownames(clust_ward) == rownames(clust_density)
rownames(clust_ward) == rownames(clust_kpro)
rownames(clust_ward) == rownames(clust_pam)
rownames(clust_ward) == rownames(clust_ward_one_hot)
rownames(clust_ward) == rownames(clust_density_one_hot)
rownames(clust_ward) == rownames(clust_kpro_one_hot)
rownames(clust_ward) == rownames(clust_pam_one_hot)

#

nonrobust=data.frame(k=seq(2,7),n=NA)

for(k in seq(2,7)){

j=1
for (i in ls()[grep("clust_",ls())]) {
  method=sub('clust_', '', i)
  if (j == 1) {
    clust.num.methods <- data.frame(vect=paste(method,mget(i)[[1]][,paste(k,"clusters",sep="")],sep="_"))
    names(clust.num.methods)[names(clust.num.methods) == "vect"] <- paste(i,k,sep="_")
  } else {
    vect=paste(i,k,sep="_")
    clust.num.methods <-
      cbind(clust.num.methods, data.frame(vect=paste(method,mget(i)[[1]][,paste(k,"clusters",sep="")],sep="_")))
    names(clust.num.methods)[names(clust.num.methods) == "vect"] <- paste(i,k,sep="_")
  }
  j=j+1
}
rownames(clust.num.methods)=rownames(clust_pam_one_hot)
# A connection data frame is a list of flows with intensity for each flow
for(i in 1:(length(colnames(clust.num.methods))-1)){
  if(i == 1){
    links<-as.data.frame(table(clust.num.methods[,c(i,(i+1))]))
    colnames(links)<-c("source","target","value")
  } else {
    
    tmp<-as.data.frame(table(clust.num.methods[,c(i,(i+1))]))
    colnames(tmp)<-c("source","target","value")
    links <-
      rbind(links, tmp)
  }
  
}

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes <- data.frame(name = c(as.character(links$source),
                             as.character(links$target)) %>% unique())

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
links$IDsource <- match(links$source, nodes$name) - 1
links$IDtarget <- match(links$target, nodes$name) - 1

links

#remove rows where values are 0
links <- links[links$value > 0, ]

# Library
library(networkD3)
library(dplyr)
# Make the Network
p <- sankeyNetwork(
  Links = links,
  Nodes = nodes,
  Source = "IDsource",
  Target = "IDtarget",
  Value = "value",
  NodeID = "name",
  sinksRight = FALSE
)

p

saveNetwork(p, paste("figures/sankey_between_all_k",k,".html",sep=""))

###
# ---- Robust combinations ----
###

#make data frame of combo frequencies
combos <- as.data.frame(table(clust.num.methods))

#remove no existant combos
combos <- combos[combos$Freq > 0, ]

#order
combos <- combos[order(combos$Freq, decreasing = T), ]

#change to strings
combos <-data.frame(lapply(combos, as.character), stringsAsFactors = FALSE)

#table of different combinations and their frequencies
head(combos)

#empty list
robust<-list()

#empty vector
robust_vect<-rep(NA,length(rownames(clust_ward)))
names(robust_vect)<-rownames(clust_ward)

#loop through ordered table to extract robust groups
#change value in loop for threshold
for(i in 1:length(combos$Freq[as.numeric(combos$Freq)>360/(10*k)])){
  foo<-rownames(clust.num.methods[clust.num.methods[, 1] == combos[i, 1] & 
                                                clust.num.methods[, 2] == combos[i, 2] &
                                                clust.num.methods[, 3] == combos[i, 3] &
                                                clust.num.methods[, 4] == combos[i, 4] &
                                                clust.num.methods[, 5] == combos[i, 5] &
                                                clust.num.methods[, 6] == combos[i, 6] &
                                                clust.num.methods[, 7] == combos[i, 7] &
                                                clust.num.methods[, 8] == combos[i, 8],])
  
  #add names of each robust group in list
  robust[[i]]<-foo
  
  #CHECK - adds number to species
  robust_vect[foo]<-i
  
}


#robust groups
robust

#complete vector of robust groups and non-robust 
robust_vect_full<-robust_vect
#saveRDS(robust_vect_full, file = here::here("outputs/11.1_robust_vect_full.rds"))
nonrobust[nonrobust$k==k,]$n=sum(is.na(robust_vect))
}
plot(nonrobust$k,nonrobust$n,type="b")

