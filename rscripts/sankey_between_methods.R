####
# Sankey plot
####
#https://r-graph-gallery.com/321-introduction-to-interactive-sankey-diagram-2.html

clust_ward <-
  readRDS(file = here::here("outputs/clust_num_k_2_7_ward.rds"))
clust_kpro <-
  readRDS(file = here::here("outputs/clust_num_k_2_7_kpro.rds"))
clust_pam <-
  readRDS(file = here::here("outputs/clust_num_k_2_7_pam.rds"))

#check orders
rownames(clust_ward) == rownames(clust_kpro)
rownames(clust_ward) == rownames(clust_pam)

#
k <- 4

# A connection data frame is a list of flows with intensity for each flow
links <-
  as.data.frame(table(as.data.frame(cbind(
    paste(clust_kpro[, k - 1], "kpro", sep = "_"),
    paste(clust_pam[, k - 1], "pam", sep = "_")
  ))))
colnames(links) <- c("source", "target", "value")

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
