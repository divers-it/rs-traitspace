#' FIGURE 9


## Import Icons ----

paths <- list.files(path = here::here("data", "icons"), pattern = "*.png$", 
                    full.names = TRUE)

files <- list.files(path = here::here("data", "icons"), pattern = "*.png$", 
                    full.names = FALSE)

all_im <- lapply(paths, png::readPNG)
names(all_im) <- gsub(".png", "", files)


## Prepare Data ----

load(file = here::here("outputs", "res_for_model.RData"))

filenames <- list.files(path = here::here("outputs"), pattern = "_res.rds$",
                        full.names = FALSE)

files <- list.files(path = here::here("outputs"), pattern = "_single.rds$",
                    full.names = TRUE)
list_res_single <- lapply(files, function(x) readRDS(x))

files <- list.files(path = here::here("outputs"), pattern = "_pcoa.rds$",
                    full.names = TRUE)
list_res_pcoa  <- lapply(files, function(x) readRDS(x))


res_for_graph_single <- na.omit(data.frame(do.call(rbind, lapply(1:length(list_res_single),
                                                                 function(i) {
  res_single <- data.frame(list_res_single[[i]]$cluster_core)
  
  # cluster_core = 1 === singleton
  res_single[res_single[ , 1] >  0, ] <- 2
  res_single[res_single[ , 1] == 0, ] <- 1
  res_single[res_single[ , 1] >  1, ] <- 0
  
  res_single$taxa <- gsub("_res\\.rds", "", filenames)[i]
  
  res_pcoa <- data.frame(list_res_pcoa[[i]]$vectors)[ , c(1:3)]
  
  res_cluster <- list_res_single[[i]]$cluster_core
  res <- cbind(res_single, res_pcoa, res_cluster)
  colnames(res) <- c("Single","taxa", "Pcoa1", "Pcoa2", "Pcoa3","cluster_ID")
  
  return(res)
}))))


res_for_graph_single$SP    <- NA 
res_for_graph_single$trait <- NA 

for (i in 1:nrow(res_for_graph_single)) { 
  
  res_for_graph_single$SP[i]    <- res_for_model[rownames(res_for_model) %in% res_for_graph_single$taxa[i], ]$S
  res_for_graph_single$trait[i] <- res_for_model[rownames(res_for_model) %in% res_for_graph_single$taxa[i], ]$Nb_trait
}


res_for_graph_single$taxa <- factor(res_for_graph_single$taxa, 
                                    levels  = unique(res_for_graph_single$taxa[order(res_for_graph_single$SP)]), 
                                    ordered = TRUE)
res_for_graph_single <-res_for_graph_single[order(res_for_graph_single$taxa, decreasing = FALSE),]

res_for_graph_single$Pcoa1  <- jitter(as.numeric(as.character(res_for_graph_single$Pcoa1)), factor = 50)
res_for_graph_single$Pcoa2  <- jitter(as.numeric(as.character(res_for_graph_single$Pcoa2)), factor = 50)
res_for_graph_single$Pcoa3  <- as.numeric(as.character(res_for_graph_single$Pcoa3))
res_for_graph_single$Single <- as.factor(as.character(res_for_graph_single$Single))

res_for_graph_single$cluster_ID[res_for_graph_single$cluster_ID != 1] <- 0


#order by species alphabetically
res_for_graph_single<-res_for_graph_single[order(rownames(res_for_graph_single)),]

#read data
dataset_filt<-readRDS(file = here::here("outputs/df_filt_trans.rds"))

#order filtered data
dataset_filt<-dataset_filt[order(rownames(dataset_filt)),]

#
cbind(res_for_graph_single,dataset_filt)

# cluster_core = 1 === singleton
p3 <- ggplot(res_for_graph_single, aes(x = Pcoa1, y = Pcoa2)) + 
  geom_point(aes(alpha = Single, shape = Single), size = 2) + 
  scale_shape_manual(values = c(4, 16)) +
  scale_alpha_manual(values = c(0.5, 0.8)) +
  ggalt::geom_encircle(s_shape = 1, expand = 0,size = 3, alpha = 0.7, show.legend = FALSE) +
  theme_bw() +
  labs(x = "PCoA axis 1") +
  labs(y = "PCoA axis 2") +
  #facet_wrap(~ taxa,ncol = 6, scales = "free") +
  #harrypotter::scale_colour_hp_d(option = "LunaLovegood", direction = 1) +
  theme(strip.background = element_blank(),
        strip.text.x     = element_blank(),
        strip.text.y     = element_blank(),
        panel.grid.major = element_blank(), 
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        #legend.position  = "none",
        axis.text.x      = element_blank(),   
        axis.text.y      = element_blank(),
        axis.title.x     = element_text(size = 14, face = "bold"),
        axis.title.y     = element_text(size = 14, face = "bold"),
        axis.ticks       = element_blank())

p3

hull  <- NULL
taxas <- unique(res_for_graph_single$taxa)

for (i in 1:length(taxas)) {
  sub <- res_for_graph_single[res_for_graph_single$taxa == taxas[i],]
  sub_hull <- sub[sub$cluster_ID == 1, ] %>%
    dplyr::slice(grDevices::chull(Pcoa1, Pcoa2)) 
  hull <- rbind(hull, sub_hull)
}

grDevices::png(file = here::here("figures", "Figure9.png")) #SAVE A4

print(p3)
dev.off()
