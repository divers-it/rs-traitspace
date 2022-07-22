#' FIGURE 5


## Import Icons ----

paths <- list.files(path = here::here("data", "icons"), pattern = "*.png$", 
                    full.names = TRUE)

files <- list.files(path = here::here("data", "icons"), pattern = "*.png$", 
                    full.names = FALSE)

all_im <- lapply(paths, png::readPNG)
names(all_im) <- gsub(".png", "", files)


## Prepare Data ----

load(file = here::here("outputs", "res_for_model.RData"))

filenames <- list.files(path = here::here("outputs"), pattern = "_miss.rds$",
                        full.names = FALSE)

files <- list.files(path = here::here("outputs"), pattern = "_miss.rds$",
                    full.names = TRUE)
list_res <- lapply(files, function(x) readRDS(x))

res_for_graph_miss <- na.omit(data.frame(do.call(rbind, lapply(1:length(list_res),
                                                               function(i) {
  res <- data.frame(list_res[[i]])
  res$taxa <- gsub("_miss\\.rds", "", filenames)[i]
  return(res)
}))))


res_for_graph_miss$SP    <- NA 
res_for_graph_miss$trait <- NA 

for (i in 1:nrow(res_for_graph_miss)){ 

  res_for_graph_miss$SP[i]    <- res_for_model[rownames(res_for_model) %in% res_for_graph_miss$taxa[i], ]$S
  res_for_graph_miss$trait[i] <- res_for_model[rownames(res_for_model) %in% res_for_graph_miss$taxa[i], ]$Nb_trait
}


res_for_graph_miss$taxa <- factor(res_for_graph_miss$taxa, 
                                  levels  = unique(res_for_graph_miss$taxa[order(res_for_graph_miss$SP)]), 
                                  ordered = TRUE)

res_for_graph_miss <- res_for_graph_miss[order(res_for_graph_miss$taxa, decreasing = FALSE), ]



p2 <- ggplot(res_for_graph_miss, aes(x = miss_percent * 100, y = AUC,
                                     fill = as.factor(miss_percent * 100))) + 
  geom_boxplot() +
  theme_bw() +
  labs(x = "Traits omission (%)", 
       y = "Quality of species trait space (AUC)", size = 14) +
  harrypotter::scale_fill_hp_d(option = "ronweasley2",direction = -1)+
  facet_wrap(~ taxa, ncol = 6) + 
  theme_bw() +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  theme(strip.background = element_blank(),
        strip.text.x     = element_blank(),
        panel.grid.major = element_blank(), 
        panel.background = element_blank(),
        axis.title.x     = element_text( size=14, face="bold"),
        axis.title.y     = element_text( size=14, face="bold"),
        legend.position  = "none")

p2

grDevices::png(file = here::here("figures", "Figure5.png")) #SAVE A4

print(p2)

dev.off()
