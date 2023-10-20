df<-read.csv("outputs/one_hot_imputed_with_phylo.csv", row.names = 1,stringsAsFactors = T)
str(df)

pc <- prcomp(df,
             center = TRUE,
             scale. = TRUE)
summary(pc)

library(ggbiplot)
g <- ggbiplot(pc
#              obs.scale = 1,
 #             var.scale = 1,
#              groups = training$Species,
#              ellipse = TRUE,
#              circle = TRUE,
#              ellipse.prob = 0.68
)
g <- g + scale_color_discrete(name = '')
g <- g + theme(legend.direction = 'horizontal',
               legend.position = 'top')
print(g)
