library(ggtern)

pointname <- c("Dioecy", "Monoecy","Monocliny:\nselfing", "Monocliny:\nmixed mating", "Monocliny:\noutcrossing")
ss <- c(90, 60, 10, 15, 15)
fi <- c(10, 20, 10, 80, 80)
sr <- c(10, 30, 90, 65, 15)

df <- data.frame(pointname,ss,fi,sr)
str(df)

# base plot
base <- ggtern(data = df, aes(ss, fi, sr))

# plot with points
base + geom_point()

# empty
empty_p <- base + 
  labs(x = "\nFloral investment per flower", y = "Sexual separation\n", z = "\nSelfing rate") + 
  theme_bw() + 
  theme_hidetitles() + 
  theme_hidelabels() + 
  theme_showarrows() + 
  theme_rotate() +
  theme(tern.axis.arrow = element_line(linewidth = 3),
        tern.axis.arrow.sep = 0.05 ) +
  theme(axis.text=element_text(size=22))

# add points
empty_p + 
  annotate(geom  = 'text',
                   x     = fi,
                   y     = ss,
                   z     = sr,
                   label = pointname)

# palette
# mypal <- c("darkred", colorRampPalette(c("navyblue", "lightblue"))(3)[2], 
#            colorRampPalette(c("navyblue", "lightblue"))(3)[1], 
#            colorRampPalette(c("navyblue", "lightblue"))(3)[3], "darkorange")

# color palette
# colour-blind friendly
# pal 1
# one colour for dioecy, one for monoecy, then three on a gradient for bisexual selfing->outcrossing, and one for NA
mypal <- c("#009E73","#CC79A7", colorRampPalette(c( "#E69F00","#D55E00"))(3), "black")

# pal 2
# one colour for dioecy, one for monoecy, then three on a gradient for bisexual selfing->outcrossing, and one for NA
# mypal <- c("#332288","#44AA99", colorRampPalette(c( "#CC6677","#882255"))(3), "black")

# coloured labels
empty_p + 
  theme(legend.position = 'none') + 
  geom_mask() + 
  geom_label_viewport(data=df,
                      hjust = "middle",
                      vjust = "center",
                      size = 6,
                      aes(x=c(0.275,0.48,0.7,0.575,0.5),
                          y=c(0.775,0.675,0.775,0.575,0.4),
                          label=pointname,
                          color=pointname),
                      inherit.aes = FALSE) +
  scale_color_manual(values=mypal[c(1,4,5,3,2)])
  
ggsave("figures/figure_5_ternary.png",
       width=9,
       height=9)
