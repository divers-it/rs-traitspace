library(ggtern)

pointname <- c("Dioecious", "Monoecious","Selfer", "Mixed mating", "Bisexual\noutcrosser")
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
  labs(x = "Floral investment per flower", y = "Sexual separation", z = "Selfing rate") + 
  theme_bw() + 
  theme_hidetitles() + 
  theme_showarrows() + 
  theme_rotate() +
  theme(tern.axis.arrow = element_line(linewidth = 3)) +
  theme(axis.text=element_text(size=12))

# add points
empty_p + 
  annotate(geom  = 'text',
                   x     = fi,
                   y     = ss,
                   z     = sr,
                   label = pointname)

# palette
mypal <- c(
           colorRampPalette(c( "navyblue","lightblue"))(3)[1],
           "darkred",
           colorRampPalette(c( "navyblue","lightblue"))(3)[2],
           "darkorange",
           colorRampPalette(c( "navyblue","lightblue"))(3)[3])

# coloured labels
empty_p + 
  theme(legend.position = 'none') + 
  geom_mask() + 
  geom_label_viewport(data=df, aes(x=c(0.2,0.4,0.75,0.65,0.5),y=c(0.8,0.75,0.8,0.6,0.3),label=pointname,color=pointname, size = 3),inherit.aes = FALSE) +
  scale_color_manual(values=mypal)
  
