require(ggplot2)
require(tidyverse)

lgw_colour_pie <- function(){

    colour_pie <- c(
    "dodgerblue2", 
    "#E31A1C", 
    "green4",
    "#6A3D9A", 
    "#FF7F00", 
    "black", 
    "white",
    "gold1",
    "skyblue2", 
    "#FB9A99", 
    "palegreen2",
    "#CAB2D6", 
    "#FDBF6F", 
    "gray70", 
    "khaki2",
    "maroon", 
    "orchid1", 
    "deeppink1", 
    "blue1", 
    "steelblue4",
    "darkturquoise", 
    "green1", 
    "yellow4", 
    "yellow3",
    "darkorange4", 
    "brown"
  )
    
    pie_data <- rep(1,26) %>% as_tibble() %>% add_column("idx" = seq(1:26))

#pie(rep(1, 26), col = colour_pie)

pie_plot <- ggplot(pie_data, aes(x="", y=as.factor(value), fill=as.factor(idx)))+
            geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start=0) +
  scale_fill_manual(values = rev(c(colour_pie))) +
  #theme_minimal() +
  geom_text(aes(y = value/2 + c(0, cumsum(value)[-length(value)]), 
                label = idx), nudge_x = 0.6, size=2)+
  labs(title = "choose colours from pie",
       x = "", y ="") +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5,
                                  size = 6),
                                  axis.line = element_blank(), 
                                  panel.grid.major = element_blank(), 
                                  panel.grid.minor = element_blank(),
                                  panel.background = element_rect(fill = "white"),
                                  axis.ticks = element_blank()
  )
  

pie_plot


}