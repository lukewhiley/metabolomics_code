require(ggplot2)
require(tidyverse)

lgw_colour_pie <- function(){

    colour_pie <- c(
      "white",
      "darkblue", 
      "darkgoldenrod2", 
      "coral3", 
      "forestgreen",
      "darkorchid4",
      "gold1",
      "steelblue2",
      "skyblue2", 
      "#FB9A99", 
      "palegreen2",
      "#CAB2D6", 
      "#FDBF6F", 
      "gray70", 
      "khaki2",
      "gray30",
      "black"
  )
    
    pie_data <- rep(1,17) %>% as_tibble() %>% add_column("idx" = seq(1:17))

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