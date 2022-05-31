lgw_colour_pie_select <- function(FUNC_OPTION_colour_choice,
                                  FUNC_OPTION_pie_label){
  
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
  
  pie_data <- rep(1,length(FUNC_OPTION_colour_choice)) %>% 
    as_tibble() %>% 
    add_column("idx" = seq(1:length(FUNC_OPTION_colour_choice)))
  
  colour_pie_selected <- colour_pie[c(FUNC_OPTION_colour_choice)]
  
  
  ##
  pie_plot <- ggplot(pie_data, aes(x="", y=as.factor(value), fill=as.factor(idx)))+
    geom_bar(width = 1, stat = "identity") +
    coord_polar("y", start=0) +
    scale_fill_manual(values = rev(c(colour_pie_selected))) +
    theme_minimal() +
    geom_text(aes(y = value/2 + c(0, cumsum(value)[-length(value)]), 
                  label = idx), nudge_x = 0.75,  size=0.75)+
    labs(title = "Project Colours",
         x = "", y ="") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5,
                                    vjust= -12,
                                    size = 15),
          axis.line = element_blank(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "lightgrey"),
          axis.ticks = element_blank(),
          plot.margin = unit(c(0,0,0,0), "mm")
          ) +
    geom_text(aes(label=rev(paste0(FUNC_OPTION_pie_label))),
              position = position_stack(vjust=0.5))
  
  #+
   # coord_fixed(ratio = 1) 
    
  
  
  pie_plot_out <- list()
  pie_plot_out$pie_plot <- pie_plot
  pie_plot_out$colour_selection <- colour_pie_selected
  
  pie_plot_out

  #pie(rep(1, length(FUNC_OPTION_colour_choice)), col = colour_pie[FUNC_OPTION_colour_choice], main = "Project colours")
  
 
  #paste0(colour_pie[FUNC_OPTION_colour_choice])
  
}