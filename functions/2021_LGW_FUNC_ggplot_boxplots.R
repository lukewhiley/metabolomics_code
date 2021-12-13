# boxplot function

lgw_ggplot_boxplot <- function(FUNC_data, 
                               FUNC_metabolite_list, 
                               FUNC_HEADER_class, 
                               FUNC_HEADER_colour,
                               FUNC_OPTION_log_data
                               ){
  
  #required packages
  require(RColorBrewer)
  require(tidyverse)
  
  bp_plotlist <- list()
  
  for(idx_metabolite in FUNC_metabolite_list){
    
    #browser()
    
    temp_plot_data <- FUNC_data %>%
      select(all_of(FUNC_HEADER_class), all_of(FUNC_HEADER_colour), all_of(idx_metabolite)) %>%
      rename(concentration = all_of(idx_metabolite), 
             plotClass = all_of(FUNC_HEADER_class))#,
             #plotColour = all_of(FUNC_HEADER_colour))
    
    y_axis_title <- "Peak Area"
    
    if(FUNC_OPTION_log_data == TRUE){
      temp_plot_data$concentration <- log(temp_plot_data$concentration+1)
      y_axis_title <- "LOG(Peak Area)"
    }
    
    
    bp <- ggplot(data=temp_plot_data,
                 aes(x=as.factor(plotClass),
                     y=as.numeric(concentration))) +
      geom_boxplot(color = "black",
                   aes(),
                   outlier.shape = NA,
                   #lwd = 1,
                   alpha=0.05)  +
      geom_jitter(aes(color = get(FUNC_HEADER_colour)), 
                  width = 0.05) +
      labs(x = paste(""), y = paste(y_axis_title)) +
      ggtitle("") +
      theme_bw() #+
      # theme(plot.title = element_text(hjust = 0.5)) +
      # theme(plot.title = element_text(size=50)) +
      # theme(axis.text.y=element_text(size = 40, margin = margin(t = 0, r = 0, b = 0, l = 2))) +
      # theme(axis.title = element_text(size = 45)) +
      # theme(axis.text.x = element_blank()) +
      # theme(legend.position = "right",
      #       legend.title = element_text(color = "black", size = 40),
      #       legend.text = element_text(color = "black", size = 40),
      #       legend.key.size = unit(2, "cm"))
    
    bp$labels$colour <- paste0(FUNC_HEADER_colour)
    bp_plotlist[idx_metabolite] <- list(bp)
  }
  
  bp_plotlist
}

