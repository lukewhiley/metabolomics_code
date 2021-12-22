# boxplot function

lgw_ggplot_boxplot <- function(FUNC_data, 
                               FUNC_metabolite_list, 
                               FUNC_HEADER_class, 
                               FUNC_HEADER_colour,
                               FUNC_OPTION_log_plot_data,
                               FUNC_OPTION_compare_means_method
                               ){
  
  #required packages
  require(RColorBrewer)
  require(tidyverse)
  
  bp_plotlist <- list()
  
  for(idx_metabolite in FUNC_metabolite_list){
    
    browser()
    
    temp_plot_data <- FUNC_data %>%
      select(all_of(FUNC_HEADER_class), all_of(FUNC_HEADER_colour), all_of(idx_metabolite)) %>%
      rename(concentration = all_of(idx_metabolite), 
             plotClass = all_of(FUNC_HEADER_class))#,
             #plotColour = all_of(FUNC_HEADER_colour))
    
    y_axis_title <- "Peak Area"
    
    if(FUNC_OPTION_log_data == TRUE){
      bp <- ggplot(data=temp_plot_data,
                   aes(x=as.factor(plotClass),
                       y=log(as.numeric(concentration)))
                       ) 
    }
    
    if(FUNC_OPTION_log_data == FALSE){
      bp <- ggplot(data=temp_plot_data,
                   aes(x=as.factor(plotClass),
                       y=as.numeric(concentration))
      ) 
    }
    
    
    bp <- ggplot(data=temp_plot_data,
                 aes(x=as.factor(plotClass),
                     y=as.numeric(concentration))) 
    bp <- bp + geom_boxplot(color = "black",
                      aes(),
                      outlier.shape = NA,
                      #lwd = 1,
                      alpha=0.05)  +
      geom_jitter(aes(color = get(FUNC_HEADER_colour)), 
                  width = 0.01,
                  size = 3)
    bp <- bp + labs(x = paste(""), y = paste(y_axis_title))
    bp <- bp +   ggtitle(idx_metabolite)
    bp <- bp +  theme_cowplot() 
    bp <- bp +  theme(plot.title = element_text(hjust = 0.5)) 
    bp <- bp +  theme(plot.title = element_text(size=15)) 
    bp <- bp +  theme(axis.text.y=element_text(size = 15, margin = margin(t = 0, r = 0, b = 0, l = 2)))
    bp <- bp + theme(axis.title = element_text(size = 15)) 
    
    #compare means
    compare_means(data = temp_plot_data,
                  concentration ~ plotClass,
                  method = FUNC_OPTION_compare_means_method,
                  p.adjust.method = "BH") %>% View()
    
    compare_means(data = temp_plot_data,
                  concentration ~ plotClass) %>% 
      View()
    
    dunn_test(data = temp_plot_data,
              concentration ~ plotClass) %>%
      View()
    
    
    #add p values
    comparisons <- NULL
    
    
    bp <- bp + stat_compare_means()
    # bp <- bp + theme(axis.text.x = element_blank()) 
    # bp <- bp + theme(legend.position = "right",
    #        legend.title = element_text(color = "black", size = 15),
    #        legend.text = element_text(color = "black", size = 15),
    #        legend.key.size = unit(2, "cm"))
    
    bp$labels$colour <- paste0(FUNC_HEADER_colour)
    bp_plotlist[idx_metabolite] <- list(bp)
  }
  
  bp_plotlist
}

