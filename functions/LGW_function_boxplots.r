# boxplot function

generic_boxplot <- function(FUNC_boxplot_data, FUNC_boxplot_class, FUNC_metabolite_list){
  require(metabom8)
  require(RColorBrewer)
  require(tidyverse)
  require(plotly)
  
  bp_plotlist <- list()
  
  for(idx_metabolite in FUNC_metabolite_list){
    
    #browser()
    
    temp_plot_data <- FUNC_boxplot_data %>%
      select(all_of(FUNC_boxplot_class), all_of(idx_metabolite)) %>%
      rename(concentration = all_of(idx_metabolite), 
             plotClass = all_of(FUNC_boxplot_class))
    
    bp <- plot_ly(temp_plot_data, 
                  x = ~plotClass, 
                  y = ~concentration, 
                  color = ~plotClass, 
                  type = "box") 
      
      
      
      
      
      # 
      # ggplot(data=temp_data, 
      #            aes(x=as.factor(paste0(FUNC_boxplot_class)), 
      #                y=as.numeric(concentration))) + 		
      # geom_boxplot(aes(color = paste0(FUNC_boxplot_class), 
      #                  fill = paste0(FUNC_boxplot_class)), 
      #              outlier.shape = NA, 
      #              lwd = 1, 
      #              alpha=0.05)  +
      # geom_jitter(aes(color = paste0(FUNC_boxplot_class)), width = 0.05) +
      # labs(x = paste(""), y = paste("Concentration")) +
      # ggtitle("") +
      # theme_bw() +
      # theme(plot.title = element_text(hjust = 0.5)) +
      # theme(plot.title = element_text(size=50)) +
      # theme(axis.text.y=element_text(size = 40, margin = margin(t = 0, r = 0, b = 0, l = 2))) +
      # theme(axis.title = element_text(size = 45)) +
      # theme(axis.text.x = element_blank()) +
      # theme(legend.position = "right",
      #       legend.title = element_text(color = "black", size = 40),
      #       legend.text = element_text(color = "black", size = 40),
      #       legend.key.size = unit(2, "cm"))
    
    bp_plotlist[idx_metabolite] <- list(bp)
  }
  
  bp_plotlist
}

