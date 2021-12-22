# boxplot function

lgw_compare_means_ggplot_boxplot <- function(FUNC_data, 
                               FUNC_metabolite_list, 
                               FUNC_HEADER_class, 
                               FUNC_HEADER_colour,
                               FUNC_OPTION_log_plot_data,
                               FUNC_OPTION_compare_means_method
                               ){
  
  #required packages
  #require(RColorBrewer)
  #require(tidyverse)
  #require(ggpubr)
  
  ##FIRST - COMPARE MEANS
  
  bp_plotlist <- list()
  #compare_means_list() <- list()
  
  for(idx_metabolite in FUNC_metabolite_list){
   
    #browser()
    
    temp_func_data <- FUNC_data %>%
      select(all_of(FUNC_HEADER_class), all_of(FUNC_HEADER_colour), all_of(idx_metabolite)) %>%
      rename(concentration = all_of(idx_metabolite), 
             plotClass = all_of(FUNC_HEADER_class))
    
    #compare means
    compare_means_result <- compare_means(data = temp_func_data,
                                          concentration ~ plotClass,
                                          method = FUNC_OPTION_compare_means_method,
                                          p.adjust.method = "BH") %>%
      rename(feature = .y.)
    
    compare_means_result$feature[1] <- idx_metabolite
    
    #post-hoc comparison
    if(FUNC_OPTION_compare_means_method == "kruskal.test"){
      dunn_test_result <- dunn_test(data = temp_func_data,
                                    concentration ~ plotClass,
                                    p.adjust.method = "BH") 
      
      #get index of significant comparisons
      dunn_test_significance_idx <- which(dunn_test_result$p.adj.signif != "ns")
      
      #get significant comparisons for plotting later - see stat_compare_means for details
      dunn_test_comparisons <- list()
      for(idx_significance in 1:length(dunn_test_significance_idx)){
        dunn_test_comparisons[idx_significance] <- list(c(dunn_test_result$group1[dunn_test_significance_idx[idx_significance]], 
                                                       dunn_test_result$group2[dunn_test_significance_idx[idx_significance]]))
      }
      
      #combine KW and Dunn result into a single line
      dunn_test_q <- dunn_test_result$p.adj %>% as.data.frame() %>% t() %>% as_tibble()
      colnames(dunn_test_q) <- paste0(dunn_test_result$group1,
                                      " - ",
                                      dunn_test_result$group2)
      
      #create output table
      compare_means_result <- bind_cols(compare_means_result,
                                        dunn_test_q) 
      
    }
    
    #SECOND box plots
    
    if(FUNC_OPTION_log_plot_data == TRUE){
      y_axis_title <- "LOG Peak Area"
      bp <- ggplot(data=temp_func_data,
                   aes(x=as.factor(plotClass),
                       y=log(as.numeric(concentration)))
                       ) 
      bp_y_max <- temp_func_data$concentration %>% max() %>% log()
    }
    
    if(FUNC_OPTION_log_plot_data == FALSE){
      y_axis_title <- "Peak Area"
      bp <- ggplot(data=temp_func_data,
                   aes(x=as.factor(plotClass),
                       y=as.numeric(concentration))
      ) 
      bp_y_max <- temp_func_data$concentration %>% max() 
    }
    
    
    # bp <- ggplot(data=temp_plot_data,
    #              aes(x=as.factor(plotClass),
    #                  y=as.numeric(concentration))) 
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
    bp$labels$colour <- paste0(FUNC_HEADER_colour) %>% str_to_title()
    
    
  #add significance to plot
  
    bp <- bp + stat_compare_means(comparisons = dunn_test_comparisons,
                                  label = "p.signif")
    # bp <- bp + stat_compare_means(method = FUNC_OPTION_compare_means_method,
    #                               label.y = ggplot_build(bp)$layout$panel_scales_y[[1]]$range$range[2]+0.5,
    #                               label.x.npc = 0.2)
    
    
    # bp <- bp + theme(axis.text.x = element_blank()) 
     bp <- bp + theme(legend.position = "right",
            legend.title = element_text(color = "black", size = 15),
           legend.text = element_text(color = "black", size = 15))
           #legend.key.size = unit(2, "cm"))
    
    bp_plotlist[[idx_metabolite]]$BP <- bp
    bp_plotlist[[idx_metabolite]]$stats <- compare_means_result
  }
  
  bp_plotlist
  #compare_means_list
}

