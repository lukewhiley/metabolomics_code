# boxplot function

lgw_compare_means_ggplot_boxplot <- function(FUNC_data, 
                               FUNC_metabolite_list, 
                               FUNC_class_to_include,
                               FUNC_HEADER_class, 
                               FUNC_HEADER_colour,
                               FUNC_OPTION_colour_choice,
                               FUNC_OPTION_log_plot_data,
                               FUNC_OPTION_compare_means_method,
                               FUNC_HEADER_paired,
                               FUNC_plot_comparisons,
                               FUNC_OPTION_plot_qc,
                               FUNC_OPTION_point_size,
                               FUNC_OPTION_plot_show_legend,
                               FUNC_OPTION_plot_outliers
                               ){
  ##FIRST - COMPARE MEANS
  bp_plotlist <- list()
  bp_plotlist$bp <- list()
  bp_plotlist$table <- list()
  
  for(idx_metabolite in FUNC_metabolite_list){
 #print(idx_metabolite)
    temp_func_data <- FUNC_data %>%
      select(all_of(FUNC_HEADER_class), all_of(FUNC_HEADER_colour), all_of(idx_metabolite)) %>%
      filter(!!as.symbol(FUNC_HEADER_class) %in% FUNC_class_to_include) %>%
      rename(concentration = all_of(idx_metabolite), 
             plot_class = all_of(FUNC_HEADER_class))
    
    
    if(!is.na(FUNC_HEADER_paired)){
    temp_func_data <- FUNC_data %>%
      select(all_of(FUNC_HEADER_class), 
             all_of(FUNC_HEADER_colour), 
             all_of(idx_metabolite),
             all_of(FUNC_HEADER_paired)) %>%
      rename(concentration = all_of(idx_metabolite), 
             plot_class = all_of(FUNC_HEADER_class),
             pair_group = all_of(FUNC_HEADER_paired))
    }
    
  
    # if function to naming fix error when FUNC_HEADER_class and FUNC_HEADER_colour are the same column
    FUNC_HEADER_temp_colour <- FUNC_HEADER_colour
    
    if(FUNC_HEADER_class == FUNC_HEADER_colour){
      FUNC_HEADER_temp_colour <- "plot_class"
    }
    
    
    #run kruskal.test on everything first - if only 2x groups returns a wilcox.test result anyway.
    
    #compare means
    compare_means_result <- compare_means(data = temp_func_data %>%
                                            filter(plot_class != "qc") %>%
                                            filter(plot_class != "QC"),
                                          concentration ~ plot_class,
                                          method = "kruskal.test",
                                          p.adjust.method = "BH") %>%
      rename(feature = .y.)
    
    compare_means_result$feature[1] <- idx_metabolite
    
    #compare_means_result$p <- compare_means_result$p %>% signif(digits = 3)
    
    ##### KRUSKAL.TEST ###########
    
    #post-hoc comparison
    if(FUNC_OPTION_compare_means_method == "kruskal.test"){
      dunn_test_result <- dunn_test(data = temp_func_data %>%
                                      filter(plot_class != "qc") %>%
                                      filter(plot_class != "QC"),
                                    concentration ~ plot_class,
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
      dunn_test_q <- dunn_test_result$p.adj %>% as.data.frame() %>% t() %>% as_tibble() %>% signif(digits = 2)
      colnames(dunn_test_q) <- paste0(dunn_test_result$group1,
                                      " - ",
                                      dunn_test_result$group2)
      
      #create output table
      compare_means_result <- bind_cols(compare_means_result,
                                        dunn_test_q) 
      
    }
  
    #####WILCOX.TEST##########
    #complete wilcox test
    if(FUNC_OPTION_compare_means_method == "wilcox.test"){
     wilcox_test_result <- compare_means(data = temp_func_data %>%
                      filter(plot_class != "qc") %>%
                      filter(plot_class != "QC"),
                      concentration ~ plot_class,
                      method = "wilcox.test",
                      p.adjust.method = "BH") %>%
        rename(feature = .y.)
     
     wilcox_test_result$comparison <- paste0(wilcox_test_result$group1, 
                                             " - ",
                                             wilcox_test_result$group2)
     
     # wilcox_test_result <- wilcox_test_result %>%
     #   filter(comparison %in% FUNC_plot_comparisons)
     
     #change name of object to match code to dunn test (above)
     dunn_test_result <- wilcox_test_result 
     
     #get index of significant comparisons
     dunn_test_significance_idx <- which(dunn_test_result$p.signif != "ns")
     
     #get significant comparisons for plotting later - see stat_compare_means for details
     dunn_test_comparisons <- list()
     for(idx_significance in 1:length(dunn_test_significance_idx)){
       dunn_test_comparisons[idx_significance] <- list(c(dunn_test_result$group1[dunn_test_significance_idx[idx_significance]], 
                                                         dunn_test_result$group2[dunn_test_significance_idx[idx_significance]]))
     }
     
     wilcox_test_result_t <- wilcox_test_result %>% 
       select(p, comparison) %>%
       t() %>% 
       as_tibble()
     wilcox_test_result_t <- setNames(wilcox_test_result_t, wilcox_test_result$comparison)
     wilcox_test_result_t <- wilcox_test_result_t[1,]
     wilcox_test_result_t <- wilcox_test_result_t %>%
       select(all_of(FUNC_plot_comparisons))
     
     wilcox_test_result_t_num <-wilcox_test_result_t %>% 
       as.numeric() %>% 
       as_tibble() %>% 
       t() %>% 
       as_tibble() %>% 
       setNames(names(wilcox_test_result_t)) 

     
     #create output table
     compare_means_result <- bind_cols(compare_means_result,
                                       wilcox_test_result_t_num)
     
     if(FUNC_OPTION_compare_means_method == "wilcox.test"){
       compare_means_result$method <- "wilcox.test"
     }
      
    }
    
    #####WILCOX.TEST.PAIRED##########

    #complete pairwise paired wilcox test
    if(FUNC_OPTION_compare_means_method == "wilcox.test.paired"){
      
      #create loop to perform pairwise tests
      wilcox_test_result <- NULL
      for(idx_str_comparison in FUNC_plot_comparisons){
        loop_comparisons <- c(gsub(" - .*","", idx_str_comparison),
                              gsub(".* - ","", idx_str_comparison)
                              )
        #create sub-frame for paired data
        temp_func_paired_data <- temp_func_data %>%
          filter(plot_class %in% loop_comparisons)
        
        temp_func_paired_data <- temp_func_paired_data %>%
          slice(which(duplicated(temp_func_paired_data$pair_group)|duplicated(temp_func_paired_data$pair_group, fromLast = TRUE))) %>%
          arrange(pair_group)
        
        #perform pairwise paired test
      loop_wilcox_test_result <-  compare_means(data = temp_func_paired_data %>%
                                             filter(plot_class != "qc") %>%
                                             filter(plot_class != "QC"),
                                           concentration ~ plot_class,
                                           method = "wilcox.test",
                                           paired = TRUE,
                                           p.adjust.method = "BH") %>%
        rename(feature = .y.)
      
    
      #create column of comparison made
      loop_wilcox_test_result$comparison <- paste0(loop_wilcox_test_result$group1, 
                                              " - ",
                                              loop_wilcox_test_result$group2)
      
      wilcox_test_result <- bind_rows(wilcox_test_result,
                                      loop_wilcox_test_result)
    }
      
      wilcox_test_result <- wilcox_test_result %>%
        filter(comparison %in% FUNC_plot_comparisons)
      
      #change name of object to match code to dunn test (above)
      dunn_test_result <- wilcox_test_result 
      
      #get index of significant comparisons
      dunn_test_significance_idx <- which(dunn_test_result$p.signif != "ns")
      
      #get significant comparisons for plotting later - see stat_compare_means for details
      dunn_test_comparisons <- list()
      for(idx_significance in 1:length(dunn_test_significance_idx)){
        dunn_test_comparisons[idx_significance] <- list(c(dunn_test_result$group1[dunn_test_significance_idx[idx_significance]], 
                                                          dunn_test_result$group2[dunn_test_significance_idx[idx_significance]]))
      }
      
      #output pairwise result
      wilcox_test_result_t <- wilcox_test_result %>% 
        select(p, comparison) %>%
        t() %>% 
        as_tibble()
      wilcox_test_result_t <- setNames(wilcox_test_result_t, wilcox_test_result$comparison)
      wilcox_test_result_t <- wilcox_test_result_t[1,]
      wilcox_test_result_t <- wilcox_test_result_t %>%
        select(all_of(FUNC_plot_comparisons))
      
      wilcox_test_result_t_num <-wilcox_test_result_t %>% 
        as.numeric() %>% 
        as_data_frame() %>% 
        t() %>% 
        as_tibble() %>% 
        setNames(names(wilcox_test_result_t)) 
      
      #create output table
      compare_means_result <- bind_cols(compare_means_result,
                                        wilcox_test_result_t_num)
      
      if(FUNC_OPTION_compare_means_method == "wilcox.test"){
        compare_means_result$method <- "wilcox.test"
      }
      }
    
  
    #SECOND box plots
    
    if(FUNC_OPTION_plot_qc == FALSE){
      temp_func_data <- temp_func_data %>% 
        filter(plot_class != "qc") %>%
        filter(plot_class != "QC")
    }
    
    temp_func_plot_data <- temp_func_data
    
    #browser()
    if(FUNC_OPTION_plot_outliers == FALSE){
      FUNC_Q1 <- quantile(temp_func_plot_data$concentration, .25)
      FUNC_Q3 <- quantile(temp_func_plot_data$concentration, .75)
      FUNC_IQR <- IQR(temp_func_plot_data$concentration)
      
      temp_func_plot_data <- subset(temp_func_plot_data, 
                            temp_func_plot_data$concentration > (FUNC_Q1 - 1.5*FUNC_IQR) & 
                              temp_func_plot_data$concentration < (FUNC_Q3 + 1.5*FUNC_IQR)
      )
    }
    
    if(FUNC_OPTION_log_plot_data == TRUE){
      y_axis_title <- "LOG Concentration (µM)"
      temp_func_plot_data$concentration <- log(temp_func_plot_data$concentration+1)
    }
    
    if(FUNC_OPTION_log_plot_data == FALSE){
      y_axis_title <- "Concentration (µM)"
    }
    
    bp_y_max <- temp_func_plot_data$concentration %>% max() 
    
    #start to build boxplot
      bp <- ggplot(data=temp_func_plot_data,
                   aes(x=as.factor(plot_class),
                       y=as.numeric(concentration))
      ) 

    bp <- bp + geom_boxplot(color = "black",
                      aes(),
                      outlier.shape = NA,
                      #lwd = 1,
                      alpha=0.05)  
    bp <- bp +  geom_jitter(aes(fill = get(FUNC_HEADER_temp_colour)), 
                  width = 0.01,
                  size = FUNC_OPTION_point_size,
                  shape = 21)
    bp <- bp + scale_fill_manual(values = c(FUNC_OPTION_colour_choice))
                  #size = 1)
    bp <- bp + labs(x = paste(""), y = paste(y_axis_title))
    bp <- bp +  ggtitle(idx_metabolite)
    bp <- bp +  theme_cowplot() 
    bp <- bp +  theme(plot.title = element_text(hjust = 0.5)) 
    bp <- bp +  theme(plot.title = element_text(size=5)) 
    bp <- bp +  theme(axis.text.y = element_text(size = 5, margin = margin(t = 0, r = 0, b = 0, l = 2)))
    bp <- bp +  theme(axis.text.x = element_text(size = 5, angle = 45, vjust = 1, hjust = 1))
    bp <- bp + theme(axis.title = element_text(size = 5)) 
    bp <- bp + ylim(0,NA)
    bp$labels$fill <- paste0(FUNC_HEADER_temp_colour) %>% str_to_title()
    
    
    #add significance to plot
    if(FUNC_OPTION_compare_means_method != "wilcox.test.paired"){
      
      #filter results to only annotate sig comparisons
      plot_coord <- dunn_test_result %>% filter(p.adj.signif != "ns") 
      
      if(nrow(plot_coord) > 0){
      
      #create a list of comparisons and thir y coordinates for annotating sig values on boxplots
      plot_comparisons <- list()
      y_positions <- NULL
      for(idx_comp in 1:nrow(plot_coord)){
        # list comparisons
        plot_comparisons[[idx_comp]] <- c(plot_coord$group1[idx_comp],
                                          plot_coord$group2[idx_comp])
        # concat coordinates
        y_positions <- c(y_positions, (bp_y_max + (idx_comp * bp_y_max *0.1)))
        x_min <- c(x_min, which(unique(bp$data$plot_class) == plot_coord$group1[idx_comp]))
        x_max <- c(x_max, which(unique(bp$data$plot_class) == plot_coord$group2[idx_comp]))
      }
      
      #add co-ordinate columns to tibble
      plot_coord  <- plot_coord %>% 
        add_column(y.position = y_positions) %>% 
        #add_x_position()
        add_column(xmin = x_min) %>%
        add_column(xmax = x_max)
   
    #add to bp
    bp <- bp + stat_pvalue_manual(plot_coord, label = "p.adj.signif", tip.length = 0.008, vjust = 0.5, size = 2.5)
      }
    }
    
    #add legend
     bp <- bp + theme(legend.position = "right",
           legend.title = element_text(color = "black", size = 5),
           legend.text = element_text(color = "black", size = 5))
           #legend.key.size = unit(2, "cm"))
     
     #show legend TRUE/FALSE
     if(FUNC_OPTION_plot_show_legend == FALSE){
       bp <- bp + theme(legend.position="none") 
     }
  
    # store plots for printing
    bp_plotlist$bp[[idx_metabolite]] <- bp
    
    #add mean and standard deviation to compare_means_result
    compare_means_result <- compare_means_result %>%
      bind_cols(
        temp_func_data %>% 
          group_by(plot_class) %>% 
          summarise(mean = mean(concentration)) %>% 
          spread(plot_class, mean) %>%
          rename_with( ~ paste0(.x, " mean")),
        temp_func_data %>% 
          group_by(plot_class) %>% 
          summarise(sd = sd(concentration)) %>% 
          spread(plot_class, sd) %>%
          rename_with( ~ paste0(.x, " sd"))
      )
    
    
    #View(compare_means_result)
    bp_plotlist$table[[idx_metabolite]] <- compare_means_result 
  }
  
  #combine table into single tibble
  bp_plotlist$table <- bp_plotlist$table %>% bind_rows()
  
  bp_plotlist
  #compare_means_list
}

