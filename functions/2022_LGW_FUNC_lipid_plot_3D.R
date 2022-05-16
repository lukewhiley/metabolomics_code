#ANPC PCA quality control visualisation

# FUNC_data = data containing individual lipid data - MUST CONTAIN sampleID column 
# FUNC_colour_by = how to colour the plot (e.g. sample class, is_ltr or cohort)
# FUNC_plot label = what to label the scores plot with (e.g. sampleID)
# FUNC_scaling = UV or Pareto

lgw_lipid_plot <- function(FUNC_data,
                           FUNC_plot_comparisons
                    #FUNC_colour_by, 
                    #FUNC_plot_label, 
                    #FUNC_title,
                    #FUNC_project_colours,
                    #FUNC_plot_class_or_sidechain
                    ){
  require(metabom8)
  require(RColorBrewer)
  require(tidyverse)
  require(plotly)
  
  lipid_plot_output <- list()
  
  #browser()
  
  #find max -log10 (p)
 max_log10_p <- FUNC_data %>% 
    select(all_of(FUNC_plot_comparisons)) %>%
    log() %>%
    abs() %>%
    max()
  
  
  for(idx_str_data in FUNC_plot_comparisons){
  
  #browser()
  
  #title_text <- FUNC_title
  
  # create plot values
  plot_Val <- FUNC_data %>% 
    select(feature_idx, sidechain, subclass, feature, p, all_of(FUNC_plot_comparisons)) %>%
    rename(lipid_class = subclass)
  
  #browser()
  
    temp_column_idx <- which(colnames(plot_Val) == idx_str_data)
    plot_Val_2 <- plot_Val %>%
      add_column(plot_Val[,temp_column_idx] %>%
                   log() %>%
                   abs() %>%
                   setNames(paste0("-Log10 (p)")))
 
  
 # plot_Val_2 <- plot_Val %>%
  
  #browser()
  
  #create factor for subclass
  plot_Val_2$lipid_class_factor <- plot_Val_2$lipid_class %>% factor(levels = unique(plot_Val_2$lipid_class %>% sort()), ordered = TRUE)
  
  plot_Val_2$lipid_sidechain_factor <- plot_Val_2$sidechain %>% factor(levels = unique(plot_Val_2$sidechain %>% sort()), ordered = TRUE)

  
  #create values for point size
  plot_Val_2 <- plot_Val_2 %>%
    add_column(point_size = rank(plot_Val_2$`-Log10 (p)`))
  
  #reset scale so that all features p<0.05 are set to size 1
  plot_Val_2$point_size <- plot_Val_2$point_size - which(plot_Val_2$`-Log10 (p)` < (log(0.05) %>% abs())) %>% length()
  
  #ser scale so that all features p<0.05 are set to size 1
  plot_Val_2$point_size[which(plot_Val_2$`-Log10 (p)` < (log(0.05) %>% abs()))] <- 1
 
  #arrange by point_size
  plot_Val_2 <- plot_Val_2 %>% arrange(`-Log10 (p)`)
  
  plot_Val_2$point_size_rank <- 1
  
  #split off all features that are significant
  temp_plot_Val <- plot_Val_2 %>%
    filter(point_size >1)
  
  #set breaks
  temp_plot_Val$point_size_rank <- as.integer(cut(temp_plot_Val$point_size, breaks = 100))
  #define size breaks
  temp_plot_Val$point_size <- 2.5
  temp_plot_Val$point_size[which(temp_plot_Val$point_size_rank>25)] <- 5
  temp_plot_Val$point_size[which(temp_plot_Val$point_size_rank>50)] <- 10
  temp_plot_Val$point_size[which(temp_plot_Val$point_size_rank>80)] <- 20
  temp_plot_Val$point_size[which(temp_plot_Val$point_size_rank>95)] <- 40
  
  #re-build plot_Val_2 table
  plot_Val_2 <- plot_Val_2 %>%
    filter(point_size == 1) %>%
    bind_rows(temp_plot_Val) %>% 
    arrange(feature_idx)
  
  #browser()
  
  #produce plot
  bp <- ggplot(data=plot_Val_2,
               aes(x=lipid_class_factor,
                   y=lipid_sidechain_factor)
               )
  
  bp <- bp + geom_point(aes(text = feature,
                            fill = `-Log10 (p)`,
                            color = `-Log10 (p)`,
                            size = point_size
                            ),
                          position=position_dodge(0.7),
                          width = 0.01,
                          #size = 0.5,
                          shape = 21
                         )
  bp <- bp + labs(x = paste("lipid class"),
                  y = paste("sidechain"))
  bp <- bp + ggtitle(paste0("Visualisation of Kruskal Wallis post-hoc Dunn's Tests: ", idx_str_data))
  bp <- bp + theme_cowplot() 
  bp <- bp + theme(plot.title = element_text(hjust = 0.5)) 
  bp <- bp + theme(plot.title = element_text(size=10)) 
  bp <- bp + theme(axis.text.y = element_text(size = 10, margin = margin(t = 0, r = 0, b = 0, l = 2)))
  bp <- bp + theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1))
  bp <- bp + theme(axis.title = element_text(size = 10)) 
  bp <- bp + theme(legend.title=element_text(size=10), 
                   legend.text=element_text(size=10))
  
  #create vertical lines to seprate classes on plot
  x_lipid_sequence <- seq(1:(length(plot_Val_2$lipid_class %>% unique())-1))+0.5
  y_lipid_sequence <- seq(1:(length(plot_Val_2$lipid_sidechain_factor %>% unique())-1))+0.5
  
  bp <- bp + geom_vline(xintercept=c(x_lipid_sequence),color="grey")
  bp <- bp + geom_hline(yintercept=c(y_lipid_sequence),color="grey")
  bp <- bp + guides(size=FALSE)
  
  bp <- bp + scale_size(limits = c(0,100))
  #bp <- bp + scale_fill_gradient()
  bp <- bp + scale_fill_viridis_c(option = "magma", limits = c(0, max_log10_p))
  bp <- bp + scale_color_viridis_c(option = "magma", limits = c(0, max_log10_p))
  #bp$labels$fill <- paste0(FUNC_HEADER_temp_colour) %>% str_to_title()
  
  lipid_plot_output[[idx_str_data]] <- bp %>% ggplotly() %>% layout(legend = list(orientation = 'h'))
  
  }
  lipid_plot_output
  
}
  

