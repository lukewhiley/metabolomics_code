# boxplot function

#log_data = TRUE/FALSE

lgw_generic_boxplot <- function(FUNC_data, FUNC_metabolite_list, FUNC_HEADER_class, FUNC_OPTION_log_data){
  require(metabom8)
  require(RColorBrewer)
  require(tidyverse)
  require(plotly)
  
  #browser()
  
  bp_plotlist <- list()
  
  for(idx_metabolite in FUNC_metabolite_list){
    
    #browser()
    # prepare plot_ly data
    temp_plot_data <- FUNC_data %>%
      select(all_of(FUNC_HEADER_class), all_of(idx_metabolite)) %>%
      rename(concentration = all_of(idx_metabolite), 
             plotClass = all_of(FUNC_HEADER_class))
    
    if(FUNC_OPTION_log_data == TRUE){
      temp_plot_data$concentration <- log(temp_plot_data$concentration+1)
    }
    
    
    # repare plot_ly attributes
    plot_colors <- RColorBrewer::brewer.pal(name = "Set2",
                                            n = length(unique(temp_plot_data$plotClass)))
    
    
    title_text <- paste0(idx_metabolite)
    
    x_axis_settings_scores <- list(
      zeroline = TRUE,
      showline = TRUE,
      linecolor = toRGB("black"),
      linewidth = 2,
      showgrid = TRUE,
      title = paste("cohort")
    )
    
    y_axis_settings_scores <- list(
      zeroline = TRUE,
      showline = TRUE,
      linecolor = toRGB("black"),
      linewidth = 2,
      showgrid = TRUE,
      title = paste("concentration")
    )
    
    if(FUNC_OPTION_log_data == TRUE){
    y_axis_settings_scores <- list(
      zeroline = TRUE,
      showline = TRUE,
      linecolor = toRGB("black"),
      linewidth = 2,
      showgrid = TRUE,
      title = paste("concentration (LOG)")
    )
    }
    
    # create plot_ly
    bp <- plot_ly(type = "box",
                  #mode = "markers",
                  temp_plot_data, 
                  x = ~plotClass, 
                  y = ~concentration, 
                  color = ~plotClass,
                  colors = c(plot_colors[1:length(unique(temp_plot_data$plotClass))]), 
                  marker = list(size = 7, 
                                opacity = 0.5,
                                line = list(color = '#000000',
                                            width = 1)
                  )
    ) %>% 
      layout(title = paste0(title_text),
             xaxis = x_axis_settings_scores,
             yaxis = y_axis_settings_scores)
    
   
    bp_plotlist[idx_metabolite] <- list(bp)
  }
  
  bp_plotlist
}

