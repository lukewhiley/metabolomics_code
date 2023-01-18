#ANPC PCA quality control visualisation

# FUNC_data = data containing individual lipid data - MUST CONTAIN sampleID column 
# FUNC_colour_by = how to colour the plot (e.g. sample class, is_ltr or cohort)
# FUNC_plot label = what to label the scores plot with (e.g. sampleID)
# FUNC_scaling = UV or Pareto

lgw_pca <- function(FUNC_data, 
                    FUNC_metabolite_list, 
                    FUNC_HEADER_colour_by, 
                    FUNC_HEADER_plot_label, 
                    FUNC_STRING_title,
                    FUNC_STRING_legend_title,
                    FUNC_OPTION_scaling,
                    FUNC_OPTION_plot_colours,
                    FUNC_OPTION_invert_y,
                    FUNC_OPTION_invert_x,
                    FUNC_OPTION_include_qc,
                    FUNC_OPTION_include_sample){
  
  pca_output <- list()
  
  #remove QCs if FUNC_OPTION_include_qc == FALSE
  if(FUNC_OPTION_include_qc == FALSE){
    FUNC_data <- FUNC_data %>% 
      filter(sample_type != "qc")
  }
  
  #remove samples if FUNC_OPTION_include_sample == FALSE
  if(FUNC_OPTION_include_sample == FALSE){
    FUNC_data <- FUNC_data %>% 
      filter(sample_type != "sample")
  }
  
  
  #create data matrix for PCA
  pca_x <- FUNC_data %>%  select(all_of(FUNC_metabolite_list)) %>% as.matrix() 
  pca_x[pca_x == 0] <- NA #remove all 0 values (above adds 1 to all values therefore anything that = 1 was a 0)
  pca_x[is.infinite(pca_x)] <- NA #remove all infinite values
  min_value <- min(pca_x, na.rm = TRUE) # find the lowest value in the matrix
  pca_x[is.na(pca_x)] <- min_value/100 # replace all NA, Inf, and 0 values with the lowest value in the matrix/100 to represent a value below limit of detection
  
  pca_x <- log(pca_x+1) #log values for plotting
  
  #create PCA model
  pca_output$pca_model <- pca(pca_x, pc=3, scale = paste(FUNC_OPTION_scaling), center = TRUE)
  
  # extract score values for plotting in plot_ly
  PC1 <- as.numeric(as.matrix(pca_output$pca_model@t[,1]))
  PC2 <- as.numeric(as.matrix(pca_output$pca_model@t[,2]))
  
  # extract loadings values for plotting in plot_ly
  plotly_loadings_data <- pca_output$pca_model@p %>% as_tibble(rownames = "variable") %>% rename(PC1 = V1, PC2 = V2)
  
  #produce plot_ly PCA scores plot
  
  # set plot attributes (controlled by FUNC_colour_by and FUNC_plot_label)
  pca_colour <- FUNC_data %>% select(all_of(FUNC_HEADER_colour_by)) #%>% as.matrix()
  colnames(pca_colour) <- "pca_colour" 
  pca_colour <- pca_colour$pca_colour
  #pca_colour[is.na(pca_colour)] <- "none"

  #scores plot label
  pca_plot_label <- FUNC_data %>% 
    select(all_of(FUNC_HEADER_plot_label)) %>% 
    as.matrix() %>% c()
  
  # create plot values
  plot_Val <- as_tibble(cbind(PC1, PC2))
  plot_Val$pca_colour <- pca_colour
  plot_Val$pca_plot_label <- pca_plot_label
  
  
  #axis settings
  x_axis_settings_scores <- list(
    zeroline = TRUE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    tickfont = list(size = 20),
    titlefont = list(size = 20),
    title = paste("PC1 (", round(pca_output$pca_model@Parameters$R2[1]*100,1), " %)", sep = "")
  )
  
  y_axis_settings_scores <- list(
    zeroline = TRUE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    tickfont = list(size = 20),
    titlefont = list(size = 20),
    title = paste("PC2 (", round(pca_output$pca_model@Parameters$R2[2]*100,1), " %)", sep = "")
  )
  
  pca_output$plot_scores <- plot_ly(type = "scatter", 
                                    mode = "markers", 
                                    data = plot_Val, 
                                    x = ~PC1, 
                                    y = ~PC2, 
                                    text = ~pca_plot_label, 
                                    color = ~pca_colour, 
                                    #colors = c(plot_colors[1:length(unique(pca_colour))]), 
                                    colors = FUNC_OPTION_plot_colours,
                                    legendgroup = ~pca_colour,
                                    showlegend = NULL,
                                    marker = list(size = 10, 
                                                  #color = '#1E90FF', 
                                                  opacity = 1,
                                                  line = list(
                                                    color = '#000000',
                                                    width = 1)
                                    )) %>% 
    layout(
      xaxis = x_axis_settings_scores,
      yaxis = y_axis_settings_scores,
      #showlegend = TRUE, 
      margin = list(l = 65, r = 50, b=65, t=85),
      title = paste0(FUNC_STRING_title, ": PCA Scores\n", nrow(plot_Val), " samples; ", nrow(plotly_loadings_data), " features")
    )
  
  if(FUNC_OPTION_invert_x == TRUE){
    pca_output$plot_scores <- pca_output$plot_scores %>%
      layout(xaxis = list(autorange = "reversed"))
  }
  
  if(FUNC_OPTION_invert_y == TRUE){
    pca_output$plot_scores <- pca_output$plot_scores %>%
      layout(yaxis = list(autorange = "reversed"))
  }
  
  
  
  # create loadings plot
  x_axis_settings_loading <- list(
    zeroline = TRUE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    tickfont = list(size = 20),
    titlefont = list(size = 20),
    title = paste("")
  )
  
  y_axis_settings_loading <- list(
    zeroline = TRUE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    tickfont = list(size = 20),
    titlefont = list(size = 20),
    title = paste("")
  )
  
  pca_output$plot_loadings <- plot_ly(type = "scatter", 
                                      mode = "markers", 
                                      data = plotly_loadings_data, 
                                      x = ~PC1, 
                                      y = ~PC2, 
                                      text = ~variable, 
                                      #color = "Metabolite", 
                                      marker = list(size = 10, color = '#808080', opacity = 0.5,
                                                    line = list(color = '#000000', width = 1)
                                      )) %>% 
    layout(
      xaxis = x_axis_settings_loading,
      yaxis = y_axis_settings_loading,
      showlegend = TRUE, 
      margin = list(l = 65, r = 50, b=65, t=85),
      title = paste0(FUNC_STRING_title,": PCA Loadings", "\n", nrow(plot_Val), " samples; ", nrow(plotly_loadings_data), " features")
    )
  
  pca_output$plot_combined <- subplot(pca_output$plot_scores, 
                                      pca_output$plot_loadings, 
                                      nrows = 1,
                                      margin = 0.05,
                                      titleX = TRUE,
                                      titleY = TRUE
  ) %>% layout(legend = list(title=list(text=sprintf("<b>%s</b>", FUNC_STRING_legend_title))),
               showlegend = TRUE, 
               margin = list(l = 65, r = 50, b=65, t=85),
               title = paste0(
                 sprintf("<b>%s</b>", FUNC_STRING_title), "\n", nrow(plot_Val), " samples; ", nrow(plotly_loadings_data), " features")
  )
  
  pca_output
  
}
