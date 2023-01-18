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
                    FUNC_OPTION_log_data,
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
  
  #create PCA model
  pca_output$pca_model <- ropls::opls(x = pca_x,
                                      y = NULL,
                                      algoC = "nipals",
                                      predI = 2,
                                      crossvalI	= 5,
                                      log10L = FUNC_OPTION_log_data,
                                      scaleC = FUNC_OPTION_scaling,
                                      subset = NULL,
                                      fig.pdfC = 'none',
                                      info.txtC = 'none'
  )
  
  # extract scores values for plotting in plot_ly
  plotly_scores_data <- pca_output$pca_model@scoreMN %>% as_tibble(rownames = "sample_idx") %>%
    add_column(pca_colour = FUNC_data[[FUNC_HEADER_colour_by]],
               pca_plot_label = FUNC_data[[FUNC_HEADER_plot_label]]
    )
  
  # extract loadings values for plotting in plot_ly
  plotly_loadings_data <- pca_output$pca_model@loadingMN %>% as_tibble(rownames = "variable")
  
  #produce plot_ly PCA scores plot
  
  #axis settings
  x_axis_settings_scores <- list(
    zeroline = TRUE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    tickfont = list(size = 20),
    titlefont = list(size = 20),
    title = paste0("PC1 (", round(pca_output$pca_model@modelDF$R2X[1]*100,1), " %)")
  )
  
  y_axis_settings_scores <- list(
    zeroline = TRUE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    tickfont = list(size = 20),
    titlefont = list(size = 20),
    title = paste0("PC2 (", round(pca_output$pca_model@modelDF$R2X[2]*100,1), " %)")
  )
  
  pca_output$plot_scores <- plot_ly(type = "scatter", 
                                    mode = "markers", 
                                    data = plotly_scores_data, 
                                    x = ~p1, 
                                    y = ~p2, 
                                    text = ~pca_plot_label, 
                                    color = ~pca_colour, 
                                    colors = FUNC_OPTION_plot_colours,
                                    legendgroup = ~pca_colour,
                                    showlegend = NULL,
                                    marker = list(size = 10,
                                                  opacity = 1,
                                                  line = list(
                                                    color = '#000000',
                                                    width = 1))) %>% 
    layout(
      xaxis = x_axis_settings_scores,
      yaxis = y_axis_settings_scores,
      #showlegend = TRUE, 
      margin = list(l = 65, r = 50, b=65, t=85),
      title = paste0(FUNC_STRING_title, ": PCA Scores\n", nrow(plotly_scores_data), " samples; ", nrow(plotly_loadings_data), " features")
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
                                      x = ~p1, 
                                      y = ~p2, 
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
      title = paste0(FUNC_STRING_title,": PCA Loadings", "\n", nrow(plotly_scores_data), " samples; ", nrow(plotly_loadings_data), " features")
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
                 sprintf("<b>%s</b>", FUNC_STRING_title), "\n", nrow(plotly_scores_data), " samples; ", nrow(plotly_loadings_data), " features")
  )
  
  pca_output
  
}
