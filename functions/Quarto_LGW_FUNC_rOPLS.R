#ANPC Lipidomics opls quality control visualisation

## REQUIRED PACKAGES

# -> metabom8 (github.com/tkimhofer)
# -> tidyverse
# -> RColorBrewer
# -> plotly


## REQUIRED ARGUMENTS

# -> FUNC_data = a tibble or data from containing data
# -> FUNC_opls_y = column name for column containing class data as y in OPLS-DA
# -> FUNC_metabolite_list = array of metabolites to use - must match appropiate column names
# -> FUNC_OPTION_colour_by = column name for column containing character string or factor to colour OPLS-DA 
# -> FUNC_OPTION_plot_label = column name for column containing character string or factor to label OPLS-DA plotly
# -> FUNC_OPTION_scaling = scaling argument for metabom8 - only use UV or Pareto
# -> FUNC_OPTION_title = title for OPLS-DA plot
# -> FUNC_OPTION_project_colours = array of colours - must match length of unique number of groups

# -> FUNC_data_predict = use if wanting to predict data to model - a tibble or data from containing data to predict. Set as FALSE if no predicition required


lgw_rOPLS <- function(FUNC_data, 
                     FUNC_metabolite_list,
                     FUNC_HEADER_class,
                     FUNC_OPLS_comparison_control,
                     FUNC_OPLS_comparison_test,
                     FUNC_OPTION_colour_by, 
                     FUNC_OPTION_plot_label, 
                     FUNC_OPTION_scaling,
                     FUNC_OPTION_log_data,
                     FUNC_OPTION_title,
                     FUNC_OPTION_project_colours,
                     FUNC_OPTION_max_orth,
                     FUNC_OPTION_invert_x,
                     FUNC_OPTION_invert_y
                     ){

  opls_output <- list()

for(FUNC_idx_str_opls in FUNC_OPLS_comparison_test){
    
    #browser()
    
  #create data matrix for opls
  opls_x <- FUNC_data %>%  
    filter(sample_type == "sample") %>% 
    filter(.data[[FUNC_HEADER_class]] == FUNC_OPLS_comparison_control|
             .data[[FUNC_HEADER_class]] == FUNC_idx_str_opls) %>%
    select(all_of(FUNC_metabolite_list)) %>% 
    as.matrix()
  
  #extract y values for opls_y
  opls_y <-  FUNC_data %>%  
    filter(sample_type == "sample") %>% 
    filter(.data[[FUNC_HEADER_class]] == FUNC_OPLS_comparison_control|
             .data[[FUNC_HEADER_class]] == FUNC_idx_str_opls) %>%
    select(.data[[FUNC_HEADER_class]]) %>% 
    as.matrix()
  
  #set random seed for reproducibility
  set.seed(123)
  #create opls model
  
      opls_output[[FUNC_idx_str_opls]]$opls_model <- ropls::opls(x = opls_x,
                                                                 y = opls_y,
                                                                 predI = 1,
                                                                 orthoI = 1,
                                                                 crossvalI	= 5,
                                                                 log10L = FUNC_OPTION_log_data,
                                                                 scaleC = FUNC_OPTION_scaling,
                                                                 subset = NULL,
                                                                 fig.pdfC = 'none',
                                                                 info.txtC = 'none'
                                                                 )
 
  # extract score values for plotting in plot_ly
  tPred <- as.numeric(as.matrix(opls_output[[FUNC_idx_str_opls]]$opls_model@scoreMN))
  tOrth <- as.numeric(as.matrix(opls_output[[FUNC_idx_str_opls]]$opls_model@orthoScoreMN))
  
#calculate cliffs delta for each metabolite
  cliffs_delta_out <- array(dim = ncol(opls_x))
  
  for(idx_cd in 1:ncol(opls_x)){
    control_idx <- which(opls_y == FUNC_OPLS_comparison_control)
    cliffs_delta_out[idx_cd] <- effsize::cliff.delta(opls_x[-control_idx, idx_cd], opls_x[control_idx, idx_cd])[["estimate"]]
  }
  
  #create eruption plot table for plotting
  plotly_eruption_data <- tibble(
    metabolite = opls_output[[FUNC_idx_str_opls]]$opls_model@vipVn %>% names(),
    CD = cliffs_delta_out,
    VIP = as.numeric(as.matrix(opls_output[[FUNC_idx_str_opls]]$opls_model@vipVn)),
    loading = as.numeric(as.matrix(opls_output[[FUNC_idx_str_opls]]$opls_model@loadingMN))
  )
  
  # produce plot_ly opls scores plot
    
  # set plot attributes (controlled by FUNC_OPTION_colour_by and FUNC_OPTION_plot_label)
  opls_colour <-  FUNC_data %>%  
    filter(sample_type == "sample") %>% 
    filter(.data[[FUNC_HEADER_class]] == FUNC_OPLS_comparison_control|
             .data[[FUNC_HEADER_class]] == FUNC_idx_str_opls) %>%
    select(all_of(FUNC_OPTION_colour_by)) #%>% as.matrix()
  colnames(opls_colour) <- "opls_colour" 
  opls_colour <- opls_colour$opls_colour
  #opls_colour[is.na(opls_colour)] <- "none"
  
  #set colours
  plot_colours <- FUNC_OPTION_project_colours
  
  #scores plot label
  opls_plot_label <- FUNC_data %>%  
    filter(sample_type == "sample") %>% 
    filter(.data[[FUNC_HEADER_class]] == FUNC_OPLS_comparison_control|
             .data[[FUNC_HEADER_class]] == FUNC_idx_str_opls) %>%
    select(all_of(FUNC_OPTION_plot_label)) %>% 
    as.matrix() %>% c()
  
  # create plot values
  plot_Val <- as_tibble(cbind(tPred, tOrth))
  plot_Val$opls_colour <- opls_colour
  plot_Val$opls_plot_label <- opls_plot_label
  
  #report if control is left (neg) or right (pos) of scores plot - helps contorl plotting later in function
  control_direction <- (plot_Val %>% filter(opls_colour == FUNC_OPLS_comparison_control))[["tPred"]] %>% median()
  
  #scores axis settings
  x_axis_settings_scores <- list(
    zeroline = TRUE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    tickfont = list(size = 20),
    titlefont = list(size = 20),
    title = paste("t_pred", sep = "")
  )
  
  y_axis_settings_scores <- list(
    zeroline = TRUE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    tickfont = list(size = 20),
    titlefont = list(size = 20),
    title = paste("t_orth", sep = "")
  )
  
  #prepare titles
  main_title <- paste0(FUNC_OPTION_title, "OPLS-DA: ", FUNC_OPLS_comparison_control, " vs ", FUNC_idx_str_opls)
  main_title_scores <- paste0(FUNC_OPTION_title, "OPLS-DA Scores: ", FUNC_OPLS_comparison_control, " vs ", FUNC_idx_str_opls)
  main_title_erupt <- paste0(FUNC_OPTION_title, "OPLS-DA Loadings: ", FUNC_OPLS_comparison_control, " vs ", FUNC_idx_str_opls)
  subtitle1 <- paste0(nrow(plot_Val), " samples; ", nrow(cliffs_delta_out), " features;")
  subtitle2 <- paste0("R2X = ", signif(opls_output[[FUNC_idx_str_opls]]$opls_model@summaryDF$`R2X(cum)`, 2),
                        "; Q2 = ", signif(opls_output[[FUNC_idx_str_opls]]$opls_model@summaryDF$`Q2(cum)`, 2))
  
  #create plotly 
  opls_output[[FUNC_idx_str_opls]]$plot_scores <- plot_ly(type = "scatter", 
                       mode = "markers", 
                       data = plot_Val, 
                       x = ~tPred, 
                       y = ~tOrth, 
                       text = ~opls_plot_label, 
                       color = ~opls_colour, 
                       colors = FUNC_OPTION_project_colours, 
                       legendgroup = ~opls_colour,
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
      margin = list(l = 50, r = 50, b = 80, t = 80, pad = 0, autoexpand = TRUE),
      title = list(text = sprintf("<b>%s</b> \n %s \n %s", main_title_scores, subtitle1, subtitle2),
                                 y = 0.97, x = 0.5)
    )
  
  #invert axis to keep control on lefthand side of the plot
  if(control_direction > 0){
    opls_output[[FUNC_idx_str_opls]]$plot_scores <- opls_output[[FUNC_idx_str_opls]]$plot_scores %>%
      layout(xaxis = list(autorange = "reversed"))
  }
  
  #options to invert axis
  if(FUNC_OPTION_invert_x == TRUE){
    opls_output[[FUNC_idx_str_opls]]$plot_scores <- opls_output[[FUNC_idx_str_opls]]$plot_scores %>%
      layout(xaxis = list(autorange = "reversed"))
  }
  
  if(FUNC_OPTION_invert_y == TRUE){
    opls_output[[FUNC_idx_str_opls]]$plot_scores <- opls_output[[FUNC_idx_str_opls]]$plot_scores %>%
      layout(yaxis = list(autorange = "reversed"))
  }
  
 # create loadings plot
  x_axis_settings_loading <- list(
    zeroline = TRUE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    title = paste("Cliff's Delta"),
    tickfont = list(size = 20),
    titlefont = list(size = 20),
    range = c(-1,1)
  )
  #browser()
  y_axis_settings_loading <- list(
    zeroline = TRUE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    title = paste("OPLS-DA VIP"),
    tickfont = list(size = 20),
    titlefont = list(size = 20),
    range = c(0,plotly_eruption_data$VIP %>% max()*1.20)
  )
  
  opls_output[[FUNC_idx_str_opls]]$plot_loadings <- plot_ly(type = "scatter", 
                                       mode = "markers", 
                                       data = plotly_eruption_data, 
                                       x = ~CD, 
                                       y = ~VIP, 
                                       text = ~metabolite,
                                       color = "Metabolite",
                                       marker = list(size = 10, color = '#808080', opacity = 0.5,
                                                     line = list(color = '#000000', width = 1)
                                       )) %>% 
    layout(
      xaxis = x_axis_settings_loading,
      yaxis = y_axis_settings_loading,
      showlegend = TRUE, 
      margin = list(l = 50, r = 50, b = 80, t = 80, pad = 0, autoexpand = TRUE),
      title = list(text = sprintf("<b>%s</b> \n %s \n %s", main_title_erupt, subtitle1, subtitle2),
                   y = 0.97, x = 0.5)
    )

  
  # #invert axis to keep control on lefthand side of the plot
  # if(control_direction > 0){
  #   opls_output[[FUNC_idx_str_opls]]$plot_loadings <- opls_output[[FUNC_idx_str_opls]]$plot_loadings %>%
  #     layout(xaxis = list(autorange = "reversed"))
  # }
  

  
  # create loadings barchart for visualisation
  opls_output[[FUNC_idx_str_opls]]$data_eruption <- plotly_eruption_data %>% arrange(desc(VIP))
  
  #create list of top 20 important features
  plotly_eruption_data_crop <- plotly_eruption_data %>%
    arrange(loading) %>%
    {bind_rows(slice_head(., n=10), slice_tail(., n=10))} 
  
  #reverse order of barchart if control plots on OPLS-DA right (pos of plot)
  if(control_direction > 0){plotly_eruption_data_crop <- plotly_eruption_data_crop %>%
    arrange(desc(loading))}
  
  #add column for plot colour
  plotly_eruption_data_crop[["colour"]] <- "pos"
  plotly_eruption_data_crop[["colour"]][1:10] <- "neg"
  
  #set metabolite names as factor for control of plotting
  plotly_eruption_data_crop$metabolite <- plotly_eruption_data_crop$metabolite %>%
    factor(levels = plotly_eruption_data_crop$metabolite, ordered = TRUE)
  
  #create bar plot
  bp <- ggplot(data=plotly_eruption_data_crop, 
               aes(x=metabolite, y=loading, fill = colour)) +
    geom_bar(stat="identity", colour = "black")
  bp <- bp + labs(x = paste(""), y = paste("OPLS-DA\nPredictive loading"))
  bp <- bp +  ggtitle("")
  bp <- bp +  theme_cowplot() 
  bp <- bp +  theme(plot.title = element_text(hjust = 0.5)) 
  bp <- bp +  theme(plot.title = element_text(size=12)) 
  bp <- bp +  theme(axis.text.y = element_text(size = 10, margin = margin(t = 0, r = 0, b = 0, l = 2)))
  bp <- bp +  theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1))
  bp <- bp + theme(axis.title = element_text(size = 12)) 
  bp <- bp + scale_fill_manual(values = FUNC_OPTION_project_colours[
    c(which(levels(FUNC_data[[FUNC_OPTION_colour_by]]) == FUNC_OPLS_comparison_control),
      which(levels(FUNC_data[[FUNC_OPTION_colour_by]]) == FUNC_idx_str_opls))])
  bp <- bp + scale_color_manual(values = c("black"))
  bp <- bp + theme(legend.position = "none") 
  
  if(control_direction > 0){bp <- bp + scale_y_reverse()}
  
  opls_output[[FUNC_idx_str_opls]]$loading_bar <- bp %>% plotly::ggplotly()
  
  
  #make a combine plot for export printing
  opls_output[[FUNC_idx_str_opls]]$plot_combined <- subplot(
    opls_output[[FUNC_idx_str_opls]]$plot_scores, 
    opls_output[[FUNC_idx_str_opls]]$plot_loadings,
    nrows = 1,
    margin = 0.05,
    titleX = TRUE,
    titleY = TRUE
  ) %>% 
    layout(
      showlegend = TRUE, 
      margin = list(l = 65, r = 50, b=65, t=100),
      title = list(text = sprintf("<b>%s</b> \n %s %s", main_title, subtitle1, subtitle2),
                   y = 0.95, x = 0.5,
                   font = list(size = 20)),
      legend = list(orientation = "h",
                    xanchor = "center",
                    x = 0.5, y=1.2,
                    font = list(size = 20),
                    itemsizing = "constant",
                    traceorder = "grouped")
    )

  }
  
  opls_output
  
}
