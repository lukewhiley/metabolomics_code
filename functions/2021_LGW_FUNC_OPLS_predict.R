#ANPC Lipidomics opls quality control visualisation

## REQUIRED PACKAGES

# -> metabom8 (github.com/tkimhofer)
# -> tidyverse
# -> RColorBrewer
# -> plotly


## REQUIRED ARGUMENTS

# -> FUNC_data = a tibble or data from containing data. Must contain headers for sample_class, sample_type
# -> FUNC_opls_class = an array of 2 string descriptors of the classes for OPLS-DA (e.g. "healthy", "control")
# -> FUNC_opls_y = column name for column containing class data as y in OPLS-DA
# -> FUNC_metabolite_list = array of metabolites to use - must match appropriate column names
# -> FUNC_colour_by = column name for column containing character string or factor to colour OPLS-DA 
# -> FUNC_plot_label = column name for column containing character string or factor to label OPLS-DA plotly
# -> FUNC_scaling = scaling argument for metabom8 - only use UV or Pareto
# -> FUNC_title = title for OPLS-DA plot
# -> FUNC_project_colours = array of colours - must match length of unique number of groups

# -> FUNC_opls_class_predict = use if wanting to predict data to model - a descriptor string of class to predict. Set as FALSE if no predicition required


lgw_opls <- function(
  # FUNC_data, 
  #                    FUNC_opls_class,
  #                    FUNC_opls_class_predict,
  #                    #FUNC_opls_y, 
  #                    FUNC_metabolite_list, 
  #                    FUNC_colour_by, 
  #                    FUNC_plot_label, 
  #                    FUNC_scaling,
  #                    FUNC_title,
  #                    FUNC_project_colours,
  #                    FUNC_max_orth
  
  FUNC_data, 
  FUNC_OPLS_comparison_control,
  FUNC_OPLS_comparison_test,
  FUNC_OPLS_predict,
  FUNC_metabolite_list,
  FUNC_colour_by,
  FUNC_plot_label,
  FUNC_scaling,
  FUNC_title,
  FUNC_project_colours = master_list$project_details$project_colours$colour_selection,
  FUNC_max_orth = 1,
  FUNC_option_invert_x = FALSE,
  FUNC_option_invert_y = FALSE
                     ){
  require(metabom8)
  require(RColorBrewer)
  require(tidyverse)
  require(plotly)
  
  opls_output <- list()
  
  TEMP_FUNC_opls_data <- FUNC_data %>%
    filter(sample_type == "sample") %>%
    filter(sample_class == FUNC_opls_class[1]| sample_class == FUNC_opls_class[2])
  
  #browser()
  
  #create data matrix for opls
  opls_x <- TEMP_FUNC_opls_data %>%
    select(all_of(FUNC_metabolite_list)) %>% 
    as.matrix()+1 
  
  
  opls_x <- log(opls_x)
  opls_x[opls_x == 0] <- NA #remove all 0 values
  opls_x[is.infinite(opls_x)] <- NA #remove all infinite values
  min_value <- min(opls_x, na.rm = TRUE) # find the lowest value in the matrix
  opls_x[is.na(opls_x)] <- min_value # replace all NA, Inf, and 0 values with the lowest value in the matrix
  
  opls_y <- TEMP_FUNC_opls_data %>%
    select("sample_class") %>% 
    as.matrix()
  
  #create opls model
  capture.output(
    suppressMessages(
      opls_output$opls_model <- 
    opls(X = opls_x, 
         Y = opls_y,
         scale = paste(FUNC_scaling),
         center = TRUE,
         maxPCo = (FUNC_max_orth +1))
  ))
  
  # extract score values for plotting in plot_ly
  PC1 <- as.numeric(as.matrix(opls_output$opls_model@t_pred))
  PC2 <- as.numeric(as.matrix(opls_output$opls_model@t_orth))
  
  # extract eruption loadings values for plotting in plot_ly
 
  capture.output(
    suppressMessages(  
      opls_output$eruption_model <- 
        eruption(mod = opls_output$opls_model, 
               pc = 1,
               p_adj = "BH")
    )
  )
    
    plotly_loadings_data <- opls_output$eruption_model$data %>% as_tibble()
    
    
  # perform OPLS-DA prediction and extract plotly loadings
    if(is.character(FUNC_opls_class_predict)){
      
      TEMP_FUNC_opls_predict_data <- FUNC_data %>%
        filter(sample_type == "sample") %>%
        filter(sample_class == FUNC_opls_class_predict)
      
    opls_x_predict <- TEMP_FUNC_opls_predict_data %>%
      select(all_of(FUNC_metabolite_list)) %>% 
      as.matrix()+1 
    
    opls_x_predict <- log(opls_x_predict)
    opls_x_predict[opls_x_predict == 0] <- NA #remove all 0 values
    opls_x_predict[is.infinite(opls_x_predict)] <- NA #remove all infinite values
    min_value <- min(opls_x, na.rm = TRUE) # find the lowest value in the origional opls_x matrix
    opls_x_predict[is.na(opls_x_predict)] <- min_value # replace all NA, Inf, and 0 values with the lowest value in the matrix
    
    opls_model_predict <- metabom8::predict_opls(opls_model =  opls_output$opls_model,
                           newdata = opls_x_predict
                          )
    
    PC1_predict <- as.numeric(as.matrix(opls_model_predict$t_pred))
    PC2_predict <- as.numeric(as.matrix(opls_model_predict$t_orth))
    }
  
    
    #produce plot_ly opls scores plot
    
  # set plot attributes (controlled by FUNC_colour_by and FUNC_plot_label)
  opls_colour <- TEMP_FUNC_opls_data %>%
    select(all_of(FUNC_colour_by)) #%>% as.matrix()
  colnames(opls_colour) <- "opls_colour" 
  opls_colour <- opls_colour$opls_colour
  opls_colour[is.na(opls_colour)] <- "none"
  
  #set colours
  plot_colours <- FUNC_project_colours
  
  #scores plot label
  opls_plot_label <- TEMP_FUNC_opls_data %>% 
    select(all_of(FUNC_plot_label)) %>% 
    as.matrix()
  
  # create plot values
  plot_Val <- as_tibble(cbind(PC1, PC2))
  plot_Val$opls_colour <- c(opls_colour)
  plot_Val$opls_plot_label <- c(opls_plot_label)
  
  #scores axis settings
  x_axis_settings_scores <- list(
    zeroline = TRUE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    title = list(text = paste("t_pred", sep = ""),
                 standoff = 5)
    #range=c(-10,10)
  )
  
  y_axis_settings_scores <- list(
    zeroline = TRUE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    title = list(text = paste("t_orth", sep = ""),
                 standoff = 5)
  )
  
  #create plotly 
  
  opls_output$plot_scores_training <- plot_ly(type = "scatter", 
                       mode = "markers", 
                       data = plot_Val, 
                       x = ~PC1, 
                       y = ~PC2, 
                       text = ~opls_plot_label, 
                       color = ~opls_colour, 
                       colors = FUNC_project_colours, 
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
      margin = list(l = 65, r = 50, b=65, t=85),
      title = paste0(FUNC_title, "\n", nrow(plot_Val), " samples; ", nrow(plotly_loadings_data), " features; R2X = ", 
                     signif(opls_output$opls_model@summary$R2X[1],2)
      )
    )
  # 
 
 
 # create loadings plot
  x_axis_settings_loading <- list(
    zeroline = TRUE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    title = paste("")
    # range(c(-10,10))
  )
  
  y_axis_settings_loading <- list(
    zeroline = TRUE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    title = paste("")
  )
  
  opls_output$plot_loadings <- plot_ly(type = "scatter", 
                                       mode = "markers", 
                                       data = plotly_loadings_data, 
                                       x = ~Cd, 
                                       y = ~p1, 
                                       text = ~id,
                                       color = "Metabolite",
                                       marker = list(size = 10, color = '#808080', opacity = 0.5,
                                                     line = list(color = '#000000', width = 1)
                                       )) %>% 
    layout(
      xaxis = x_axis_settings_loading,
      yaxis = y_axis_settings_loading,
      showlegend = TRUE, 
      margin = list(l = 65, r = 50, b=65, t=85),
      title = paste0(FUNC_title, "\n", nrow(plot_Val), " samples; ", nrow(plotly_loadings_data), " features; R2X = ", 
                     signif(opls_output$opls_model@summary$R2X[1],2)
      )
    )
  
  ##################
  ##################
  
  # #create plot values for prediction data
  # 

  
  if(is.character(FUNC_opls_class_predict)){
    opls_colour_predict <- TEMP_FUNC_opls_predict_data %>%
      select(all_of(FUNC_colour_by)) #%>% as.matrix()
    colnames(opls_colour_predict) <- "opls_colour"
    opls_colour_predict <- opls_colour_predict$opls_colour
    opls_colour_predict[is.na(opls_colour_predict)] <- "none"
    #
    # #set colours
    # plot_colours <- FUNC_project_colours

    #scores plot label
    opls_plot_label_predict <- TEMP_FUNC_opls_predict_data %>%
      select(all_of(FUNC_plot_label)) %>%
      as.matrix()

    # create plot values
    plot_Val_predict <- as_tibble(cbind(PC1_predict, PC2_predict))
    plot_Val_predict$opls_colour <- c(opls_colour_predict)
    plot_Val_predict$opls_plot_label <- c(opls_plot_label_predict)
    colnames(plot_Val_predict) <- colnames(plot_Val)

    plot_Val_predict <- rbind(plot_Val,
                      plot_Val_predict)
    
    
    
    opls_output$plot_scores_predict <- plot_ly(type = "scatter", 
                                                mode = "markers", 
                                                data = plot_Val_predict, 
                                                x = ~PC1, 
                                                y = ~PC2, 
                                                text = ~opls_plot_label, 
                                                color = ~opls_colour, 
                                                colors = FUNC_project_colours, 
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
        margin = list(l = 65, r = 50, b=65, t=85),
        title = paste0(FUNC_title, "\n", nrow(plot_Val), " samples; ", nrow(plotly_loadings_data), " features; R2X = ", 
                       signif(opls_output$opls_model@summary$R2X[1],2)
        )
      )
    
    

  }
  
  ##################
  ##################
  
  if(is.logical(FUNC_opls_class_predict)){
    opls_output$plot_combined <- subplot(
      opls_output$plot_scores_training,
      opls_output$plot_loadings,
      nrows = 2,
      margin = 0.05,
      titleX = TRUE,
      titleY = TRUE
    ) %>%
      layout(
        showlegend = TRUE,
        margin = list(l = 65, r = 50, b=65, t=85),
        title = paste0(FUNC_title, "\n", nrow(plot_Val), " samples; ", nrow(plotly_loadings_data), " features; R2X = ",
                       signif(opls_output$opls_model@summary$R2X[1],2)
                       )
      )
    
    #set widget height
    opls_output$plot_combined$height <- 500
  
  }
  
  
  if(is.character(FUNC_opls_class_predict)){
    
    #remove axis legend for all but 1 scores plot, to avoid repeat legend display
    
    opls_output$plot_scores_training$x$attrs[[1]]$showlegend <- FALSE
    
    #set common axis for all subplots between taining and predict data (10% buffer around min and max values in training and predict data)
    
    #find xaxis min and max
    axis_x_max <- max(
      max(c(plot_Val$PC1, plot_Val_predict$PC1)),
      abs(min(min(c(plot_Val$PC1, plot_Val_predict$PC1))))
      )
    axis_x_range <- c(-(axis_x_max * 1.1), axis_x_max * 1.1) %>% 
      signif(digits = 2)
    
    axis_y_max <- max(
      max(c(plot_Val$PC2, plot_Val_predict$PC2)),
      abs(min(min(c(plot_Val$PC2, plot_Val_predict$PC2))))
    )
    axis_y_range <- c(-(axis_y_max * 1.1), axis_y_max * 1.1) %>% 
      signif(digits = 2)
      
   
    
    #set axis on training and prediction plots 
    #training
    opls_output$plot_scores_training$x$layoutAttrs[[1]]$xaxis$range <- axis_x_range
    opls_output$plot_scores_training$x$layoutAttrs[[1]]$yaxis$range <- axis_y_range
    
    #prediction
    opls_output$plot_scores_predict$x$layoutAttrs[[1]]$xaxis$range <- axis_x_range
    opls_output$plot_scores_predict$x$layoutAttrs[[1]]$yaxis$range <- axis_y_range
    
    #eruption
    opls_output$plot_loadings$x$layoutAttrs[[1]]$xaxis$range <- c(-1,1)
    
    #create subplot
  opls_output$plot_combined <- subplot(
    opls_output$plot_scores_training,
    opls_output$plot_scores_predict,
    opls_output$plot_loadings,
    nrows = length(FUNC_opls_class_predict)+2,
    margin = 0.05,
    titleX = TRUE,
    titleY = TRUE,
    heights = c(0.3,0.35,0.3)
  ) %>%
    layout(
      showlegend = TRUE,
     # margin = list(l = 65, r = 50, b=65, t=85),
      title = paste0(FUNC_title, "\n", nrow(plot_Val), " samples; ", nrow(plotly_loadings_data), " features; R2X = ",
                     signif(opls_output$opls_model@summary$R2X[1],2)
                     )
    )

  #set overall widget output height
  opls_output$plot_combined$height <- 250 * (length(FUNC_opls_class_predict)+2)
  }
  
  opls_output$data_eruption <- plotly_loadings_data %>% arrange(desc(p1))
  
  opls_output
  
}




## backup code


# opls_output$plot_combined <- subplot(
#   opls_output$plot_scores, 
#   opls_output$plot_loadings, 
#   nrows = 1,
#   margin = 0.05,
#   titleX = TRUE,
#   titleY = TRUE
# ) %>% 
#   layout(
#     showlegend = TRUE, 
#     margin = list(l = 65, r = 50, b=65, t=85),
#     title = paste0(FUNC_title, "\n", nrow(plot_Val), " samples; ", nrow(plotly_loadings_data), " features; R2X = ", 
#                    signif(opls_output$opls_model@summary$R2X[1],2)
#                    )
#   )













