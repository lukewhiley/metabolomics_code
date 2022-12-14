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


lgw_opls_predict <- function(FUNC_data, 
                             FUNC_metabolite_list,
                             FUNC_HEADER_class,
                             FUNC_OPLS_comparison_control,
                             FUNC_OPLS_comparison_test,
                             FUNC_OPLS_predict,
                             FUNC_OPTION_colour_by, 
                             FUNC_OPTION_plot_label, 
                             FUNC_OPTION_scaling,
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
    
    #impute missing values
    opls_x[opls_x == 0] <- NA #remove all 0 values (above adds 1 to all values therefore anything that = 1 was a 0)
    opls_x[is.infinite(opls_x)] <- NA #remove all infinite values
    min_value <- min(opls_x, na.rm = TRUE) # find the lowest value in the matrix
    opls_x[is.na(opls_x)] <- min_value/100 # replace all NA, Inf, and 0 values with the lowest value in the matrix/100 to represent a value below limit of detection
    
    #log data
    opls_x <- log(opls_x+1) #log values for model
    
    opls_y <-  FUNC_data %>%  
      filter(sample_type == "sample") %>% 
      filter(.data[[FUNC_HEADER_class]] == FUNC_OPLS_comparison_control|
               .data[[FUNC_HEADER_class]] == FUNC_idx_str_opls) %>%
      select(.data[[FUNC_HEADER_class]]) %>% 
      as.matrix()
    
    #set random seed for reproducibility
    set.seed(123)
    #create opls model
    capture.output(
      suppressMessages(
        opls_output[[FUNC_idx_str_opls]]$opls_model <- 
          metabom8:: opls(X = opls_x, 
                          Y = opls_y,
                          scale = paste(FUNC_OPTION_scaling),
                          center = TRUE,
                          maxPCo = (FUNC_OPTION_max_orth +1))
      ))
    
    # extract score values for plotting in plot_ly
    PC1 <- as.numeric(as.matrix(opls_output[[FUNC_idx_str_opls]]$opls_model@t_pred))
    PC2 <- as.numeric(as.matrix(opls_output[[FUNC_idx_str_opls]]$opls_model@t_orth))
    
    # extract loadings values for plotting in plot_ly
    
    capture.output(
      suppressMessages(  
        opls_output[[FUNC_idx_str_opls]]$eruption_model <- 
          eruption(mod = opls_output[[FUNC_idx_str_opls]]$opls_model, 
                   pc = 1,
                   p_adj = "BH")
      )
    )
    
    plotly_loadings_data <- opls_output[[FUNC_idx_str_opls]]$eruption_model$data %>% as_tibble() 
  
  
    
  ######OPLS-DA PREDICTION OF TEST DATA######  
  # perform OPLS-DA prediction and extract plotly loadings
    
    
    opls_x_predict <- FUNC_data %>%  
      filter(sample_type == "sample") %>% 
      filter(.data[[FUNC_HEADER_class]] == FUNC_OPLS_predict) %>%
      select(all_of(FUNC_metabolite_list)) %>% 
      as.matrix()

    #impute missing values
    opls_x_predict[opls_x_predict == 0] <- NA #remove all 0 values (above adds 1 to all values therefore anything that = 1 was a 0)
    opls_x_predict[is.infinite(opls_x_predict)] <- NA #remove all infinite values
    min_value <- min(opls_x_predict, na.rm = TRUE) # find the lowest value in the matrix
    opls_x_predict[is.na(opls_x_predict)] <- min_value/100 # replace all NA, Inf, and 0 values with the lowest value in the matrix/100 to represent a value below limit of detection
    
    #log data
    opls_x_predict <- log(opls_x_predict+1) #log values for model
    
    ###### PREPARE DATA FOR TRAINING PLOTLY PLOTS ######   
    
    #PREPARE DATA
    #TRAINING
    # set plot attributes (controlled by FUNC_OPTION_colour_by and FUNC_OPTION_plot_label)
    opls_colour <- FUNC_data %>%  
      filter(sample_type == "sample") %>% 
      filter(.data[[FUNC_HEADER_class]] == FUNC_OPLS_comparison_control|
               .data[[FUNC_HEADER_class]] == FUNC_idx_str_opls) %>%
      select(all_of(FUNC_OPTION_colour_by))
    
    colnames(opls_colour) <- "opls_colour"
    opls_colour <- opls_colour$opls_colour
    
    #set colours
    plot_colours <- FUNC_OPTION_project_colours
    
    #scores plot label
    opls_plot_label <- FUNC_data %>%  
      filter(sample_type == "sample") %>% 
      filter(.data[[FUNC_HEADER_class]] == FUNC_OPLS_comparison_control|
               .data[[FUNC_HEADER_class]] == FUNC_idx_str_opls) %>%
      select(all_of(FUNC_OPTION_plot_label)) %>% 
      as.matrix() %>% c()
    
    # create plot values for training data
    plot_Val_training <- as_tibble(cbind(PC1, PC2))
    plot_Val_training$opls_colour <- opls_colour
    plot_Val_training$opls_plot_label <- opls_plot_label
  
    
    #PREDICT
    
    #make prediction model
    opls_model_predict <- metabom8::predict_opls(opls_model =  opls_output[[FUNC_OPLS_comparison_test]]$opls_model,
                                                 newdata = opls_x_predict
    )
    
    
    #extract prediction scores data for plotting in plotly
    PC1_predict <- as.numeric(as.matrix(opls_model_predict$t_pred))
    PC2_predict <- as.numeric(as.matrix(opls_model_predict$t_orth))
    
    
    ###### PREPARE DATA FOR PREDICTION PLOTLY PLOTS ######  
    
    # set plot attributes (controlled by FUNC_OPTION_colour_by and FUNC_OPTION_plot_label)
    opls_colour_predict <- FUNC_data %>%  
      filter(sample_type == "sample") %>% 
      filter(.data[[FUNC_HEADER_class]] == FUNC_OPLS_predict) %>%
      select(all_of(FUNC_OPTION_colour_by)) 
    
    colnames(opls_colour_predict) <- "opls_colour" 
    opls_colour_predict <- opls_colour_predict$opls_colour
  
    
    #scores plot label prediction data
    opls_plot_label_predict <- FUNC_data %>%  
      filter(sample_type == "sample") %>% 
      filter(.data[[FUNC_HEADER_class]] == FUNC_OPLS_predict) %>%
      select(all_of(FUNC_OPTION_plot_label)) %>% 
      as.matrix() %>% c()
     
    # create plot values for predict data
    plot_Val_predict <- as_tibble(cbind(PC1_predict, PC2_predict)) %>%
      rename(PC1 = PC1_predict, PC2 = PC2_predict)
    plot_Val_predict$opls_colour <- opls_colour_predict
    plot_Val_predict$opls_plot_label <- opls_plot_label_predict
    

    
    #COMBINED PLOT_DATA
    plot_Val <- bind_rows(plot_Val_training, plot_Val_predict)
    
    
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
    
    #create plotly 
    opls_output[[FUNC_idx_str_opls]]$plot_scores <- plot_ly(type = "scatter", 
                                                            mode = "markers", 
                                                            data = plot_Val, 
                                                            x = ~PC1, 
                                                            y = ~PC2, 
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
        margin = list(l = 65, r = 50, b=65, t=85),
        title = paste0(FUNC_OPTION_title, ": OPLS-DA Scores: ", FUNC_OPLS_comparison_control, " vs ", FUNC_idx_str_opls, "\n", nrow(plot_Val), " samples; ", nrow(plotly_loadings_data), " features; R2X = ", 
                       signif(opls_output[[FUNC_idx_str_opls]]$opls_model@summary$R2X[1],2)
        )
      )
    
    if(FUNC_OPTION_invert_x == TRUE){
      pca_output$plot_scores <- pca_output$plot_scores %>%
        layout(xaxis = list(autorange = "reversed"))
    }
    
    if(FUNC_OPTION_invert_y == TRUE){
      pca_output$plot_scores <- pca_output$plot_scores %>%
        layout(yaxis = list(autorange = "reversed"))
    }
    
    # 
    
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
      title = paste("|Ppred|"),
      tickfont = list(size = 20),
      titlefont = list(size = 20),
      range = c(0,plotly_loadings_data$p1 %>% max()*1.20)
    )
    
    opls_output[[FUNC_idx_str_opls]]$plot_loadings <- plot_ly(type = "scatter", 
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
        title = paste0(FUNC_OPTION_title, ": OPLS-DA Loadings: ", FUNC_OPLS_comparison_control, " vs ", FUNC_idx_str_opls, "\n", nrow(plot_Val), " samples; ", nrow(plotly_loadings_data), " features; R2X = ", 
                       signif(opls_output[[FUNC_idx_str_opls]]$opls_model@summary$R2X[1],2)
        )
      )
    
    
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
        margin = list(l = 65, r = 50, b=65, t=85),
        title = paste0(master_list$project_details$project_name, ": OPLS-DA: ", FUNC_OPLS_comparison_control, " vs ", FUNC_idx_str_opls,  "\n", nrow(plot_Val), " samples; ", nrow(plotly_loadings_data), " features; R2X = ", 
                       signif(opls_output[[FUNC_idx_str_opls]]$opls_model@summary$R2X[1],2)
        )
      )
    
    opls_output[[FUNC_idx_str_opls]]$data_eruption <- plotly_loadings_data %>% arrange(desc(p1))
    
  }
  
  opls_output
  
}

    
    