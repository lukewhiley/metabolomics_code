#ANPC Lipidomics opls quality control visualisation

lipids_opls_2_class <- function(individual_multivariate_data, family_multivariate_data, multivariate_class, plot_label){
  #browser()
  lipid <- individual_multivariate_data %>% select(contains("(")) %>% colnames()
  lipid_class <- sub("\\(.*", "", lipid) %>% unique()
  lipid_class <- lipid_class[!grepl("sampleID", lipid_class)] %>% as_tibble()
  
  multivariate_data_list <- list(individual_multivariate_data, family_multivariate_data)
  
  opls_plot_list <- lapply(multivariate_data_list, function(func_list){
   #browser()
    multivariate_data <- func_list
    column_length <- multivariate_data %>% select(contains("(")) %>% ncol()
    
    if(column_length > 0){ 
      opls_x <- multivariate_data %>%  select(all_of(lipid)) %>% as.matrix()+1 
      opls_x <- log(opls_x)
      title_text <- "individual lipid species"
    }
    if(column_length == 0){ 
      opls_x <- multivariate_data %>%  select(all_of(lipid_class$value)) %>% as.matrix()+1 
      opls_x <- log(opls_x)
      title_text <- "lipid class"
    }
    
    opls_x[opls_x == 0] <- NA #remove all 0 values
    opls_x[is.infinite(opls_x)] <- NA #remove all infinite values
    min_value <- min(opls_x, na.rm = TRUE) # find the lowest value in the matrix
    opls_x[is.na(opls_x)] <- min_value # replace all NA, Inf, and 0 values with the lowest value in the matrix
    
    opls_class <- multivariate_data %>% select(multivariate_class) %>% as.matrix()
    opls_plot_label <- multivariate_data %>% select(plot_label) %>% as.matrix()
    opls_class[is.na(opls_class)] <- "none"
    sampleID <- multivariate_data %>% select(sampleID)
    
    opls_model <- opls(opls_x, opls_class, 
                       scale = "UV", 
                       center = TRUE, 
                       cv = list(method = "MC", k = 3, split = 2/3)
    )
    opls_stats <- opls_model@summary
    t_pred <- as.numeric(as.matrix(opls_model@t_pred[,1]))
    t_orth <- as.numeric(as.matrix(opls_model@t_orth[,1]))
    
    #browser()
    #produce plot_ly opls plot
    plot_Val <- as_tibble(cbind(t_pred, t_orth))
    plot_Val$sampleID <- sampleID$sampleID
    plot_Val$sample_group <- c(opls_class)
    plot_Val$opls_plot_label <- c(opls_plot_label)
    
    #plotly_loadings_data <- opls_model@p %>% as_tibble(rownames = "lipid") %>% rename(t_pred = V1, t_orth = V2)
    
    x_axis_settings_scores <- list(
      zeroline = TRUE,
      showline = TRUE,
      linecolor = toRGB("black"),
      linewidth = 2,
      showgrid = TRUE,
      title = paste("t_pred")
    )
    
    y_axis_settings_scores <- list(
      zeroline = TRUE,
      showline = TRUE,
      linecolor = toRGB("black"),
      linewidth = 2,
      showgrid = TRUE,
      title = paste("t_orth")
    )
    
    plotly_opls <- plot_ly(type = "scatter", mode = "markers", data = plot_Val, x = ~t_pred, y = ~t_orth, text =~opls_plot_label, color = ~sample_group, colors = c("#1E90FF", "#FF0000"), 
                          marker = list(size = 7, 
                                        opacity = 0.5,
                                        line = list(
                                          color = '#000000',
                                          width = 1)
                          )) %>% 
      layout(title = paste(project_name, " Plotly opls - ", title_text, sep = ""),
             xaxis = x_axis_settings_scores,
             yaxis = y_axis_settings_scores
      )
    
    #browser()
    loadings_data <- eruption(opls_model)
    
    plotly_loadings_data <- loadings_data$data
    
    
    x_axis_settings_loading <- list(
      zeroline = TRUE,
      showline = TRUE,
      linecolor = toRGB("black"),
      linewidth = 2,
      showgrid = TRUE,
      title = paste("Cliff's Delta"),
      range = c(-1,1)
    )
    
    y_axis_settings_loading <- list(
      zeroline = TRUE,
      showline = TRUE,
      linecolor = toRGB("black"),
      linewidth = 2,
      showgrid = TRUE,
      title = paste(""),
      range = c(0, max(plotly_loadings_data$p1)*1.1)
    )
    
 
    #browser()
    plotly_loadings <- plot_ly(type = "scatter", mode = "markers", data = plotly_loadings_data, x = ~Cd, y = ~p1, text = ~id, 
                               color = "OPLS loadings",
                               marker = list(size = 7, color = '#808080', opacity = 0.5,
                                             line = list(color = '#000000', width = 1)
                               )) %>% 
      layout(title = paste(project_name, " OPLS eruption plot - ", title_text, sep = ""),
             xaxis = x_axis_settings_loading,
             yaxis = y_axis_settings_loading,
             showlegend = FALSE
      )
    
    combined_plotly <- subplot(plotly_opls, plotly_loadings, 
                               margin = c(0.05, 0.05, 0.02, 0.02),
                               titleX = TRUE,
                               titleY = TRUE
    ) %>% layout(showlegend = TRUE, title =  "")
    
    
    #opls_plot_list <- c(list(combined_plotly))
    #print(opls_plot_list)
    list(combined_plotly,
         plotly_loadings_data,
      opls_stats
    )
  })
  
  #opls_plot_list <- c(opls_plot_list)
                     
  
  opls_plot_list
  
  
  
}
