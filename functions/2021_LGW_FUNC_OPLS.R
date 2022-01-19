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
# -> FUNC_colour_by = column name for column containing character string or factor to colour OPLS-DA 
# -> FUNC_plot_label = column name for column containing character string or factor to label OPLS-DA plotly
# -> FUNC_scaling = scaling argument for metabom8 - only use UV or Pareto
# -> FUNC_title = title for OPLS-DA plot
# -> FUNC_project_colours = array of colours - must match length of unique number of groups

# -> FUNC_data_predict = use if wanting to predict data to model - a tibble or data from containing data to predict. Set as FALSE if no predicition required


lgw_opls <- function(FUNC_data, 
                     FUNC_opls_y, 
                     FUNC_metabolite_list, 
                     FUNC_colour_by, 
                     FUNC_plot_label, 
                     FUNC_scaling,
                     FUNC_title,
                     FUNC_project_colours,
                     FUNC_max_orth
                     ){
  require(metabom8)
  require(RColorBrewer)
  require(tidyverse)
  require(plotly)
  
  opls_output <- list()
  
  #browser()
  
  #create data matrix for opls
  opls_x <- FUNC_data %>%  select(all_of(FUNC_metabolite_list)) %>% as.matrix()+1 
  opls_x <- log(opls_x)
  opls_x[opls_x == 0] <- NA #remove all 0 values
  opls_x[is.infinite(opls_x)] <- NA #remove all infinite values
  min_value <- min(opls_x, na.rm = TRUE) # find the lowest value in the matrix
  opls_x[is.na(opls_x)] <- min_value # replace all NA, Inf, and 0 values with the lowest value in the matrix
  
  opls_y = FUNC_opls_y
  
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
  
  # extract loadings values for plotting in plot_ly
 
  capture.output(
    suppressMessages(  
      opls_output$eruption_model <- 
        eruption(mod = opls_output$opls_model, 
               pc = 1,
               p_adj = "BH")
    )
  )
    
    plotly_loadings_data <- opls_output$eruption_model$data %>% as_tibble() 
  
  #produce plot_ly opls scores plot
    
  # set plot attributes (controlled by FUNC_colour_by and FUNC_plot_label)
  opls_colour <- FUNC_data %>% select(all_of(FUNC_colour_by)) #%>% as.matrix()
  colnames(opls_colour) <- "opls_colour" 
  opls_colour <- opls_colour$opls_colour
  opls_colour[is.na(opls_colour)] <- "none"
  
  #set colours
  plot_colours <- FUNC_project_colours
  
  #scores plot label
  opls_plot_label <- FUNC_data %>% 
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
    title = paste("t_pred", sep = "")
  )
  
  y_axis_settings_scores <- list(
    zeroline = TRUE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    title = paste("t_orth", sep = "")
  )
  
  #create plotly 
  opls_output$plot_scores <- plot_ly(type = "scatter", 
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

  
  opls_output$plot_combined <- subplot(
    opls_output$plot_scores, 
    opls_output$plot_loadings, 
    nrows = 1,
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
  
  opls_output$data_eruption <- plotly_loadings_data %>% arrange(desc(p1))
  
  opls_output
  
}
