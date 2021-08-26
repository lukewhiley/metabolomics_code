#ANPC Lipidomics opls quality control visualisation

## REQUIRED PACKAGES

# -> metabom8 (github.com/tkimhofer)
# -> tidyverse
# -> RColorBrewer
# -> plotly


## REQUIRED ARGUMENTS

# -> FUNC_individual_multivariate_data = a tibble or data fram containng data
# -> FUNC_opls_y = column name for column containing class data as y in OPLS-DA
# -> FUNC_metabolite_list = array of metabolites to use - must match appropiate column names
# -> FUNC_colour_by = column name for column containing ccharacter string or factor to colour OPLS-DA 
# -> FUNC_plot_label = column name for column containing ccharacter string or factor to label OPLS-DA plotly
# -> FUNC_scaling = scaling argument for metabom8 - only use UV or Pareto


lgw_opls <- function(FUNC_individual_multivariate_data, FUNC_opls_y, FUNC_metabolite_list, FUNC_colour_by, FUNC_plot_label, FUNC_scaling){
  require(metabom8)
  require(RColorBrewer)
  require(tidyverse)
  require(plotly)
  
  #browser()
  
  #create data matrix for opls
  opls_x <- FUNC_individual_multivariate_data %>%  select(all_of(FUNC_metabolite_list)) %>% as.matrix()+1 
  opls_x <- log(opls_x)
  title_text <- "opls"
  opls_x[opls_x == 0] <- NA #remove all 0 values
  opls_x[is.infinite(opls_x)] <- NA #remove all infinite values
  min_value <- min(opls_x, na.rm = TRUE) # find the lowest value in the matrix
  opls_x[is.na(opls_x)] <- min_value # replace all NA, Inf, and 0 values with the lowest value in the matrix
  
  #create opls model
  opls_model <- opls(X = opls_x, 
                     Y = FUNC_opls_y,
                     scale = paste(FUNC_scaling), 
                     center = TRUE)
  
  # extract score values for plotting in plot_ly
  PC1 <- as.numeric(as.matrix(opls_model@t_pred))
  PC2 <- as.numeric(as.matrix(opls_model@t_orth))
  
  # extract loadings values for plotting in plot_ly
 
    eruption_model <- eruption(mod = opls_model, 
                             pc = 1,
                             p_adj = "BH")
    
    plotly_loadings_data <- eruption_model$data %>% as_tibble() 
  
  #produce plot_ly opls scores plot
  
  # set plot attributes (controlled by FUNC_colour_by and FUNC_plot_label)
  opls_colour <- FUNC_individual_multivariate_data %>% select(all_of(FUNC_colour_by)) %>% as.matrix()
  opls_colour[is.na(opls_colour)] <- "none"
  opls_plot_label <- FUNC_individual_multivariate_data %>% 
    select(all_of(FUNC_plot_label)) %>% 
    as.matrix()
  
  # create plot values
  plot_Val <- as_tibble(cbind(PC1, PC2))
  plot_Val$opls_colour <- c(opls_colour)
  plot_Val$opls_plot_label <- c(opls_plot_label)
  
  plot_colors <- RColorBrewer::brewer.pal(name = "Set2",
                                          n = length(unique(opls_colour)))
  
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
  
 plotly_opls <- plot_ly(type = "scatter", 
                       mode = "markers", 
                       data = plot_Val, 
                       x = ~PC1, 
                       y = ~PC2, 
                       text = ~opls_plot_label, 
                       color = ~opls_colour, 
                       colors = c(plot_colors[1:length(unique(opls_colour))]), 
                        marker = list(size = 10, 
                                      #color = '#1E90FF', 
                                      opacity = 0.5,
                                      line = list(
                                        color = '#000000',
                                        width = 1)
                        )) %>% 
    layout(title = paste(" Plotly opls - ", title_text, sep = ""),
           xaxis = x_axis_settings_scores,
           yaxis = y_axis_settings_scores)
  
 
 
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
  
  plotly_loadings <- plot_ly(type = "scatter", 
                             mode = "markers", 
                             data = plotly_loadings_data, 
                             x = ~Cd, 
                             y = ~p1, 
                             text = ~id, 
                             marker = list(size = 10, color = '#808080', opacity = 0.5,
                                           line = list(color = '#000000', width = 1)
                             )) %>% 
    layout(title = paste(" Plotly opls - ", title_text, sep = ""),
           xaxis = x_axis_settings_loading,
           yaxis = y_axis_settings_loading
    )
  
  combined_plotly <- subplot(plotly_opls, plotly_loadings, 
                             margin = c(0.05, 0.05, 0.01, 0.01),
                             titleX = TRUE,
                             titleY = TRUE
  ) %>% layout(showlegend = TRUE, title =  "")
  
  plotly_loadings_data %>% arrange(desc(p1)) %>% kable() %>% print()
  
  combined_plotly
  
  

}