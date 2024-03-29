#ANPC Lipidomics PCA quality control visualisation

lipids_pca_universal <- function(individual_multivariate_data, family_multivariate_data, multivariate_class, plot_label){
  #browser()
   scale_answer <- "blank"
  while(scale_answer != "UV" & scale_answer != "Pareto"){
    scale_answer <- dlgInput("What scaling do you want to apply to the PCA?", "UV/Pareto")$res
  }
  
   colour_answer <- "blank"
   groups_total <- nrow(unique(select(individual_multivariate_data, multivariate_class)))
   while(colour_answer == "blank" & length(colour_answer) < groups_total){
     colour_answer <- dlgInput(paste("You have ", paste(groups_total), " groups.  Please enter enough colours for the plot here and seperate by a comma, no spaces.  Colours should be in hexidecimal format"), 
                               "#1E90FF,#FF0000,#000000, etc")$res
   }
   
  colour_answer <- as.list(strsplit(colour_answer, ',')[[1]]) %>%
    lapply(function(x)gsub('\\s+', '',x)) %>% 
    unlist()
   
  lipid <- individual_multivariate_data %>% select(contains("(")) %>% colnames()
  lipid_class <- sub("\\(.*", "", lipid) %>% unique()
  lipid_class <- lipid_class[!grepl("sampleID", lipid_class)] %>% as_tibble()
  
  multivariate_data_list <- list(individual_multivariate_data, family_multivariate_data)
  
  pca_plot_list <- lapply(multivariate_data_list, function(func_list){
    #browser()
    multivariate_data <- func_list
    column_length <- multivariate_data %>% select(contains("(")) %>% ncol()
    
    if(column_length > 0){ 
      pca_x <- multivariate_data %>%  select(all_of(lipid)) %>% as.matrix()+1 
      pca_x <- log(pca_x)
      title_text <- "individual lipid species"
    }
    if(column_length == 0){ 
      pca_x <- multivariate_data %>%  select(all_of(lipid_class$value)) %>% as.matrix()+1 
      pca_x <- log(pca_x)
      title_text <- "lipid class"
    }
    
    pca_x[pca_x == 0] <- NA #remove all 0 values
    pca_x[is.infinite(pca_x)] <- NA #remove all infinite values
    min_value <- min(pca_x, na.rm = TRUE) # find the lowest value in the matrix
    pca_x[is.na(pca_x)] <- min_value # replace all NA, Inf, and 0 values with the lowest value in the matrix
    
    pca_class <- multivariate_data %>% select(multivariate_class) %>% as.matrix()
    pca_plot_label <- multivariate_data %>% select(plot_label) %>% as.matrix()
    pca_class[is.na(pca_class)] <- "none"
    sampleID <- multivariate_data %>% select(sampleID)
    
    pca_model <- pca(pca_x, scale = paste(scale_answer), center = TRUE)
    PC1 <- as.numeric(as.matrix(pca_model@t[,1]))
    PC2 <- as.numeric(as.matrix(pca_model@t[,2]))
    
    #browser()
    #produce plot_ly PCA plot
    plot_Val <- as_tibble(cbind(PC1, PC2))
    plot_Val$sampleID <- sampleID$sampleID
    plot_Val$sample_group <- c(pca_class)
    plot_Val$pca_plot_label <- c(pca_plot_label)
    # plot_Val_samples <- plot_Val %>% filter(sample_group == "sample")
    # plot_Val_ltr <- plot_Val %>% filter(sample_group == "LTR")
    
    plotly_loadings_data <- pca_model@p %>% as_tibble(rownames = "lipid") %>% rename(PC1 = V1, PC2 = V2)
    
    x_axis_settings_scores <- list(
      zeroline = TRUE,
      showline = TRUE,
      linecolor = toRGB("black"),
      linewidth = 2,
      showgrid = TRUE,
      title = paste("PC1 (", round(pca_model@Parameters$R2[1]*100,1), " %)", sep = "")
    )
    
    y_axis_settings_scores <- list(
      zeroline = TRUE,
      showline = TRUE,
      linecolor = toRGB("black"),
      linewidth = 2,
      showgrid = TRUE,
      title = paste("PC2 (", round(pca_model@Parameters$R2[2]*100,1), " %)", sep = "")
      )
    
    
    
    plotly_pca <- plot_ly(type = "scatter", mode = "markers", data = plot_Val, x = ~PC1, y = ~PC2, text =~pca_plot_label, color = ~sample_group, colors = c(colour_answer), 
                          marker = list(size = 7, 
                                        opacity = 0.5,
                                        line = list(
                                          color = '#000000',
                                          width = 1)
                                        )) %>% 
      layout(title = paste(project_name, " Plotly PCA - ", title_text, sep = ""),
             xaxis = x_axis_settings_scores,
             yaxis = y_axis_settings_scores
             )
    
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
    
    
    plotly_loadings <- plot_ly(type = "scatter", mode = "markers", data = plotly_loadings_data, x = ~PC1, y = ~PC2, text = ~lipid, 
                               color = "PCA loadings",
                               marker = list(size = 7, color = '#808080', opacity = 0.5,
                                             line = list(color = '#000000', width = 1)
                               )) %>% 
      layout(title = paste(project_name, " Plotly PCA - ", title_text, sep = ""),
             xaxis = x_axis_settings_loading,
             yaxis = y_axis_settings_loading,
             showlegend = FALSE
      )
    
    combined_plotly <- subplot(plotly_pca, plotly_loadings, 
                               margin = c(0.01, 0.01, 0.2, 0.01),
                               titleX = TRUE,
                               titleY = TRUE
                               ) %>% layout(showlegend = TRUE, title =  "Plotly PCA")
    
    pca_plot_list <- c(list(combined_plotly))
    #print(pca_plot_list)
    
  })
  
  pca_plot_list <- c(pca_plot_list,
                     list(scale_answer)
                     )
  pca_plot_list
  
  
  
}
