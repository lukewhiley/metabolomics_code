#ANPC Lipidomics PCA quality control visualisation

lipids_pca <- function(individual_multivariate_data, family_multivariate_data, multivariate_class, plot_label){
  
  temp_answer <- "blank"
  while(temp_answer != "UV" & temp_answer != "Pareto"){
    temp_answer <- dlgInput("What scaling do you want to apply to the PCA?", "UV/Pareto")$res
  }
  
  lipid <- individual_multivariate_data %>% select(contains("(")) %>% colnames()
  lipid_class <- sub("\\(.*", "", lipid) %>% unique()
  lipid_class <- lipid_class[!grepl("sampleID", lipid_class)] %>% as_tibble()
  
  multivariate_data_list <- list(individual_multivariate_data, family_multivariate_data)
  
  pca_plot_list <- lapply(multivariate_data_list, function(func_list){
    multivariate_data <- func_list
    column_length <- multivariate_data %>% select(contains("(")) %>% ncol()
    
    if(column_length > 0){ 
      pca_x <- multivariate_data %>%  select(all_of(lipid)) %>% as.matrix()
      title_text <- "individual lipid species"
    }
    if(column_length == 0){ 
      pca_x <- multivariate_data %>%  select(all_of(lipid_class$value)) %>% as.matrix()
      title_text <- "lipid family"
    }
    
    pca_x[pca_x == 0] <- NA #remove all 0 values
    pca_x[is.infinite(pca_x)] <- NA #remove all infinite values
    min_value <- min(pca_x, na.rm = TRUE) # find the lowest value in the matrix
    pca_x[is.na(pca_x)] <- min_value # replace all NA, Inf, and 0 values with the lowest value in the matrix
    
    pca_class <- multivariate_data %>% select(multivariate_class) %>% as.matrix()
    pca_plot_label <- multivariate_data %>% select(plot_label) %>% as.matrix()
    pca_class[is.na(pca_class)] <- "none"
    sampleID <- multivariate_data %>% select(sampleID)
    
    pca_model <- pca(pca_x, scale = paste(temp_answer), center = TRUE)
    PC1 <- as.numeric(as.matrix(pca_model@t[,1]))
    PC2 <- as.numeric(as.matrix(pca_model@t[,2]))
    
    #produce plot_ly PCA plot
    plot_Val <- as_tibble(cbind(PC1, PC2))
    plot_Val$sampleID <- sampleID$sampleID
    plot_Val$sample_group <- c(pca_class)
    plot_Val$pca_plot_label <- c(pca_plot_label)
    
    plotly_loadings_data <- pca_model@p %>% as_tibble(rownames = "lipid") %>% rename(PC1 = V1, PC2 = V2)
    
    plotly_pca <- plot_ly(type = "scatter", plot_Val, x = ~PC1, y = ~PC2, text =~pca_plot_label, color = ~sample_group, colors = plot_colours, marker = list(size = 6)) %>% 
      layout(title = paste(project_name, " Plotly PCA - ", title_text, sep = ""))
    
    plotly_loadings <- plot_ly(type = "scatter", plotly_loadings_data, x = ~PC1, y = ~PC2, text = ~lipid, marker = list(color = "black")) %>% 
      layout(title = paste(project_name, " Plotly PCA - ", title_text, sep = ""))
    
    combined_plotly <- subplot(plotly_pca, plotly_loadings, 
                               margin = c(0.01, 0.01, 0.2, 0.01)) %>% layout(showlegend = FALSE, title =  "Plotly PCA")
    
    pca_plot_list <- c(list(combined_plotly))
    print(pca_plot_list)
    
  })
  
  pca_plot_list
}