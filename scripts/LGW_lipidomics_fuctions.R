# create family data


create_family_data_summed <- function(individual_lipid_data){
  lipid_class <- as_tibble(c("CE",  "CER", "DAG",  "DCER", "FFA", "HCER", "LCER", "LPC", "LPE", "LPG", "PC", "PE", "PG", "PI", "PS", "SM", "TAG"))
  #lipid_family <- individual_lipid_data %>% select(all_of(lipid_class)) %>% unique()
  apply(lipid_class, 1, function(func_family){
  #browser()
  family_targets <- which(sub("\\(.*", "", colnames(individual_lipid_data)) == func_family) # find the columns in each family
  temp_family_data <- individual_lipid_data %>% select(family_targets)  %>% mutate(rowsum = rowSums(.)) %>% select(rowsum)
  colnames(temp_family_data) <- func_family
  temp_family_data
}) %>% bind_cols()

#colnames(temp_family_data) <- c(lipid_family$lipid_class)
}


######### PCA analysis ##########

# pca plot using metaboanalyst

lipids_pca <- function(multivariate_data, multivariate_class){
  if(length(grep("TAG", colnames(multivariate_data))) > 1){ 
    pca_x <- multivariate_data %>%  select(c(lipid)) %>% as.matrix()
    title_text <- "individual lipid species"
  }
  if(length(grep("TAG", colnames(multivariate_data))) == 1){ 
    pca_x <- multivariate_data %>%  select(c(lipid_family$lipid_class)) %>% as.matrix()
    title_text <- "lipid family"
  }
  
  pca_x[pca_x == 0] <- NA #remove all 0 values
  pca_x[is.infinite(pca_x)] <- NA #remove all infinite values
  min_value <- min(pca_x, na.rm = TRUE) # find the lowest value in the matrix
  pca_x[is.na(pca_x)] <- min_value # replace all NA, Inf, and 0 values with the lowest value in the matrix
  
  pca_class <- multivariate_data %>% select(multivariate_class) %>% as.matrix()
  pca_class[is.na(pca_class)] <- "none"
  sampleID <- multivariate_data %>% select(sampleID)
  
  pca_model <- pca(pca_x, scale = "UV", center = TRUE)
  PC1 <- as.numeric(as.matrix(pca_model@t[,1]))
  PC2 <- as.numeric(as.matrix(pca_model@t[,2]))
  
  #produce PCA plot
  plot_Val <- as_tibble(cbind(PC1, PC2))
  plot_Val$sampleID <- sampleID$sampleID
  plot_Val$sample_group <- c(pca_class)
  
  pca_plot <- ggplot() +
    geom_vline(xintercept = 0, colour="black", linetype = "longdash", alpha = 0.4)+
    geom_point(data=plot_Val, aes(x =  PC1, y =  PC2, color = as.factor(sample_group)),  size = 3.0)+
    theme_bw() +
    scale_color_manual(values=c(plot_colours))+
    xlab("PC1")+
    ylab("PC2")+
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size=30)) +
    theme(axis.text.y=element_text(size = 22, margin = margin(t = 0, r = 0, b = 0, l = 2))) +
    theme(axis.text.x=element_text(size = 22)) +
    theme(axis.title = element_text(size = 25)) +
    theme(legend.text = element_text(size = 18)) +
    theme(legend.title  = element_text(size = 18))+
    theme(legend.position = "right") + 
    labs(color = "Sample group") +
    ggtitle(paste(project_name, " PCA - ", title_text, sep = ""))
  
  if(label_sampleIDs == TRUE){
    pca_plot <- pca_plot + geom_text(data=plot_Val, aes(x =  PC1, y =  PC2, label = sampleID))
  }

  
  pca_plot_loadings <- plotload_cat(pca_model, pc = c(1,2), an = list("","","")) +
    scale_colour_manual(values='black') +
    theme_bw() +
    xlab("PC1")+
    ylab("PC2")+
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size=30)) +
    theme(axis.text.y=element_text(size = 22, margin = margin(t = 0, r = 0, b = 0, l = 2))) +
    theme(axis.text.x=element_text(size = 22)) +
    theme(axis.title = element_text(size = 25)) +
    theme(legend.position = "none") +
    ggtitle(paste(project_name, " PCA - ", title_text, sep = ""))
  
  #pca_plot_list <- list(pca_plot, pca_plot_loadings)
  #print(pca_plot_list)
  
  plotly_loadings_data <- pca_model@p %>% as_tibble(rownames = "lipid") %>% rename(PC1 = V1, PC2 = V2)
  
  plotly_pca <- plot_ly(type = "scatter", plot_Val, x = ~PC1, y = ~PC2, text =~sampleID, color = ~sample_group, colors = plot_colours, marker = list(size = 12)) %>% 
    layout(title = paste(project_name, " Plotly PCA - ", title_text, sep = ""))
  
  plotly_loadings <- plot_ly(type = "scatter", plotly_loadings_data, x = ~PC1, y = ~PC2, text = ~lipid, marker = list(color = "black")) %>% 
    layout(title = paste(project_name, " Plotly PCA - ", title_text, sep = ""))
  
  combined_plotly <- subplot(plotly_pca, plotly_loadings) %>% layout(showlegend = FALSE, title =  "Plotly PCA")
  #pca_plotly_list <- list(plotly_pca, plotly_loadings, combined_plotly)
  #print(pca_plotly_list) 
  
  pca_plot_list <- list(pca_plot, pca_plot_loadings, plotly_pca, plotly_loadings, combined_plotly)
  pca_plot_list
}


#interactive plotly PCA plots

lipids_pca_plotly <- function(multivariate_data, multivariate_class){
 
  if(length(grep("TAG", colnames(multivariate_data))) > 1){ 
    pca_x <- multivariate_data %>%  select(c(lipid)) %>% as.matrix()
    title_text <- "individual lipid species"
  }
  
  if(length(grep("TAG", colnames(multivariate_data))) == 1){ 
    pca_x <- multivariate_data %>%  select(c(lipid_family$lipid_class)) %>% as.matrix()
    title_text <- "lipid family"
  }
  
  pca_class <- multivariate_data %>% select(class) %>% as.matrix()
  sampleID <- multivariate_data %>% select(sampleID)
  
  pca_model <- pca(pca_x, scale = "UV", center = TRUE)
  PC1 <- as.numeric(as.matrix(pca_model@t[,1]))
  PC2 <- as.numeric(as.matrix(pca_model@t[,2]))
  
  #produce PCA plot
  plot_Val <- as_tibble(cbind(PC1, PC2))
  plot_Val$sampleID <- sampleID$sampleID
  plot_Val$class <- pca_class
  
  plotly_loadings_data <- pca_model@p %>% as_tibble(rownames = "lipid") %>% rename(PC1 = V1, PC2 = V2)

  plotly_pca <- plot_ly(plot_Val, x = ~PC1, y = ~PC2, text =~sampleID, color = ~class, colors = plot_two_group_colours, marker = list(size = 12)) %>% 
    layout(title = paste(project_name, " Plotly PCA - ", title_text, sep = ""))
  plotly_loadings <- plot_ly(plotly_loadings_data, x = ~PC1, y = ~PC2, text = ~lipid, marker = list(color = "black")) %>% 
    layout(title = paste(project_name, " Plotly PCA - ", title_text, sep = ""))
  combined_plotly <- subplot(plotly_pca, plotly_loadings) %>% layout(showlegend = FALSE, title =  "Plotly PCA")
  pca_plotly_list <- list(plotly_pca, plotly_loadings, combined_plotly)
  print(pca_plotly_list) 
}


###### OPLS-DA function #####

lipids_opls <- function(multivariate_data){
  # browser()
  if(length(grep("TAG", colnames(multivariate_data))) > 1){ 
    pca_x <- multivariate_data %>%  select(c(lipid)) %>% as.matrix()
  }
  
  if(length(grep("TAG", colnames(multivariate_data))) == 1){ 
    pca_x <- multivariate_data %>%  select(c(lipid_family$lipid_class)) %>% as.matrix()
  }

  opls_x <-  multivariate_data %>% select(c(lipid)) %>% as.matrix() 
  opls_y <- multivariate_data %>% select(class) %>% as.matrix()
  sampleID <- multivariate_data %>% select(sampleID)
  
  opls_model <- opls(log10(opls_x+1), opls_y, center = TRUE, scale = "UV")
  plot_t_pred <- as.numeric(as.matrix(opls_model@t_pred))
  plot_t_orth <- as.numeric(as.matrix(opls_model@t_orth))
  
  #produce OPLS-DA plot
  plot_Val <- as_tibble(cbind(plot_t_pred, plot_t_orth))
  plot_Val$class <- opls_y[,'class']
  plot_Val$sampleID <- sampleID$sampleID
  
  opls_plot <- ggplot() +
    geom_vline(xintercept = 0, colour="black", linetype = "longdash", alpha = 0.4)+
    geom_point(data=plot_Val, aes(x =  plot_t_pred, y =  plot_t_orth, color = class),  size = 6.0)+
    theme_bw() +
    scale_color_manual(values=plot_two_group_colours)+
    xlab("t_pred")+
    ylab("t_orth")+
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size=30)) +
    theme(axis.text.y=element_text(size = 22, margin = margin(t = 0, r = 0, b = 0, l = 2))) +
    theme(axis.text.x=element_text(size = 22)) +
    theme(axis.title = element_text(size = 25)) +
    theme(legend.position = 'none') +
    ggtitle(paste(project_name, " OPLS-DA training model - individual lipid species", sep = ""))
  
  if(label_sampleIDs == TRUE){
    opls_plot <- opls_plot + geom_text(data=plot_Val, aes(x =  plot_t_pred, y =  plot_t_orth, label = sampleID))
  }
  
  opls_plot_loadings <- eruption(opls_model)
  
  opls_plot_list <- list(opls_plot, opls_plot_loadings)
  print(opls_plot_list)
}
