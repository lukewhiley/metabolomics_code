# create summed lipid class data


create_lipid_class_data_summed <- function(individual_lipid_data){
  #browser()
  lipid_class <- individual_lipid_data %>% select(contains("(")) %>% colnames() 
  lipid_class <- sub("\\(.*", "", lipid_class) %>% unique()
  lipid_class <- lipid_class[!grepl("sampleID", lipid_class)] %>% as_tibble()
  
  temp_class_data <- apply(lipid_class, 1, function(func_lipid_class){
  #browser()
  class_targets <- which(sub("\\(.*", "", colnames(individual_lipid_data)) == func_lipid_class) # find the columns in each lipid class
  temp_class_data <- individual_lipid_data %>% select(all_of(class_targets))  %>% mutate(rowsum = rowSums(.)) %>% select(rowsum)
  colnames(temp_class_data) <- func_lipid_class
  temp_class_data
}) %>% bind_cols()

temp_class_data <- cbind(individual_lipid_data$sampleID, temp_class_data) %>% rename("sampleID" = "individual_lipid_data$sampleID")
temp_class_data
}

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

######### PCA analysis ##########

# pca plot using metabomate

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
  
  #produce static PCA plot
  plot_Val <- as_tibble(cbind(PC1, PC2))
  plot_Val$sampleID <- sampleID$sampleID
  plot_Val$sample_group <- c(pca_class)
  plot_Val$pca_plot_label <- c(pca_plot_label)
  
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
  
  plotly_pca <- plot_ly(type = "scatter", plot_Val, x = ~PC1, y = ~PC2, text =~pca_plot_label, color = ~sample_group, colors = plot_colours, marker = list(size = 6)) %>% 
    layout(title = paste(project_name, " Plotly PCA - ", title_text, sep = ""))
  
  plotly_loadings <- plot_ly(type = "scatter", plotly_loadings_data, x = ~PC1, y = ~PC2, text = ~lipid, marker = list(color = "black")) %>% 
    layout(title = paste(project_name, " Plotly PCA - ", title_text, sep = ""))
  
  combined_plotly <- subplot(plotly_pca, plotly_loadings, 
                             margin = c(0.01, 0.01, 0.2, 0.01)) %>% layout(showlegend = FALSE, title =  "Plotly PCA")
  #pca_plotly_list <- list(plotly_pca, plotly_loadings, combined_plotly)
  #print(pca_plotly_list) 
  
  pca_plot_list <- c(list(#pca_plot, 
                          #pca_plot_loadings, 
                          #plotly_pca, 
                          #plotly_loadings, 
                          combined_plotly))
  print(pca_plot_list)
  
 })
  
    pca_plot_list
  
  # interactive_pca <- subplot(pca_plot_list[[5]],
  #                               #pca_plot_list[[10]], 
  #                               nrows = 2, margin = c(0.1, 0.1, 0.1, 0.01)) %>%
  #   layout(annotations = list(
  #     list(x = 0.1, y = 1.025, text = "Individual lipids scores plot", showarrow = F, xref='paper', yref='paper'),
  #     list(x = 0.85 , y = 1.025, text = "Individual lipids loadings plot", showarrow = F, xref='paper', yref='paper'),
  #     list(x = 0.12 , y = 0.45, text = "Lipid family scores plot", showarrow = F, xref='paper', yref='paper'),
  #     list(x = 0.85 , y = 0.45, text = "Lipid family loadings plot", showarrow = F, xref='paper', yref='paper'))
  #   )
  
  #pca_output <- c(pca_plot_list)
  
 # pca_output
  }
  #pca_plot_list
#}

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

###### OPLS-DA function #####

lipids_opls <- function(individual_multivariate_data, family_multivariate_data, multivariate_class){
  #browser()
  lipid_class <- as_tibble(c("CE",  "CER", "DAG",  "DCER", "FFA", "HCER", "LCER", "LPC", "LPE", "LPG", "PC", "PE", "PG", "PI", "PS", "SM", "TAG"))
  multivariate_data_list <- list(individual_multivariate_data, family_multivariate_data)
  pca_plot_list <- NULL
  
  for(idx_list in 1:2){
    multivariate_data <- multivariate_data_list[[idx_list]]

  if(length(grep("TAG", colnames(multivariate_data))) > 1){ 
    pca_x <- multivariate_data %>%  select(all_of(lipid)) %>% as.matrix()
  }
  
  if(length(grep("TAG", colnames(multivariate_data))) == 1){ 
    pca_x <- multivariate_data %>%  select(all_of(lipid_class$value)) %>% as.matrix()
  }

  pca_x[pca_x == 0] <- NA #remove all 0 values
  pca_x[is.infinite(pca_x)] <- NA #remove all infinite values
  min_value <- min(pca_x, na.rm = TRUE) # find the lowest value in the matrix
  pca_x[is.na(pca_x)] <- min_value # replace all NA, Inf, and 0 values with the lowest value in the matrix
  
  opls_x <- pca_x
  opls_y <- multivariate_data %>% select(multivariate_class) %>% as.matrix()
  sampleID <- multivariate_data %>% select(sampleID)
  
  opls_model <- opls(log10(opls_x+1), 
                     opls_y, 
                     center = TRUE, 
                     scale = "UV", 
                     cv = list(method = "k-fold_stratified", k = 7, split = 2/3))
  
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
    scale_color_manual(values=c(plot_colours))+
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
  }
  opls_plot_list
}

