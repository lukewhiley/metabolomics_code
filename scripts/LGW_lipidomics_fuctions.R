######### PCA analysis ##########

# pca plot using metaboanalyst

lipids_pca <- function(multivariate_data){
  if(length(grep("TAG", colnames(multivariate_data))) > 1){ 
  pca_x <- multivariate_data %>%  select(c(lipid)) %>% as.matrix()
}

if(length(grep("TAG", colnames(multivariate_data))) == 1){ 
   pca_x <- multivariate_data %>%  select(c(lipid_family$lipid_class)) %>% as.matrix()
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
  
  pca_plot <- ggplot() +
    geom_vline(xintercept = 0, colour="black", linetype = "longdash", alpha = 0.4)+
    geom_point(data=plot_Val, aes(x =  PC1, y =  PC2, color = class),  size = 6.0)+
    theme_bw() +
    scale_color_manual(values=plot_two_group_colours)+
    xlab("PC1")+
    ylab("PC2")+
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size=30)) +
    theme(axis.text.y=element_text(size = 22, margin = margin(t = 0, r = 0, b = 0, l = 2))) +
    theme(axis.text.x=element_text(size = 22)) +
    theme(axis.title = element_text(size = 25)) +
    theme(legend.text = element_text(size = 18)) +
    theme(legend.title  = element_text(size = 18))+
    theme(legend.position = "bottom") + 
    ggtitle(paste(project_name, " PCA - individual lipid species", sep = ""))
  
  if(label_sampleIDs == TRUE){
    pca_plot <- pca_plot + geom_text(data=plot_Val, aes(x =  PC1, y =  PC2, label = sampleID))
  }

  
  pca_plot_loadings <- plotload_cat(pca_model, pc = c(1,2), an = list("","","")) +
    scale_colour_manual(values='black') +
    theme_bw() +
    xlab("PC1")+
    ylab("PC2")+
    ggtitle("PCA Loadings - individual lipid species")+
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size=30)) +
    theme(axis.text.y=element_text(size = 22, margin = margin(t = 0, r = 0, b = 0, l = 2))) +
    theme(axis.text.x=element_text(size = 22)) +
    theme(axis.title = element_text(size = 25)) +
    theme(legend.position = "none"); pca_plot_loadings
  
  pca_plot_list <- list(pca_plot, pca_plot_loadings)
  print(pca_plot_list)
  

}


#interactive plotly PCA plots

lipids_pca_plotly <- function(multivariate_data){
    if(length(grep("TAG", colnames(multivariate_data))) > 1){ 
      pca_x <- multivariate_data %>%  select(c(lipid)) %>% as.matrix()
    }
    
    if(length(grep("TAG", colnames(multivariate_data))) == 1){ 
      pca_x <- multivariate_data %>%  select(c(lipid_family$lipid_class)) %>% as.matrix()
    }
  pca_class <- multivariate_data %>% select(class) %>% as.matrix()
  sampleID <- multivariate_data %>% select(sampleID)
  pca_model <- pca(pca_x, scale = "UV", center = TRUE)
  
  plotly_loadings_data <- pca_model@p %>% as_tibble(rownames = "lipid") %>% rename(PC1 = V1, PC2 = V2)
  
  plotly_pca <- plot_ly(plot_Val, x = ~PC1, y = ~PC2, text = ~sampleID, color = ~class, colors = plot_two_group_colours, marker = list(size = 12)); plotly_pca
  plotly_loadings <- plot_ly(plotly_loadings_data, x = ~PC1, y = ~PC2, text = ~lipid, colors = "black");plotly_loadings
  combined_plotly <- subplot(plotly_pca, plotly_loadings) %>% layout(showlegend = FALSE)
  pca_plotly_list <- list(plotly_pca, plotly_loadings, combined_plotly)
  print(pca_plotly_list)
}
