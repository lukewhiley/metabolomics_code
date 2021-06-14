# produces a PCA of the normalised sampels

#source PCA QC function from guthub
lipidomics_PCA_QC_function <- GET(url = "https://raw.githubusercontent.com/lukewhiley/metabolomics_code/main/ANPC_lipidomics_tools/functions/2021-LGW-lipidomics-PCA_QC_checkR_function.r") %>% content(as = "text")
eval(parse(text = lipidomics_PCA_QC_function), envir = .GlobalEnv)
rm(lipidomics_PCA_QC_function)

#label data
final_individual_lipid_data$sample_class <- "sample"
final_individual_lipid_data$sample_class[grep("LTR", final_individual_lipid_data$sampleID)] <- "LTR"

final_class_lipid_data$sample_class <- "sample"
final_class_lipid_data$sample_class[grep("LTR", final_class_lipid_data$sampleID)] <- "LTR"

#run function
pca_check_status <- "change"
while(pca_check_status == "change"){

  plotlist <- apply(lipid_class_list %>% select(value), 1, function(lipidClass){
    #browser()
    PCA_class_plot_data <- final_individual_lipid_data %>% select(sampleID, sample_class, starts_with(paste0(lipidClass, "(")))
    
  if(ncol(PCA_class_plot_data) > 3){
    pca_p <- lipids_pca_ltr(PCA_class_plot_data, final_class_lipid_data, multivariate_class = "sample_class", plot_label = "sampleID", scaling = "UV")
    pca_p[[1]][[1]]
  }
  })
  
  plotlist <- plotlist[-which(sapply(plotlist, is.null))]
  
  PCA_class_plot <- subplot(plotlist, nrows = 4, titleX = FALSE, margin = c(0.015,0.015, 0.05,0.05))
  
  saveWidget(PCA_class_plot, file = paste(project_dir_html, "/", project_name, "_", user_name, "_PCA_class_check_plot.html", sep=""))# save plotly widget
  browseURL(paste(project_dir_html, "/", project_name, "_", user_name, "_PCA_class_check_plot.html", sep="")) #open plot_ly widget in internet browser
}


