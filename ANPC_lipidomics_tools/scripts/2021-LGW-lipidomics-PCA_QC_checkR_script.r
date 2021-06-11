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

pca_p <- lipids_pca_ltr(final_individual_lipid_data, final_class_lipid_data, multivariate_class = "sample_class", plot_label = "sampleID", scaling = "option")

saveWidget(pca_p[[1]][[1]], file = paste(project_dir_html, "/", project_name, "_", user_name, "_QC_PCA_all_lipids.html", sep=""))# save plotly widget
browseURL(paste(project_dir_html, "/", project_name, "_", user_name, "_QC_PCA_all_lipids.html", sep="")) #open plotly widget in internet browser
saveWidget(pca_p[[2]][[1]], file = paste(project_dir_html, "/", project_name, "_", user_name, "_QC_PCA_lipid_class.html", sep=""))# save plotly widget
browseURL(paste(project_dir_html, "/", project_name, "_", user_name, "_QC_PCA_lipid_class.html", sep="")) #open plotly widget in internet browser


pca_check_status <- dlgInput("Check the PCA plots. Are you happy to continue? or do wish to change the scalling type?", "continue/change")$res
while(pca_check_status != "continue" & pca_check_status != "change"){
  pca_check_status <- dlgInput("Check the PCA plots. Are you happy to continue? or do wish to change the scalling type?", "continue/change")$res
}
}

scale_used <- pca_p[[3]][1]


