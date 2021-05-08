# produces a PCA of the normalised sampels



individual_lipid_data$sample_class <- "sample"
individual_lipid_data$sample_class[grep("LTR", individual_lipid_data$sampleID)] <- "LTR"

class_lipid_data$sample_class <- "sample"
class_lipid_data$sample_class[grep("LTR", class_lipid_data$sampleID)] <- "LTR"

plot_colours <- c("red", "lightblue3")
label_sampleIDs <- FALSE

pca_check_status <- "change"
while(pca_check_status == "change"){

pca_p <- lipids_pca(individual_lipid_data, class_lipid_data, multivariate_class = "sample_class", plot_label = "sampleID")

if(!dir.exists(paste(project_dir, "/html_files", sep=""))){dir.create(paste(project_dir, "/html_files", sep=""))} # create a new directory to store html widgets
saveWidget(pca_p[[1]][[1]], file = paste(project_dir, "/html_files/",project_name, "_", user_name, "_QC_PCA_all_lipids.html", sep=""))# save plotly widget
browseURL(paste(project_dir, "/html_files/",project_name, "_", user_name, "_QC_PCA_all_lipids.html", sep="")) #open plotly widget in internet browser
saveWidget(pca_p[[2]][[1]], file = paste(project_dir, "/html_files/",project_name, "_", user_name, "_QC_PCA_lipid_class.html", sep=""))# save plotly widget
browseURL(paste(project_dir, "/html_files/",project_name, "_", user_name, "_QC_PCA_lipid_class.html", sep="")) #open plotly widget in internet browser


pca_check_status <- dlgInput("Check the PCA plots. Are you happy to continue? or do wish to change the scalling type?", "continue/change")$res
while(tic_check_status != "continue" & sil_check_status != "change"){
  pca_check_status <- dlgInput("Check the PCA plots. Are you happy to continue? or do wish to change the scalling type?", "continue/change")$res
}
}