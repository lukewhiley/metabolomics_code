# produces a PCA of the normalised sampels

individual_lipid_data$sample_class <- "sample"
individual_lipid_data$sample_class[grep("LTR", individual_lipid_data$sampleID)] <- "LTR"

class_lipid_data$sample_class <- "sample"
class_lipid_data$sample_class[grep("LTR", class_lipid_data$sampleID)] <- "LTR"

plot_colours <- c("red", "lightblue3")
label_sampleIDs <- FALSE

lipids_pca(individual_lipid_data, class_lipid_data, multivariate_class = "sample_class", plot_label = "sampleID")
