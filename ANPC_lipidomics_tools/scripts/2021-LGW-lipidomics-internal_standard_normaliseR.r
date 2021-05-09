################################################
###### Internal standard normalization #########
#################################################

# this section normalizes each SRM lipid target with the appropriate internal standard
# Requires a template guide with internal standard transition for each target lipid SRM transition

dlg_message("Time for normalization using the internal standards :-)", type = 'ok')
dlg_message("REQUIRES - a template csv listing each lipid target and the assigned internal standard", type = 'ok')

#import transition report 3
filtered_data <- individual_lipid_data_tic_filtered %>% filter(!grepl("conditioning", sampleID))
filtered_data[is.na(filtered_data)] <- 0

dlg_message("Please select this template file now.", type = 'ok')

sil_target_list <- read_csv(file = file.choose(.)) %>% clean_names
lipid_data <- filtered_data %>% select(-sampleID, - plate_id) %>% select(!contains("SIL"))
sil_data <- filtered_data %>% select(-sampleID, - plate_id) %>% select(contains("SIL"))

normalisation_factor <- lapply(filtered_data$sampleID, function(sampleNORM){
  #browser()
  func_data <- filtered_data %>% filter(sampleID == sampleNORM) %>% select(contains("SIL")) %>% rowSums(na.rm = TRUE)
}) %>% unlist

ratio_data <- apply(as_tibble(colnames(lipid_data)), 1, function(LIPID){
  #browser()
  func_data <- lipid_data %>% select(all_of(LIPID))
  sil_to_use <- sil_target_list$note[which(sil_target_list$precursor_name==LIPID)]
  func_data_sil <- sil_data %>% select(sil_to_use)
  normalised_data <- func_data/func_data_sil
  #normalised_data <- (func_data/((func_data+func_data_sil)))/normalisation_factor
  #normalised_data <- (func_data/func_data_sil)+normalisation_factor
  #normalised_data <- func_data/(func_data+normalisation_factor)
  #normalised_data <- func_data/func_data_sil
  normalised_data
}) %>% bind_cols() %>% add_column(filtered_data$sampleID, filtered_data$plate_id, .before = 1)

colnames(ratio_data) <- c("sampleID", "plateID", colnames(lipid_data))

ltr_rsd <- apply(as_tibble(colnames(lipid_data)), 1, function(RSD){
  #browser()
  func_data <- ratio_data %>% filter(grepl("LTR", sampleID)) %>% select(all_of(RSD))
  (sd(func_data$value)*100)/mean(func_data$value)
}) %>% as_tibble() %>% add_column(colnames(lipid_data), .before = 1)

colnames(ltr_rsd) <- c("lipid", "RSD")
dlg_message(paste("number of feature ratios with with an LTR RSD of <30% =", length(which(ltr_rsd$RSD < 30))), type = 'ok')
dlg_message(paste("number of feature ratios with with an LTR RSD of <20% =", length(which(ltr_rsd$RSD < 20))), type = 'ok')
dlg_message(paste("number of feature ratios with with an LTR RSD of <15% =", length(which(ltr_rsd$RSD < 15))), type = 'ok')
dlg_message(paste("number of feature ratios with with an LTR RSD of <10% =", length(which(ltr_rsd$RSD < 10))), type = 'ok')

lipid_keep_list <- ltr_rsd %>% filter(RSD < 30)

final_dataset <- ratio_data %>% select(sampleID, plateID, all_of(lipid_keep_list$lipid))
final_dataset[is.na(final_dataset)] <- 0

# visualisation of normalied data
# first - produce a plot of all normalized features to see if there are any overall trends in the data

total_summed_ratio <- apply(final_dataset %>% select(sampleID), 1, function(summedTIC){
  #browser()
  temp_data <- final_dataset %>% filter(sampleID == summedTIC) %>% select(-sampleID, -plateID) %>% rowSums(na.rm = TRUE)
}) %>% c() %>% as_tibble() %>%  add_column(final_dataset$sampleID, .before = 1) %>% 
  rename(summed_TIC = value, sampleID = "final_dataset$sampleID")
total_summed_ratio$sample_idx <- c(1:nrow(total_summed_ratio))

sd(total_summed_ratio$summed_TIC*100)/mean(total_summed_ratio$summed_TIC)

total_summed_ratio$sample <- "sample"
total_summed_ratio$sample[grep("LTR", total_summed_ratio$sampleID)] <- "LTR"

p <- plot_ly(type = "scatter", total_summed_ratio, x = ~sample_idx, y = ~summed_TIC, text = ~sampleID, color = ~sample, colors = c("red", "lightblue"));p

plate_number <- unique(plate_id) %>% substr(14,14) %>% unique()
plate_idx <- lapply(unique(plateid), function(plateID){grep(plateID, total_summed_sil$sampleID)[1]}) %>% unlist()
for (idx_line in 2:length(plate_idx)){
  p <- add_trace(p, x = plate_idx[idx_line], type = 'scatter', mode = 'lines', color = paste("plate_", plate_number[idx_line], sep=""), line = list(color = "grey", dash = "dash"), showlegend = FALSE)
}
normalized_check_p <- p

#create html widget and display it in the users internet browser
if(!dir.exists(paste(project_dir, "/html_files", sep=""))){dir.create(paste(project_dir, "/html_files", sep=""))} # create a new directory to store html widgets
saveWidget(normalized_check_p, file = paste(project_dir, "/html_files/",project_name, "_", user_name, "_normalized_check_plot.html", sep=""))# save plotly widget
browseURL(paste(project_dir, "/html_files/",project_name, "_", user_name, "_normalized_check_plot.html", sep="")) #open plotly widget in internet browser

dlg_message("Check plot for summed all normilized features. Press OK to continue", type = 'ok')



# second - produce a plot of normalized features, summed by class to see if there are any overall trends in the lipid class data

dlg_message("Now we are going to look at the data summed by lipid class", type = 'ok')

final_individual_lipid_data <- final_dataset

final_class_lipid_data <- create_lipid_class_data_summed(final_individual_lipid_data)

lipid_class_list <- final_individual_lipid_data %>% select(contains("(")) %>% colnames() 
lipid_class_list <- sub("\\(.*", "", lipid_class_list) %>% unique()
lipid_class_list <- lipid_class_list[!grepl("sampleID", lipid_class_list)] %>% as_tibble()

#add ltr TRUE/FALSE column

final_class_lipid_data$is_ltr <- "sample"
final_class_lipid_data$is_ltr[grep("LTR", final_class_lipid_data$sampleID)] <- "LTR"

plotlist <- apply(lipid_class_list %>% select(value), 1, function(lipidClass){
  #browser()
  plot_data <- final_class_lipid_data %>% select("sampleID", "is_ltr", all_of(lipidClass)) %>% rename(ms_response = value) 
  plate_id <- str_extract(plot_data$sampleID, "PLIP.*")
  plate_id <- substr(plate_id, 0,15)
  plot_data$sample_index <- paste(plate_id, sub(".*\\_", "", plot_data$sampleID), sep="_")
  plot_data <- plot_data %>% arrange(sample_index)
  plot_data$idx <- 1:nrow(plot_data)
  
  plot_data_ltr <- plot_data %>% filter(is_ltr == "LTR")
  plot_data <- plot_data %>% filter(is_ltr == "sample")
  
  plate_idx <- lapply(unique(plate_id), function(plateID){
    grep(plateID, plot_data$sampleID)[1]
  }) %>% unlist
  
  p <- plot_ly(type = 'scatter', mode   = 'markers', plot_data, x = ~idx, y = ~log(ms_response+1), text = ~sampleID, color = ~is_ltr,  colors = c("red", "lightblue3"), showlegend = FALSE) %>% 
    layout(xaxis = list(title = paste(lipidClass)))
  
  p <- add_trace(p, data = plot_data_ltr, x = ~idx, y = ~log(ms_response+1), text = ~sampleID, color = ~is_ltr,  colors = c("red", "lightblue3"), showlegend = FALSE)
  
  plate_number <- unique(plate_id) %>% substr(14,14)
  
  plot_limits <- log(c(min(plot_data$ms_response+1), max(plot_data$ms_response+1)))
  
  for (idx_line in 2:length(plate_idx)){
    p <- p %>% add_segments(x = plate_idx[idx_line], xend = plate_idx[idx_line], y = plot_limits[1], yend = plot_limits[2], line = list(color = "grey", dash = "dash"), showlegend = FALSE, text = plateid[plate_idx[idx_line]])
  }
  p
})

normalized_check_class_p <- subplot(plotlist, nrows = 4, titleX = TRUE, margin = c(0.01,0.01,0.05,0.05))

if(!dir.exists(paste(project_dir, "/html_files", sep=""))){dir.create(paste(project_dir, "/html_files", sep=""))} # create a new directory to store html widgets
saveWidget(normalized_check_class_p, file = paste(project_dir, "/html_files/",project_name, "_", user_name, "_normalized_check_class_plot.html", sep=""))# save plotly widget
browseURL(paste(project_dir, "/html_files/",project_name, "_", user_name, "_normalized_check_class_plot.html", sep="")) #open plotly widget in internet browser

dlg_message("Check plot for summed lipid class normilized features. Press OK to continue", type = 'ok')

