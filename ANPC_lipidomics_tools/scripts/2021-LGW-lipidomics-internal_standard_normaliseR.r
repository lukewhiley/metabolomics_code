################################################
###### Internal standard normalization #########
#################################################

# this section normalizes each SRM lipid target with the appropriate internal standard
# Requires a template guide with internal standard transition for each target lipid SRM transition

dlg_message("Time for normalization using the internal standards :-)", type = 'ok')
dlg_message("REQUIRES - a template csv listing each lipid target and the assigned internal standard. Column headings required are 'Precursor Name' containing the lipid target  and 'Note' containing the IS name", type = 'ok')

#import transition report 3
filtered_data <- individual_lipid_data_tic_intensity_filtered %>% filter(!grepl("conditioning", sampleID))
filtered_data[is.na(filtered_data)] <- 0

dlg_message("Please select this template file now.", type = 'ok')

temp_answer <- "change"
while(temp_answer == "change"){
sil_target_list <- read_csv(file = file.choose(.)) %>% clean_names
lipid_data <- filtered_data %>% select(-sampleID, - plate_id) %>% select(!contains("SIL"))
sil_data <- filtered_data %>% select(-sampleID, - plate_id) %>% select(contains("SIL"))

# this section checks each of the SIL IS used in the target list template in the LTRS. It evaluates if:
##  a: is the internal standard present in the LTR samples? Some batches of IS do not contain every IS availible. This alos prevents user error if the IS batch has not been made correctly.
##  b: the RSD of the internal standard signal intensity in LTRs. If the SIL IS is suitable for use in the dataset its signal should be stable in the LTRs. (raises a warning for >30%)

dlg_message("Checks to see if all internal standards are present in the SIL internal standard mix", type = 'ok')

sil_data_check <- individual_lipid_data_sil_tic_filtered %>% select(sampleID, plate_id, contains("SIL")) %>% filter(grepl("LTR", sampleID))

sil_list <- sil_target_list %>% filter(grepl("SIL", note)) %>% select(note) %>% unique() 

sil_sum <- lapply(sil_list$note, function(FUNC_SIL){
  #browser()
  temp_func_data_sum <- sil_data_check %>% select(all_of(FUNC_SIL)) %>% as.matrix() %>% sum() %>% log()
}) %>% c() %>% unlist() %>% as_tibble() %>% rename(SIL_SUM = value) %>% add_column(sil_list, .before = 1) %>% arrange(SIL_SUM)

sil_sum_q1 <- quantile(sil_sum$SIL_SUM, 0.25) %>% as.numeric()
inter_quantile_range <- as.numeric(quantile(sil_sum$SIL_SUM, 0.75)) - as.numeric(quantile(sil_sum$SIL_SUM, 0.25))
sil_sum_lower_threshold <- sil_sum_q1 - inter_quantile_range

#create a list of IS that fail the test
sil_list_warning <- sil_sum$note[which(sil_sum$SIL_SUM < sil_sum_lower_threshold)]

dlg_message(paste("Warning! These internal standards have a very low signal and may not be present in the mixture:", 
                  paste(sil_list_warning, collapse = ", "), 
                  " - double check skyline!"), 
            type = 'ok')

sil_rsd <- lapply(sil_list$note, function(FUNC_SIL){
  #browser()
  temp_func_data <- sil_data_check %>% select(all_of(FUNC_SIL))
  temp_func_data_mean <- temp_func_data %>% as.matrix() %>% mean()
  temp_func_data_sd <- temp_func_data %>% as.matrix() %>% sd()
  temp_func_data_rsd <- (temp_func_data_sd/temp_func_data_mean)*100
  temp_func_data_rsd
}) %>% c() %>% unlist() %>% as_tibble() %>% rename(SIL_RSD = value) %>% add_column(sil_list, .before = 1)

sil_list_warning_2 <- sil_rsd %>% filter(SIL_RSD > 30) 

dlg_message(paste( "############################################", 
                   "Warning 1! These internal standards have a very low signal and may not be present in the mixture:",
                   paste(sil_list_warning, collapse = ";     ."),
                   "############################################", 
                   "Warning 2! These internal standards have a %RSD >30%: ",
                   paste(sil_list_warning_2$note, collapse = ";     ."),
                   "############################################",
                   "double check skyline!"), 
            type = 'ok')


sil_list_warning <- c(sil_list_warning, unlist(sil_list_warning_2$sil_list))

temp_answer <- dlgInput("Do you wish to continue or use different internal standards?", "continue/change")$res

#add in a check in case the user enters the incorrect entry. It must be "continue" or "change" to continue
while(temp_answer != "continue" & temp_answer!="change"){
  temp_answer <- dlgInput("Do you wish to continue or use different internal standards?", "continue/change")$res
}

# if change has been selected the user now has an opportunity to change the internal standard import file template
if(temp_answer == "change"){
  dlg_message("OK - please edit and select a new internal standard template file now")
}
}

# this section creates a response ratio by dividing the signal area for each target lipid by the peak area from the appropriate SIL IS metabolite. As defined in the imported template above.

ratio_data <- apply(as_tibble(colnames(lipid_data)), 1, function(LIPID){
  #browser()
  func_data <- lipid_data %>% select(all_of(LIPID))
  sil_to_use <- sil_target_list$note[which(sil_target_list$precursor_name==LIPID)]
  func_data_sil <- sil_data %>% select(sil_to_use)
  normalised_data <- func_data/func_data_sil
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

