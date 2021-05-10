###################################
###### summed SIL checkR #########
##################################

#this script fits into the lipidomics lipid exploreR. 
# the script will look for samples that have been got preparation errors
# the script sums the total intensity from SIL internal standard
# Samples are removed if the summed intensity is either > or < x standard deviations (SD) above or below the mean summed SIL internal standard signal of the dataset. 
# users define the SD cuttoff (x)
# Sample SIL intensity > x SD indicates high concentration of SIL, so likely as a result of excess SIL IS added
# Sample SIL intensity < x SD below mean indicates too little volume added of SIL internal standard added
# users can select if failed samples are removed

dlg_message("Internal standard check. This next step will assess the internal standards accross all of the samples. If internal standards have been incorrectly added the summed signal intensity will be too low/high.", type = 'ok')

total_summed_sil <- apply(individual_lipid_data %>% select(sampleID), 1, function(summedSIL){
  temp_data <- individual_lipid_data %>% filter(sampleID == summedSIL) %>% select(-sampleID) %>% select(contains("SIL")) %>% rowSums(na.rm = TRUE)
}) %>% c() %>% as_tibble() %>%  add_column(individual_lipid_data$sampleID, .before = 1) %>% 
  rename(SIL_TIC = value, sampleID = "individual_lipid_data$sampleID")

plateid <- str_extract(individual_lipid_data$sampleID, "PLIP.*")
plateid <- substr(plateid, 0,15)
plate_id <- paste(plateid, sub(".*\\_", "", individual_lipid_data$sampleID), sep="_")

total_summed_sil <- total_summed_sil %>% add_column(plate_id, .before = 2) %>% arrange(plate_id)
total_summed_sil$sample_idx <- c(1:nrow(total_summed_sil))
total_summed_sil$LOG_SIL_TIC <- log(total_summed_sil$SIL_TIC)

# while loop here
sil_check_status <- "change"
while(sil_check_status == "change"){

#flag samples with SIL x standard deviations below mean
temp_answer <- dlgInput("What do you wish to set for the fail cut off filter.  x number of median absolute deviations from the median", "e.g. recommended default x = 3")$res
while(is.na(as.numeric(temp_answer))){
  temp_answer <- dlgInput("You did not enter a numeric value.  What do you wish to set for the fail cut off filter.  x number of median absolute deviations from the median", "e.g. recommended default x = 3")$res
}

median_sil_tic <- median(total_summed_sil$SIL_TIC)
mad_sil_tic <- mad(total_summed_sil$SIL_TIC)
sil_cut_off_lower <- median_sil_tic - (as.numeric(temp_answer)*mad_sil_tic)
sil_cut_off_upper <- median_sil_tic + (as.numeric(temp_answer)*mad_sil_tic)

#create lists of which samples have failed the SIL internalk standard check
sil_qc_fail <- total_summed_sil$sampleID[which(total_summed_sil$SIL_TIC < sil_cut_off_lower | total_summed_sil$SIL_TIC > sil_cut_off_upper)] %>% as_tibble %>% rename(sampleID = value)
sil_qc_fail$fail_point <- "sil"
sil_qc_fail_ltr <- sil_qc_fail %>% filter(grepl("LTR", sampleID))
sil_qc_fail_samples <- sil_qc_fail %>% filter(!grepl("LTR", sampleID))

#visualise for reports
total_summed_sil$removed <- "pass_qc"
total_summed_sil$removed[total_summed_sil$sampleID %in% sil_qc_fail$sampleID] <- "removed"

p <- plot_ly(type = "scatter", mode = "markers", total_summed_sil, x = ~sample_idx, y = ~LOG_SIL_TIC, text = ~sampleID, color = ~removed, colors = c("lightblue3", "red"))

plate_number <- unique(plate_id) %>% substr(14,14) %>% unique()
plate_idx <- lapply(unique(plateid), function(plateID){grep(plateID, total_summed_sil$sampleID)[1]}) %>% unlist()
for (idx_line in 2:length(plate_idx)){
  p <- add_trace(p, x = plate_idx[idx_line], type = 'scatter', mode = 'lines', color = paste("plate_", plate_number[idx_line], sep=""), line = list(color = "grey", dash = "dash"), showlegend = FALSE)
}
p <- add_trace(p, y = log(sil_cut_off_lower), type = 'scatter', mode = 'lines', color = "SIL QC threshold", line = list(color = "red", dash = "dash"))
p <- add_trace(p, y = log(median_sil_tic), type = 'scatter', mode = 'lines', color = "SIL QC threshold", line = list(color = "black", dash = "dash"))
p <- add_trace(p, y = log(sil_cut_off_upper), type = 'scatter', mode = 'lines', color = "SIL QC threshold", line = list(color = "red", dash = "dash"))

#create html widget and display it in the users internet browser
sil_check_p <- p
if(!dir.exists(paste(project_dir, "/html_files", sep=""))){dir.create(paste(project_dir, "/html_files", sep=""))} # create a new directory to store html widgets
saveWidget(sil_check_p, file = paste(project_dir, "/html_files/",project_name, "_", user_name, "_SIL_check_plot.html", sep=""))# save plotly widget
browseURL(paste(project_dir, "/html_files/",project_name, "_", user_name, "_SIL_check_plot.html", sep="")) #open plotly widget in internet browser

#sil_qc_fail - ask the user if they wish to continue or change the threshold
sil_check_status <- dlgInput(paste(nrow(sil_qc_fail), "samples FAILED the SIL QC check.  Check the plot.  Are you happy to continue? or change the exclusion threshold?"), "continue/change")$res
while(sil_check_status != "continue" & sil_check_status != "change"){
  sil_check_status <- dlgInput(paste(nrow(sil_qc_fail), "samples FAILED the SIL QC check.  Check the plot.  Are you happy to continue? or change the exclusion threshold?"), "continue/change")$res
}
}

#sil_qc_fail - ask the user if they wish to remove all/none/samples/LTR which failed the QC check
temp_answer <- "blank"
while(temp_answer != "all" & temp_answer != "none" & temp_answer != "samples" & temp_answer != "LTR"){
  temp_answer <- dlgInput(paste(nrow(sil_qc_fail), "samples FAILED the SIL QC check.  ",  nrow(sil_qc_fail_ltr),"were LTRs.  Do you want to remove failed samples?"), "all/none/samples/LTR")$res
  if(temp_answer == "all"){individual_lipid_data_sil_filtered <- individual_lipid_data %>% filter(!sampleID %in% sil_qc_fail$sampleID)}
  if(temp_answer == "samples"){individual_lipid_data_sil_filtered <- individual_lipid_data %>% filter(!sampleID %in% sil_qc_fail_samples$sampleID)}
  if(temp_answer == "LTR"){individual_lipid_data_sil_filtered <- individual_lipid_data %>% filter(!sampleID %in% sil_qc_fail_ltr$sampleID)}
  if(temp_answer == "none"){individual_lipid_data_sil_filtered <- individual_lipid_data}
}

#tidy up environment

remove_list <- c("idx_line", "mad_sil_tic", "median_sil_tic", "plate_number", "sil_check_status", "sil_cut_off_lower", "sil_cut_off_upper", "summed_SIL_checkR_script", "remove_list",
                 "p", "total_summed_sil")
rm(list = remove_list)
      

