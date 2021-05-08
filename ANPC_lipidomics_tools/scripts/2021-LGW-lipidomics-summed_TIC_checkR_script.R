###################################
###### summed TIC checkR #########
##################################

#this script fits into the lipidomics lipid exploreR. 
# the script will look for samples that have been got preparation errors
# the script sums the total intensity from all endogenous targets (e.g. not internal standards)
dlg_message("Summed TIC check. This next step will assess the summed TIC accross all of the samples. If samples have been incorrectly prepared the summed TIC intensity will be too low/high. Used to remove obvious outliers only.", type = 'ok')

total_summed_tic <- apply(individual_lipid_data_sil_filtered %>% select(sampleID), 1, function(summedTIC){
  #browser()
  temp_data <- individual_lipid_data_sil_filtered %>% filter(sampleID == summedTIC) %>% select(-sampleID, -plate_id) %>% select(!contains("SIL")) %>% rowSums(na.rm = TRUE)
}) %>% c() %>% as_tibble() %>%  add_column(individual_lipid_data_sil_filtered$sampleID, .before = 1) %>% 
  rename(summed_TIC = value, sampleID = "individual_lipid_data_sil_filtered$sampleID")

plateid <- str_extract(individual_lipid_data_sil_filtered$sampleID, "PLIP.*")
plateid <- substr(plateid, 0,15)
plate_id <- paste(plateid, sub(".*\\_", "", individual_lipid_data_sil_filtered$sampleID), sep="_")

total_summed_tic <- total_summed_tic %>% add_column(plate_id, .before = 2) %>% arrange(plate_id)
total_summed_tic$sample_idx <- c(1:nrow(total_summed_tic))
total_summed_tic$LOG_summed_TIC <- log(total_summed_tic$summed_TIC)

#while loop here
tic_check_status <- "change"
while(tic_check_status == "change"){

#flag samples with SIL x number of standard deviations below mean
temp_answer <- dlgInput("What do you wish to set for the fail cut off filter.  x number of median absolute deviations from the median", "e.g. recommended default x = 3")$res
while(is.na(as.numeric(temp_answer))){
  temp_answer <- dlgInput("You did not enter a numeric value.  What do you wish to set for the fail cut off filter.  x number of median absolute deviations from the median", "e.g. recommended default x = 3")$res
}

temp_answer <- as.numeric(temp_answer)
median_summed_tic <- median(total_summed_tic$summed_TIC)
mad_summed_tic <- mad(total_summed_tic$summed_TIC)
tic_cut_off_lower <- median_summed_tic - (as.numeric(temp_answer)*mad_summed_tic)

tic_qc_fail <- total_summed_tic$sampleID[which(total_summed_tic$summed_TIC < tic_cut_off_lower)] %>% as_tibble %>% rename(sampleID = value)
tic_qc_fail$fail_point <- "tic"
tic_qc_fail_ltr <- tic_qc_fail %>% filter(grepl("LTR", sampleID))# create tibble of failed LTRs
tic_qc_fail_samples <- tic_qc_fail %>% filter(!grepl("LTR", sampleID)) # create tibble of failed samples (not LTRS)

temp_answer <- "blank"
while(temp_answer != "all" & temp_answer != "none" & temp_answer != "samples" & temp_answer != "LTR"){
  temp_answer <- dlgInput(paste(nrow(tic_qc_fail), "samples FAILED the summed TIC QC check.  ",  nrow(tic_qc_fail_ltr),"were LTRs.  Do you want to remove failed samples?"), "all/none/samples/LTR")$res
  if(temp_answer == "all"){individual_lipid_data_tic_filtered <- individual_lipid_data_sil_filtered %>% filter(!sampleID %in% tic_qc_fail$sampleID)}
  if(temp_answer == "samples"){individual_lipid_data_tic_filtered <- individual_lipid_data_sil_filtered %>% filter(!sampleID %in% tic_qc_fail_samples$sampleID)}
  if(temp_answer == "LTR"){individual_lipid_data_tic_filtered <- individual_lipid_data_sil_filtered %>% filter(!sampleID %in% tic_qc_fail_ltr$sampleID)}
  if(temp_answer == "none"){individual_lipid_data_tic_filtered <- individual_lipid_data_sil_filtered}
}

#visualise for reports
total_summed_tic$removed <- "pass_qc"
total_summed_tic$removed[total_summed_tic$sampleID %in% tic_qc_fail$sampleID] <- "removed"

p <- plot_ly(type = "scatter", mode = "markers", total_summed_tic, x = ~sample_idx, y = ~LOG_summed_TIC, text = ~sampleID, color = ~removed, colors = c("lightblue3", "red"))


plate_number <- unique(plate_id) %>% substr(14,14) %>% unique()
plate_idx <- lapply(unique(plateid), function(plateID){grep(plateID, total_summed_tic$sampleID)[1]}) %>% unlist()
for (idx_line in 2:length(plate_idx)){
  p <- add_trace(p, x = plate_idx[idx_line], type = 'scatter', mode = 'lines', color = paste("plate_", plate_number[idx_line], sep=""), line = list(color = "grey", dash = "dash"), showlegend = FALSE)
}
p <- add_trace(p, y = log(tic_cut_off_lower), type = 'scatter', mode = 'lines', color = "TIC QC threshold", line = list(color = "red", dash = "dash"))
p <- add_trace(p, y = log(median_summed_tic), type = 'scatter', mode = 'lines', color = "TIC QC threshold", line = list(color = "black", dash = "dash"))


tic_check_p <- p
if(!dir.exists(paste(project_dir, "/html_files", sep=""))){dir.create(paste(project_dir, "/html_files", sep=""))} # create a new directory to store html widgets
saveWidget(tic_check_p, file = paste(project_dir, "/html_files/",project_name, "_", user_name, "_TIC_check_plot.html", sep=""))# save plotly widget
browseURL(paste(project_dir, "/html_files/",project_name, "_", user_name, "_TIC_check_plot.html", sep="")) #open plotly widget in internet browser

tic_check_status <- dlgInput("Check the plot. Are you happy to continue? or do wish to change the exclusion threshold?", "continue/change")$res
while(tic_check_status != "continue" & sil_check_status != "change"){
  tic_check_status <- dlgInput("Error. Check the plot. Are you happy to continue? or do wish to change the exclusion threshold?", "continue/change")$res
}
}

