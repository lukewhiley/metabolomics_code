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

sil_check_status <- "change"

while(sil_check_status == "change"){
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

#flag samples with SIL x standard deviations below mean
temp_answer <- dlgInput("What do you wish to set for the fail cut off filter.  x number of standard deviations from the mean", "e.g.   x = 2")$res

mean_sil_tic <- mean(total_summed_sil$SIL_TIC)
sd_sil_tic <- sd(total_summed_sil$SIL_TIC)
sil_sd_mean_cut_off_lower <- mean_sil_tic - (as.numeric(temp_answer)*sd_sil_tic)
sil_sd_mean_cut_off_upper <- mean_sil_tic + (as.numeric(temp_answer)*sd_sil_tic)

sil_qc_fail <- total_summed_sil$sampleID[which(total_summed_sil$SIL_TIC < sil_sd_mean_cut_off_lower | total_summed_sil$SIL_TIC > sil_sd_mean_cut_off_upper)] %>% as_tibble %>% rename(sampleID = value)
sil_qc_fail$fail_point <- "sil"

#sil_qc_fail
temp_answer <- dlgInput(paste(nrow(sil_qc_fail), "sample FAILED the SIL QC check.  Do you want to remove failed samples?"), "yes/no")$res
if(temp_answer == "yes"){individual_lipid_data_sil_filtered <- individual_lipid_data %>% filter(!sampleID %in% sil_qc_fail$sampleID)}

#visualise for reports
total_summed_sil$removed <- "pass_qc"
total_summed_sil$removed[total_summed_sil$sampleID %in% sil_qc_fail$sampleID] <- "removed"

p <- plot_ly(type = "scatter", total_summed_sil, x = ~sample_idx, y = ~LOG_SIL_TIC, text = ~sampleID, color = ~removed, colors = c("lightblue3", "red"))

plate_number <- unique(plate_id) %>% substr(14,14) %>% unique()
plate_idx <- lapply(unique(plateid), function(plateID){grep(plateID, total_summed_sil$sampleID)[1]}) %>% unlist()
for (idx_line in 2:length(plate_idx)){
  p <- add_trace(p, x = plate_idx[idx_line], type = 'scatter', mode = 'lines', color = paste("plate_", plate_number[idx_line], sep=""), line = list(color = "grey", dash = "dash"), showlegend = FALSE)
}
p <- add_trace(p, y = log(sil_sd_mean_cut_off_lower), type = 'scatter', mode = 'lines', color = "SIL QC threshold", line = list(color = "red", dash = "dash"))
p <- add_trace(p, y = log(sil_sd_mean_cut_off_upper), type = 'scatter', mode = 'lines', color = "SIL QC threshold", line = list(color = "red", dash = "dash"));p

sil_check_p <- p
if(!dir.exists(paste(project_dir, "/html_files", sep=""))){dir.create(paste(project_dir, "/html_files", sep=""))} # create a new directory to store html widgets
saveWidget(sil_check_p, file = paste(project_dir, "/html_files/",project_name, "_", user_name, "_SIL_check_plot.html", sep=""))# save plotly widget
browseURL(paste(project_dir, "/html_files/",project_name, "_", user_name, "_SIL_check_plot.html", sep="")) #open plotly widget in internet browser

sil_check_status <- dlgInput("Check the plot. Are you happy to continue? or do wish to change the exclusion threshold?", "continue/change")$res
while(sil_check_status != "continue" & sil_check_status != "change"){
  sil_check_status <- dlgInput("Error. Check the plot. Are you happy to continue? or do wish to change the exclusion threshold?", "continue/change")$res
}
}
