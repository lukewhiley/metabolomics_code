# intensity threshold filter to remove lipids that are not present in the dataset

individual_lipid_data_intensity <- lapply(master_lipid_data$lipid_target %>% unique(), function(FUNC_INTENSITY){
  #browser()
  sampleID <- master_lipid_data$replicate %>% unique() %>% as_tibble() # create list of sample IDs
  temp_data <- master_lipid_data %>% filter(lipid_target == FUNC_INTENSITY) %>% select(replicate, height)
  colnames(temp_data) <- c("value", FUNC_INTENSITY) 
  temp_data <- left_join(sampleID, temp_data, by = "value") 
}) %>% bind_cols() %>% select(all_of(lipid)) %>% add_column(sampleID, .before = 1) %>% filter(sampleID %in% individual_lipid_data_sil_tic_filtered$sampleID)

individual_lipid_data_intensity <- individual_lipid_data_intensity %>% select(!contains("SIL"))

lipid_intensity_list <- individual_lipid_data_intensity %>% select(contains("(")) %>% colnames()

intensity_threshold_fail_action <- "change"
while (intensity_threshold_fail_action == "change") {
  
# user choice on intensity threshold to be applied
intensity_threshold <- NA
while(is.na(intensity_threshold)) {
  intensity_threshold <- dlgInput("What signal threshold do you want to set for the intensity cut-off filter", "e.g. recommended default = 2500")$res %>% as.numeric()
}

# user choice on % of samples that must pass intensity threshold
intensity_threshold_percentage <- NA
while(is.na(intensity_threshold_percentage)) {
  intensity_threshold_percentage <- dlgInput("what % of samples do you want over this threshold?", "e.g. recommended default = 50%")$res %>% as.numeric()
}

# user choice should filter intensity be set for samples/LTR/both
intensity_threshold_ltr <- "blank"
while(intensity_threshold_ltr != "samples" & intensity_threshold_ltr != "LTR" & intensity_threshold_ltr != "both") {
  intensity_threshold_ltr <- dlgInput("Do you want to apply the filtering using samples/LTR/both. Recommended default is LTR - allows for better QC control.", "samples/LTR/both")$res
}


#apply filters defined by the user
lipid_intensity_filter_fail <- lapply(lipid_intensity_list, function(FUNC_INTENSITY){
#browser()
  #first filter data according to user input - check intensity in LTRs, samples, or both
  if(intensity_threshold_ltr == "LTR"){
    individual_lipid_data_intensity_temp <- individual_lipid_data_intensity %>% filter(grepl("LTR", sampleID))
  }
  if(intensity_threshold_ltr == "samples"){
    individual_lipid_data_intensity_temp <- individual_lipid_data_intensity %>% filter(!grepl("LTR", sampleID))
  }
  if(intensity_threshold_ltr == "both"){
    individual_lipid_data_intensity_temp <- individual_lipid_data_intensity 
  }
  
  #then find features that are greater than the user defined intensity filter cut-off
  temp_filter <- individual_lipid_data_intensity_temp %>% select(sampleID, all_of(FUNC_INTENSITY)) 
  temp_idx <- which(temp_filter[,2] > intensity_threshold)
  
  #calculate how many targets are above the threshold % set by the user 
  threshold_percentage <- (100/nrow(individual_lipid_data_intensity_temp))*length(temp_idx)
  if(threshold_percentage < intensity_threshold_percentage){
  failed_lipid <- FUNC_INTENSITY
  failed_lipid
  }
 }) %>% c() %>% unlist()

intensity_threshold_fail_action <- dlgInput(paste(length(lipid_intensity_filter_fail), " lipids failed the intensity check. Do you want to remove/keep or change the threshold settings?"), "remove/keep/change")$res
while(intensity_threshold_fail_action != "remove" & intensity_threshold_fail_action != "keep" & intensity_threshold_fail_action != "change") {
  intensity_threshold_fail_action <- dlgInput(paste(length(lipid_intensity_filter_fail), " lipids failed the intensity check. Do you want to remove/keep or change the threshold settings?"), "remove/keep/change")$res
}
}

if(intensity_threshold_fail_action == "remove"){individual_lipid_data_sil_tic_intensity_filtered <- individual_lipid_data_sil_tic_filtered %>% select(!all_of(lipid_intensity_filter_fail))}

if(intensity_threshold_fail_action == "keep"){individual_lipid_data_sil_tic_intensity_filtered <- individual_lipid_data_sil_tic_filtered}

