################################################
###### Internal standard normalization #########
#################################################

# this section normalizes each SRM lipid target with the appropriate internal standard
# Requires a template guide with internal standard transition for each target lipid SRM transition

dlg_message("Time for normalization using the internal standards :-)", type = 'ok')
dlg_message("REQUIRES - a template csv listing each lipid target and the assigned internal standard", type = 'ok')

#import transition report 3
filtered_data <- individual_lipid_data_tic_filtered %>% filter(!grepl("conditioning", sampleID))

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
paste("number of feature ratios with with an LTR RSD of <30% =", length(which(ltr_rsd$RSD < 30)))
paste("number of feature ratios with with an LTR RSD of <20% =", length(which(ltr_rsd$RSD < 20)))
paste("number of feature ratios with with an LTR RSD of <15% =", length(which(ltr_rsd$RSD < 15)))
paste("number of feature ratios with with an LTR RSD of <10% =", length(which(ltr_rsd$RSD < 10)))

ltr_rsd <- ltr_rsd %>% arrange(RSD)
plot(ltr_rsd$RSD)

lipid_keep_list <- ltr_rsd %>% filter(RSD < 30)

final_dataset <- ratio_data %>% select(sampleID, plateID, all_of(lipid_keep_list$lipid))


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

p