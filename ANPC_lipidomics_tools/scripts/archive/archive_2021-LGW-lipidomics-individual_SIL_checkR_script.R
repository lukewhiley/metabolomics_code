# check to see if all the IS are present

dlg_message("Checks to see if all internal standards are present in the SIL internal standard mix", type = 'ok')

sil_data_check <- individual_lipid_data_sil_tic_filtered %>% select(sampleID, plate_id, contains("SIL")) %>% filter(grepl("LTR", sampleID))

sil_list <- sil_data_check %>% select(contains("SIL")) %>% colnames() 

sil_sum <- lapply(sil_list, function(FUNC_SIL){
  #browser()
  temp_func_data <- sil_data_check %>% select(all_of(FUNC_SIL))
  temp_func_data_sum <- temp_func_data %>% as.matrix() %>% sum() %>% log()
}) %>% c() %>% unlist() %>% as_tibble() %>% rename(SIL_SUM = value) %>% add_column(sil_list, .before = 1) %>% arrange(SIL_SUM)

sil_sum_q1 <- quantile(sil_sum$SIL_SUM, 0.25) %>% as.numeric()
inter_quantile_range <- as.numeric(quantile(sil_sum$SIL_SUM, 0.75)) - as.numeric(quantile(sil_sum$SIL_SUM, 0.25))
sil_sum_lower_threshold <- sil_sum_q1 - inter_quantile_range

sil_list_warning <- sil_sum$sil_list[which(sil_sum$SIL_SUM < sil_sum_lower_threshold)]

dlg_message(paste("Warning! These internal standards do not look like they are present in the samples:", 
                  paste(sil_list_warning, collapse = ", "), 
                  " - double check skyline!"), 
            type = 'ok')

sil_rsd <- lapply(sil_list, function(FUNC_SIL){
  #browser()
  temp_func_data <- sil_data_check %>% select(all_of(FUNC_SIL))
  temp_func_data_mean <- temp_func_data %>% as.matrix() %>% mean()
  temp_func_data_sd <- temp_func_data %>% as.matrix() %>% sd()
  temp_func_data_rsd <- (temp_func_data_sd/temp_func_data_mean)*100
  temp_func_data_rsd
}) %>% c() %>% unlist() %>% as_tibble() %>% rename(SIL_RSD = value) %>% add_column(sil_list, .before = 1)

sil_list_warning_2 <- sil_rsd %>% filter(SIL_RSD > 30) 

dlg_message(paste("Warning! These internal standards have a %RSD >30% in LTRs", 
                  paste(sil_list_warning_2$sil_list, collapse = ", "),
                  " - double check skyline!"), 
                  type = 'ok')


sil_list_warning <- c(sil_list_warning, unlist(sil_list_warning_2$sil_list))


