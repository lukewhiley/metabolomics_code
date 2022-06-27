# RT findeR
function(FUNC_data_path,
         FUNC_mrm_guide,
         FUNC_OPTION_qc_type,
         FUNC_OPTION_max_qc_replicates){
  require(MSnbase)
  #list mzML files
  mzML_filelist <- list.files(FUNC_data_path, pattern = ".mzML") %>% as_tibble() %>% filter(grepl(paste0(FUNC_qc_type), value)) %>% filter(!grepl("conditioning", value)) %>% filter(!grepl("blank", value))
  
  #if number of mzML qcs > FUNC_OPTION_max_qc_replicates samples in the mzML filelist then take a quartile sample of total LTRs
  if(length(mzML_filelist$value) > FUNC_OPTION_max_qc_replicates) {
    mzML_filelist_idx <- c(seq(1, nrow(mzML_filelist), by = floor(nrow(mzML_filelist)/FUNC_OPTION_max_qc_replicates)), nrow(mzML_filelist))
    mzML_filelist_crop <- mzML_filelist[mzML_filelist_idx,]
  }
  
  #if number of mzML qcs < | == FUNC_OPTION_max_qc_replicates samples in the mzML filelist then use all of the available samples for retention time optimisation
  if(length(mzML_filelist$value) < FUNC_OPTION_max_qc_replicates | length(mzML_filelist$value) == FUNC_OPTION_max_qc_replicates){
    mzML_filelist_crop <- mzML_filelist
  }
  
  # find RT
  for(mzML_idx in mzML_filelist_crop$value){
    test_spectra <- MSnbase::readSRMData(paste0(FUNC_data_path, "/", mzML_idx)) #read in mzML file [mzML_idx]
    
    rt_find <- NULL
    #for each mrm transtion in the transition data
    for (transition_idx in 1:nrow(test_spectra)){
      # find precursor reference
      precursor_mz <- test_spectra[transition_idx,]@precursorMz[1] %>% 
        round(2) 
      #find product ion reference
      product_mz <- test_spectra[transition_idx,]@productMz[1] %>% 
        round(2) 
      # find median max intensity scan from all of the QCs
      max_intensity_idx <- which(test_spectra[transition_idx,]@intensity == max(test_spectra[transition_idx,]@intensity)) %>% 
        median()
      #find # of scans in window
      scans_in_window <- length(test_spectra[transition_idx,]@rtime) 
      #find apex scan of RT
      apex_rt <- test_spectra[transition_idx,]@rtime[max_intensity_idx] 
      # find lipid name from tempalte
      precursor_name <- FUNC_mrm_guide$precursor_name[which(round(FUNC_mrm_guide$precursor_mz,2) == precursor_mz & round(FUNC_mrm_guide$product_mz,2) == product_mz)] %>% 
        unique() 
      
      if(length(precursor_name) == 1){
        FUNC_mrm_guide$explicit_retention_time_2 <- FUNC_mrm_guide$explicit_retention_time
        FUNC_mrm_guide$explicit_retention_time_2[which(FUNC_mrm_guide$precursor_name == precursor_name)] <- round(apex_rt, 2)
      }
    
      detach("package:MSnbase", unload = TRUE)
      
    }
  }
  FUNC_mrm_guide 
}
    
    
    
    
  #   if(mzML_idx == 1){rt_find_master <-  rt_find %>% as_tibble() %>% setNames(c("lipid", paste0(FUNC_qc_type,"_", mzML_idx)))}
  #   if(mzML_idx > 1){rt_find_master <-  rt_find %>% as_tibble() %>% setNames(c("lipid", paste0(FUNC_qc_type, "_", mzML_idx))) %>% left_join(rt_find_master, by = "lipid")}
  #   
  #   #rt_find_master <-  rt_find %>% as_tibble() %>% setNames(c("lipid", paste0("LTR_", mzML_idx))) 
  # }
  # 
  # rt_find_master[,grepl(paste0(qc_type), colnames(rt_find_master))] <- sapply(rt_find_master[,grepl(paste0(qc_type), colnames(rt_find_master))], as.numeric)
  # rt_find_master$median_rt <- NA
  # for(lipid_idx in 1:nrow(rt_find_master)){
  #   rt_find_master$median_rt[lipid_idx] <- rt_find_master[lipid_idx, grepl(paste0(qc_type), colnames(rt_find_master))] %>% as.matrix() %>% median(na.rm = TRUE)
  #   transition_metadata$explicit_retention_time[which(transition_metadata$precursor_name == rt_find_master$lipid[lipid_idx])] <- rt_find_master$median_rt[lipid_idx]
  # }
  
