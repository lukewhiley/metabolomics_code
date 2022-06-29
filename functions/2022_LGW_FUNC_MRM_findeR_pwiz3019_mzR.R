# RT findeR
mzR_mrm_findR <- function(FUNC_mzR,
                      FUNC_mrm_guide,
                      FUNC_OPTION_qc_type,
                      FUNC_OPTION_max_qc_replicates){
  #browser()
  #list mzR objects
  mzML_filelist <- names(FUNC_mzR)
  #list mzR objects matching qc tyope
  mzML_filelist_qc <- mzML_filelist[grep(FUNC_OPTION_qc_type, mzML_filelist)]
  
  #if number of mzR[qc_type] > FUNC_OPTION_max_qc_replicates in the mzR filelist then take a filter subset of mzR to perform RT find
  if(length(mzML_filelist_qc) > FUNC_OPTION_max_qc_replicates) {
    mzML_filelist_idx <- c(seq(1, length(mzML_filelist_qc), 
                               by = floor(length(mzML_filelist_qc)/
                                            FUNC_OPTION_max_qc_replicates)), 
                           length(mzML_filelist_qc))
    mzML_filelist_crop <- mzML_filelist_qc[mzML_filelist_idx]
  }
  
  #if number of mzR[qc_type] < | == FUNC_OPTION_max_qc_replicates in the mzR filelist then use all of the available mzR for RT optimisation
  if(length(mzML_filelist_qc) < FUNC_OPTION_max_qc_replicates | length(mzML_filelist_qc) == FUNC_OPTION_max_qc_replicates){
    mzML_filelist_crop <- mzML_filelist_qc
  }
  
#browser()    
rt_find <- NULL
    #for each mrm transtion in the transition data
for (idx_mrm in 1:nrow(FUNC_mrm_guide)){
      # find precursor reference
      precursor_mz <- FUNC_mrm_guide$precursor_mz[idx_mrm] %>% 
        round(2) 
      #find product ion reference
      product_mz <- FUNC_mrm_guide$product_mz[idx_mrm] %>% 
        round(2) 
      
      #find transition in each mzML file and find median peak apex
      mzml_rt_apex_out <- NULL
      for(idx_mzML in mzML_filelist_crop){
        #find the data channel in the mzml file containing the data
        idx_mrm_channel <- which(
          mzR::chromatogramHeader(FUNC_mzR[[idx_mzML]])$precursorIsolationWindowTargetMZ == precursor_mz &
            mzR::chromatogramHeader(FUNC_mzR[[idx_mzML]])$productIsolationWindowTargetMZ == product_mz
        )

      #only complete the below if idx_mrm_channel finds a single unique match
      if(length(idx_mrm_channel) ==1){
      #find scan index of max intensity within mrm channel
      mzml_max_intensity <- which.max(
        mzR::chromatograms(FUNC_mzR[[idx_mzML]])[[idx_mrm_channel]][,2]
      )
      #find rt_peak_apex
      mzml_rt_apex <- mzR::chromatograms(FUNC_mzR[[idx_mzML]])[[idx_mrm_channel]]$time[mzml_max_intensity] %>% round(2)
      #c() mzml rt apex from all mzml
      mzml_rt_apex_out <- c(mzml_rt_apex_out, mzml_rt_apex)
      }
      }
      if(length(mzml_rt_apex_out) > 0){
      mzml_median_rt <- median(mzml_rt_apex_out)
      FUNC_mrm_guide$explicit_retention_time[idx_mrm] <- round(mzml_median_rt, 2)
      }
}
  #output final table
  FUNC_mrm_guide
}