# RT findeR
mrm_findR <- function(FUNC_mzml_data_path,
         FUNC_mrm_guide,
         FUNC_OPTION_qc_type,
         FUNC_OPTION_max_qc_replicates){
  force(FUNC_mzml_data_path)
  force(FUNC_mrm_guide)
  force(FUNC_OPTION_qc_type)
  force(FUNC_OPTION_max_qc_replicates)
  #browser()
  #list mzML files
  mzML_filelist <- list.files(FUNC_mzml_data_path, pattern = ".mzML") %>% 
    as_tibble() %>% 
    filter(grepl(paste0(FUNC_OPTION_qc_type), value)) %>% 
    filter(!grepl("conditioning", value)) %>% 
    filter(!grepl("blank", value))
  #if number of mzML qcs > FUNC_OPTION_max_qc_replicates samples in the mzML filelist then take a quartile sample of total LTRs
  if(length(mzML_filelist$value) > FUNC_OPTION_max_qc_replicates) {
    mzML_filelist_idx <- c(seq(1, nrow(mzML_filelist), 
                               by = floor(nrow(mzML_filelist)/
                                            FUNC_OPTION_max_qc_replicates)), 
                           nrow(mzML_filelist))
    mzML_filelist_crop <- mzML_filelist[mzML_filelist_idx,]
  }
  
  #if number of mzML qcs < | == FUNC_OPTION_max_qc_replicates samples in the mzML filelist then use all of the available samples for retention time optimisation
  if(length(mzML_filelist$value) < FUNC_OPTION_max_qc_replicates | length(mzML_filelist$value) == FUNC_OPTION_max_qc_replicates){
    mzML_filelist_crop <- mzML_filelist
  }
  
  #read in mzml files into a list
  FUNC_spectra <- list()
  #browser()
  # #due to a package conflict with MSnBase detach and reattach packages
  # detachAllPackages <- function() {
  #   basic.packages.blank <- c("stats", "graphics", "grDevices", "utils", "datasets", "methods", "base")    
  #   basic.packages <- paste("package:", basic.packages.blank, sep = "")   
  #   package.list <- search()[ifelse(unlist(gregexpr("package:", search())) == 1, TRUE, FALSE)]   
  #   package.list <- setdiff(package.list, basic.packages)   
  #   if (length(package.list) > 0) {   
  #     for (package in package.list) {   
  #       detach(package, character.only = TRUE)   
  #     }   
  #   } 
  #   package.list
  # } 
  # 
 # # browser()
 #  package_list <- detachAllPackages()
 #  #browser()
 #  library(MSnbase)
  
  for(idx_mzML in mzML_filelist_crop$value){
    #read in mzML file [mzML_idx])
    FUNC_spectra[[idx_mzML]] <- q3ML::openFile(paste0(FUNC_mzml_data_path, "/", idx_mzML)) #read in mzML file [mzML_idx]
  }
  detach("package:MSnbase", unload=TRUE)
  #browser()
  for(idx_package in package_list){
library(package = sub("package:", "", idx_package), character.only = TRUE)
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
      for(idx_mzml in names(FUNC_spectra)){
        #find the data channel in the mzml file containing the data
      idx_mrm_channel <- which(
        FUNC_spectra[[idx_mzml]]@featureData@data$precursorIsolationWindowTargetMZ == precursor_mz &
          FUNC_spectra[[idx_mzml]]@featureData@data$productIsolationWindowTargetMZ == product_mz)
      #only complete the below if idx_mrm_channel finds a single unique match
      if(length(idx_mrm_channel) ==1){
      #find scan index of max intensity within mrm channel
      mzml_max_intensity <- which.max(FUNC_spectra[[idx_mzml]][idx_mrm_channel]@intensity)
      #find rt_peak_apex
      mzml_rt_apex <- FUNC_spectra[[idx_mzml]][idx_mrm_channel]@rtime[mzml_max_intensity] %>% round(2)
      #c() mzml rt apex from all mzml
      mzml_rt_apex_out <- c(mzml_rt_apex_out, mzml_rt_apex)
      mzml_median_rt <- median(mzml_rt_apex_out)
      FUNC_mrm_guide$explicit_retention_time[idx_mrm] <- round(mzml_median_rt, 2)
      }
      }
}
  #output finalk table
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
  
