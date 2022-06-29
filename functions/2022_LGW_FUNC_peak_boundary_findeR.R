# RT findeR
peak_boundary_findR <- function(FUNC_data,
                          FUNC_OPTION_qc_type){

FUNC_filenames <- FUNC_data$replicate %>% unique()
FUNC_metabolite <- unique(FUNC_data$peptide)

#filter only LTRs or PQC
if(FUNC_OPTION_qc_type == "LTR"| FUNC_OPTION_qc_type == "PQC"){
  FUNC_data_qc <- FUNC_data %>% filter(grepl(paste0(FUNC_OPTION_qc_type), replicate))
}

#if no LTR or PQC is present use all samples
if(FUNC_OPTION_qc_type == "none"){
  FUNC_data_qc <- FUNC_data
}

FUNC_data_qc$area <- sapply(FUNC_data_qc$area, as.numeric) #ensure area column is numeric
rt_boundary_output <- lapply(FUNC_metabolite, function(FUNC_LIPID){
  #browser()
  #print(FUNC_LIPID)
  rt_boundary <-  filter(FUNC_data_qc, peptide == FUNC_LIPID) %>% filter(!is.na(area)) %>% filter(!grepl("conditioning", replicate)) %>% arrange(retention_time)
  if(nrow(rt_boundary) > 0) {
    start_time <- rt_boundary %>% select(start_time) %>% sapply(as.numeric) %>% min()
    end_time <- rt_boundary %>% select(end_time) %>% sapply(as.numeric) %>% max()
   
    rt_boundary_filelist <- FUNC_filenames %>% as_tibble() %>% rename(FileName = value)
    rt_boundary_filelist$FullPeptideName <- FUNC_LIPID
    rt_boundary_filelist$MinStartTime <- round(start_time,2)
    rt_boundary_filelist$MaxEndTime <- round(end_time, 2)
    rt_boundary_filelist
  }
}) %>% bind_rows

rt_boundary_output <- rt_boundary_output %>% filter(!is.na(MinStartTime)) %>% filter(!is.na(MaxEndTime)) #remove any NAs that are in the peak boundary guide

rt_boundary_output
}