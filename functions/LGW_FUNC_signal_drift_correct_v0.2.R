#ANPC signal batch correction

## REQUIRED PACKAGES

# -> metabom8 (github.com/tkimhofer)
# -> tidyverse
# -> RColorBrewer
# -> plotly
# -> statTarget


## REQUIRED ARGUMENTS

# -> FUNC_project_directory: directory to create output files
# -> FUNC_data: data to be batch corrected
# -> FUNC_metabolite_list: list of metabolite targets
# -> FUNC_header_sample_id: string header of column name containing sample_id info
# -> FUNC_header_batch: string header of column name containing batch ID
# -> FUNC_header_sample_type: string header of column name containing sample_type (sample or qc)
# -> FUNC_header_run_order: string header of column name containing sample run order idx
# -> FUNC_option_method: string; RF or loess
# -> FUNC_option_coCV: %RSD in qc cut off. default is 30

lgw_signal_correction <- function(FUNC_project_directory,
                                  FUNC_data,
                                  FUNC_metabolite_list,
                                  FUNC_header_sample_id,
                                  FUNC_header_batch,
                                  FUNC_header_sample_type,
                                  FUNC_header_run_order,
                                  FUNC_option_method,
                                  FUNC_option_coCV
                                  ){

require(statTarget)
  #browser()
#create directories for use later
#dir.create(paste(FUNC_project_directory))

dir.create(paste(FUNC_project_directory, "/", Sys.Date(), "_signal_correction_results", sep=""))
setwd(paste(FUNC_project_directory, "/", Sys.Date(), "_signal_correction_results", sep=""))  

#create data list 
FUNC_list <- list()
#set master data for function
FUNC_list$master_data <- FUNC_data %>% as_tibble()
#fill infinite values created at the ratio step with a small value
FUNC_list$master_data[sapply(FUNC_list$master_data, is.infinite)] <- 1e-5

# SECTION 1 ---------------------------------------------------------------
#create the required metadata file (PhenoFile) for statTarget::shiftCor

FUNC_list$PhenoFile <- list()

# build PhenoFile file template
FUNC_list$PhenoFile$template <- FUNC_list$master_data %>% 
  select(all_of(FUNC_header_sample_id)) %>%
  rename(sample = all_of(FUNC_header_sample_id)) %>%
  add_column(FUNC_list$master_data %>% 
               select(all_of(FUNC_header_sample_id))) %>%
  add_column(FUNC_list$master_data %>%
               select(all_of(FUNC_header_batch))) %>%
  add_column(FUNC_list$master_data %>%
               select(all_of(FUNC_header_sample_type))) %>%
  rename(class = all_of(FUNC_header_sample_type)) %>% 
  add_column(FUNC_list$master_data %>%
               select(all_of(FUNC_header_sample_type))) %>%
  add_column(FUNC_list$master_data %>%
               select(all_of(FUNC_header_run_order)))

#### arrange by run order
FUNC_list$PhenoFile$template <- FUNC_list$PhenoFile$template %>%
  arrange_at(FUNC_header_run_order)

#### QC placement
#stat target needs a qc at postion 1 and last position in the run. Default run order takes this into account, but some datasets this is not the case
#in this instance this section of code artificially moves the first and last qc into position. It is completed for each batch.

FUNC_list$PhenoFile$template_qc_order <- NULL
qc_idx <- NULL
for(idx_batch in FUNC_list$PhenoFile$template %>% 
    select(all_of(FUNC_header_batch)) %>%
    unique() %>%
    as.matrix() %>%
    c()
){

  #create a temp tibble batch specific
  loop_temp_data <- FUNC_list$PhenoFile$template %>%
    filter(!!as.symbol(FUNC_header_batch) == idx_batch)
  
#ensure a sample_type "qc" is "first" and "last" in the worklist order. Required for statTarget::shiftCor
loop_qc_idx <- which(loop_temp_data %>% 
                       select(all_of(FUNC_header_sample_type)) 
                      == "qc")
# browser()
#if qc is not run before the samples - artificially move first qc to run order position 1. This is required for statTarget
if(loop_qc_idx[1] > 1){
  loop_temp_data <- loop_temp_data %>%
    slice(loop_qc_idx[1],1:nrow(loop_temp_data)) %>%
    slice(-(loop_qc_idx[1]+1))
}

#create last qc
if(loop_qc_idx[length(loop_qc_idx)] < nrow(loop_temp_data)){
  loop_temp_data <- loop_temp_data %>%
    slice(1:nrow(loop_temp_data), loop_qc_idx[length(loop_qc_idx)]) %>%
    slice(-loop_qc_idx[length(loop_qc_idx)])
}

#create total qc_idx for use later
qc_idx <- c(qc_idx,
            loop_qc_idx)

FUNC_list$PhenoFile$template_qc_order <- bind_rows(FUNC_list$PhenoFile$template_qc_order,
                                                    loop_temp_data)

}

#set sample column for statTarget requires "QC" in QC rows, and sample name in sample rows
FUNC_list$PhenoFile$template_qc_order$sample[which(FUNC_list$PhenoFile$template_qc_order %>% 
                                                      select(all_of(FUNC_header_sample_type)) == "qc")] <- paste0("QC", rep(1:length(qc_idx)))
FUNC_list$PhenoFile$template_qc_order$sample[which(FUNC_list$PhenoFile$template_qc_order %>% 
                                                      select(all_of(FUNC_header_sample_type)) == "sample")] <- paste0("sample", 
                                                                                                                       rep(1:(nrow(FUNC_list$PhenoFile$template_qc_order)-length(qc_idx))))
#set NA for class column in rows that are NA
FUNC_list$PhenoFile$template_qc_order$class[which(FUNC_list$PhenoFile$template_qc_order %>% 
                                                     select(all_of(FUNC_header_sample_type)) == "qc")] <- NA

#write out for later (requirement of statTarget)
FUNC_list$PhenoFile$template_sample_id <- FUNC_list$PhenoFile$template_qc_order %>% 
  rename(sample_id = all_of(FUNC_header_sample_id),
         batch = all_of(FUNC_header_batch),
         order = all_of(FUNC_header_run_order)) %>%
  select(sample, batch, class, order, sample_id)

FUNC_list$PhenoFile$template_sample_id$order <- c(1:nrow(FUNC_list$PhenoFile$template_sample_id))
FUNC_list$PhenoFile$template_sample_id$batch <- gsub(".*?([0-9]+).*", "\\1", FUNC_list$PhenoFile$template_sample_id$batch) %>%
  as.numeric()

#final Phenofile
FUNC_list$PhenoFile$PhenoFile <- FUNC_list$PhenoFile$template_sample_id %>%
  select(-sample_id)

# write out as csv (requirement for statTarget::shiftCor)
write_csv(x = FUNC_list$PhenoFile$PhenoFile,
          file = paste(FUNC_project_directory, "/", Sys.Date(), "_signal_correction_results", "/PhenoFile.csv", sep="")
          )


# SECTION 2 - create data for statTarget::shiftCor  -----------------------------------

FUNC_list$ProfileFile <- list()

#must have samples in columns and metabolites in rows
FUNC_list$ProfileFile$template  <- FUNC_list$master_data %>%
  select(all_of(FUNC_header_sample_id),
         all_of(FUNC_metabolite_list)) %>%
  rename(sample_id = all_of(FUNC_header_sample_id))

#match run order to sil_trend_cor_meta
FUNC_list$ProfileFile$template_qc_order <- FUNC_list$PhenoFile$template_sample_id %>%
  select(sample, sample_id) %>%
  left_join(FUNC_list$ProfileFile$template, by = "sample_id") %>%
  select(-sample_id)

#transpose tibble for statTarget
FUNC_list$ProfileFile$ProfileFile <- FUNC_list$ProfileFile$template_qc_order  %>%
    gather(key = name, value = value, 2:ncol(FUNC_list$ProfileFile$template_qc_order)) %>% 
    spread(key = names(FUNC_list$ProfileFile$template_qc_order)[1],value = 'value') %>%
  select(name, FUNC_list$ProfileFile$template_qc_order$sample)

#create a metabolite list and create metabolite code
FUNC_list$ProfileFile$metabolite_list <- FUNC_list$ProfileFile$ProfileFile %>%
  select(name) %>%
  add_column(metabolite_code = paste0("M", rep(1:nrow(FUNC_list$ProfileFile$ProfileFile))))

#add metabolite code to data
FUNC_list$ProfileFile$ProfileFile$name <- FUNC_list$ProfileFile$metabolite_list$metabolite_code

# write out as csv (requirement for statTarget::shiftCor)
write_csv(x = FUNC_list$ProfileFile$ProfileFile, 
          file = paste(FUNC_project_directory, "/", Sys.Date(), "_signal_correction_results", "/ProfileFile.csv", sep="")
          )

#files from previous work
#samPeno <- "/Users/lukegraywhiley/Library/CloudStorage/OneDrive-MurdochUniversity/projects/Ryan - CABIN/lipids/data/batch_correction/test_files/sil_trend_cor_meta.csv"
#samFile <- "/Users/lukegraywhiley/Library/CloudStorage/OneDrive-MurdochUniversity/projects/Ryan - CABIN/lipids/data/batch_correction/test_files/sil_trend_cor_data.csv"

#script files
samPeno <- paste(FUNC_project_directory, "/", Sys.Date(), "_signal_correction_results", "/PhenoFile.csv", sep="")
samFile <- paste(FUNC_project_directory, "/", Sys.Date(), "_signal_correction_results",  "/ProfileFile.csv", sep="")

#example files from online
#samPeno <- "/Users/lukegraywhiley/Library/CloudStorage/OneDrive-MurdochUniversity/projects/Ryan - CABIN/lipids/data/batch_correction/test_files/Data_example/PhenoFile.csv"
#samFile <- "/Users/lukegraywhiley/Library/CloudStorage/OneDrive-MurdochUniversity/projects/Ryan - CABIN/lipids/data/batch_correction/test_files/Data_example/ProfileFile.csv"

#browser()
if(FUNC_option_method == "RF"){
  capture.output(
shiftCor(samPeno = samPeno,
         samFile =  samFile,
         Frule = 0.8,
         ntree = 500,
         MLmethod = 'QCRFSC',
         QCspan = 0,
         imputeM = "minHalf",
         plot = FALSE,
         coCV = FUNC_option_coCV
         )
)
}

if(FUNC_option_method == "loess"){
  capture.output(
  shiftCor(samPeno = samPeno,
           samFile =  samFile,
           Frule = 0.8,
           MLmethod = 'QCRLSC',
           QCspan = 0,
           imputeM = "minHalf",
           plot = FALSE,
           coCV = FUNC_option_coCV
  )
  )
}

#browser()
FUNC_list$corrected_data$data <- read_csv(paste(FUNC_project_directory, "/", Sys.Date(), "_signal_correction_results", "/statTarget/shiftCor/After_shiftCor/shift_all_cor.csv", sep=""),
                                          show_col_types = FALSE)

#create new list of meatbolie codes (statTarget might throw some away based on cut off CV values (e.g >30% variation))
FUNC_list$ProfileFile$metabolite_list_update <- FUNC_list$ProfileFile$metabolite_list %>% 
  filter(metabolite_code %in% (FUNC_list$corrected_data$data %>% 
              select(contains("M", ignore.case = FALSE)) %>%
              names()))

#recombine with sample filenames
FUNC_list$corrected_data$data_with_sample_id <- FUNC_list$PhenoFile$template_sample_id %>% 
  select(-class) %>% 
  right_join(FUNC_list$corrected_data$data, by = 'sample') %>%
  rename_with(~all_of(FUNC_header_sample_id),  "sample_id")

#add annotation data
FUNC_list$corrected_data$data_with_meta <- FUNC_list$master_data %>%
  select(-all_of(FUNC_metabolite_list)) %>%
  left_join(FUNC_list$corrected_data$data_with_sample_id,
            by = FUNC_header_sample_id) %>%
  select(contains("sample"), all_of(FUNC_list$ProfileFile$metabolite_list_update$metabolite_code))

#set column names to metabolite
FUNC_list$corrected_data$data_with_header <- FUNC_list$corrected_data$data_with_meta %>%
  setNames(c(FUNC_list$corrected_data$data_with_meta %>% 
                      select(-contains("M", ignore.case = FALSE)) %>% 
               names(),
             FUNC_list$ProfileFile$metabolite_list_update$name))

# SECTION 3 - concentration adjustment ------------------------------------
#because the StatTarget correction changes the output concentrations of the lipids, this next section re-scales the values based on the change (ratio) between pre and post corrected signal mean in the QCs

# #create empty list to store results
# FUNC_list$corrected_data$data_ratio_adjusted <- list()

#step one - get mean value for each metabolite in the QC samples - pre-single drift corrected data 
FUNC_list$corrected_data$qc_means <- FUNC_list$master_data %>%
  filter(!!as.symbol(FUNC_header_sample_type) == "qc") %>%
  #select(all_of(FUNC_metabolite_list)) %>%
  select(all_of(FUNC_list$ProfileFile$metabolite_list_update$name)) %>%
  colMeans() %>%
  as_tibble() %>%
  rename(original_mean = value) %>%
  add_column(metabolite = FUNC_list$ProfileFile$metabolite_list_update$name, 
             .before = "original_mean")

#step two - get mean value for each metabolite in the QC samples - post-single drift corrected data 
FUNC_list$corrected_data$qc_means <- FUNC_list$corrected_data$qc_means %>%  
  add_column(corrected_mean = FUNC_list$corrected_data$data_with_header %>%
               filter(!!as.symbol(FUNC_header_sample_type) == "qc") %>%
               select(all_of(FUNC_list$corrected_data$qc_means$metabolite)) %>%
               colMeans())

#step three - create ratio factor for concentration adjustment
FUNC_list$corrected_data$qc_means$correction_ratio <- FUNC_list$corrected_data$qc_means$corrected_mean/
  FUNC_list$corrected_data$qc_means$original_mean

#step 4 - adjust data concentrations
FUNC_list$corrected_data$data_qc_mean_adjusted <- FUNC_list$corrected_data$data_with_header
for(idx_metabolite in FUNC_list$ProfileFile$metabolite_list_update$name){
  FUNC_list$corrected_data$data_qc_mean_adjusted[[idx_metabolite]] <- FUNC_list$corrected_data$data_qc_mean_adjusted[[idx_metabolite]]/
    FUNC_list$corrected_data$qc_means %>% 
    filter(metabolite == idx_metabolite) %>% 
    select(correction_ratio) %>% 
    as.numeric()
}
 #browser()
output_list <- list()
output_list$scaled_data <- FUNC_list$corrected_data$data_with_header %>%
  select(-sample)
output_list$concentration_data <- FUNC_list$corrected_data$data_qc_mean_adjusted %>% 
  select(-sample)

output_list$concentration_data
}

#drop extra column "sample" that is added in the function and then export final corrected data


























####################################
# archive
############


### correlate data with non-signal drift ratio data to prevent over correction
###only keep features where correlate r >0.85

# FUNC_list$corrected_data$data_filtered <- NULL
# 
# for(idx_metabolite in FUNC_metabolite_list){
#   temp_loop_data <- cor(FUNC_list$corrected_data$ratio_corrected_data[[idx_metabolite]],
#                           FUNC_list$master_data[[idx_metabolite]])
#   
#   
#   
#   
#   
#   
#   #extract data from ratio data (before signal correction)
#   loop_data <- lipid_exploreR_data$individual_lipid_data_sil_tic_intensity_filtered_ratio %>%
#     select(sampleID, idx)
#   
#     processed_data <- signal_drift_corrected_data %>% 
#       select(sampleID, plateID, all_of(idx))
#     
#     test_data <- processed_data %>%
#       left_join(loop_data, by = "sampleID")
#      
#     corr_result <- cor(test_data[,3], test_data[,4]) %>% as_tibble(rownames = "lipid")
#     colnames(corr_result)[2] <- "r_value"
#     
#     corr_out <- rbind(corr_out, 
#                       corr_result)
#     
#   }
# 
# corr_out$lipid <- gsub(".x", "", corr_out$lipid)
# 
# 
# corr_filter_list <- corr_out %>%
#   filter(r_value > 0.75)
#   
# #select corrected data to those lipids that are correlated with > 0.85 to ratio data alone 
# signal_drift_corrected_data <- signal_drift_corrected_data %>% 
#   select(sampleID, plateID, all_of(corr_filter_list$lipid))
# 
# signal_drift_corrected_class_data <- create_lipid_class_data_summed(signal_drift_corrected_data)

