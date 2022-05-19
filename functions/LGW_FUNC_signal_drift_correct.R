#ANPC signal batch correction

## REQUIRED PACKAGES

# -> metabom8 (github.com/tkimhofer)
# -> tidyverse
# -> RColorBrewer
# -> plotly
# -> statTarget


## REQUIRED ARGUMENTS

# -> FUNC_data = a tibble or data from containing data
# -> FUNC_opls_y = column name for column containing class data as y in OPLS-DA
# -> FUNC_metabolite_list = array of metabolites to use - must match appropiate column names
# -> FUNC_colour_by = column name for column containing character string or factor to colour OPLS-DA 
# -> FUNC_plot_label = column name for column containing character string or factor to label OPLS-DA plotly
# -> FUNC_scaling = scaling argument for metabom8 - only use UV or Pareto
# -> FUNC_title = title for OPLS-DA plot
# -> FUNC_project_colours = array of colours - must match length of unique number of groups

# -> FUNC_data_predict = use if wanting to predict data to model - a tibble or data from containing data to predict. Set as FALSE if no predicition required


lgw_signal_correction <- function(FUNC_project_directory,
                                  FUNC_data, 
                                  
                                  
){

require(statTarget)

dir.create(paste(project_dir, "/", Sys.Date(), "_signal_correction_results", sep=""))
setwd(paste(project_dir, "/", Sys.Date(), "_signal_correction_results", sep=""))

#fill infinite values created at the ratio step with a small value
sil_trend <- lipid_exploreR_data$individual_lipid_data_sil_tic_intensity_filtered_ratio
sil_trend[sapply(sil_trend, is.infinite)] <- 1e-5


#this section creates the required metadata file for statTarget::shiftCor

sil_trend_cor_meta <- sil_trend %>% select(sampleID) 
sil_trend_cor_meta$batch <- 1
sil_trend_cor_meta$class <- 1
sil_trend_cor_meta$class[grep("LTR", sil_trend_cor_meta$sampleID)] <- NA
sil_trend_cor_meta$order <- c(1:nrow(sil_trend_cor_meta))
sil_trend_cor_meta <- as_tibble(sil_trend_cor_meta)

#ensure an LTR is "first" and "last" in the worklist order. Required for statTarget::shiftCor
LTR_locations <- grep("LTR", sil_trend_cor_meta$sampleID)


#create first LTR
if(LTR_locations[1] > 1){
sil_trend_cor_meta$order[LTR_locations[1]] <- 1
sil_trend_cor_meta$order[1:(LTR_locations[1]-1)] <- sil_trend_cor_meta$order[1:(LTR_locations[1]-1)] +1 # re-label the samples
}

#create last LTR
if(LTR_locations[length(LTR_locations)] < nrow(sil_trend_cor_meta)){
sil_trend_cor_meta$order[LTR_locations[length(LTR_locations)]] <- nrow(sil_trend_cor_meta)
sil_trend_cor_meta$order[(LTR_locations[length(LTR_locations)]+1):nrow(sil_trend_cor_meta)] <- sil_trend_cor_meta$order[(LTR_locations[length(LTR_locations)]+1):nrow(sil_trend_cor_meta)]-1
}


#create new labels simply QCx and samples
sil_trend_cor_meta$sample <- NA
#QC
sil_trend_cor_meta$sample[is.na(sil_trend_cor_meta$class)] <- paste("QC", 1:length(sil_trend_cor_meta$sample[is.na(sil_trend_cor_meta$class)]), sep="")
#Sample
sil_trend_cor_meta$sample[!is.na(sil_trend_cor_meta$class)] <- paste("sample", 1:length(sil_trend_cor_meta$sample[!is.na(sil_trend_cor_meta$class)]), sep="")

#order
sil_trend_cor_meta <- sil_trend_cor_meta %>% arrange(`order`)
sil_trend_cor_meta_2 <- sil_trend_cor_meta %>% select(sample, batch, class, order)

# write out as csv (requirement for statTarget::shiftCor)
write_csv(x = sil_trend_cor_meta_2,
          file = paste(project_dir, "/", Sys.Date(), "_signal_correction_results", "/sil_trend_cor_meta.csv", sep="")
          )

#create data for statTarget::shiftCor
sil_trend_cor_data <- sil_trend
colnames(sil_trend_cor_data) <- gsub("[():]", "_", colnames(sil_trend_cor_data))
sil_trend_cor_data <- sil_trend_cor_meta %>% select(sampleID, sample) %>% right_join(sil_trend_cor_data, by = 'sampleID') %>% select(sample, contains("_"))
sil_trend_cor_data_2 <- as_tibble(cbind(nms = names(sil_trend_cor_data), t(sil_trend_cor_data)))
colnames(sil_trend_cor_data_2) <- sil_trend_cor_data_2[1,]
sil_trend_cor_data_2 <- sil_trend_cor_data_2[-1,]
sil_trend_cor_data_2[,2:ncol(sil_trend_cor_data_2)] <- sapply(sil_trend_cor_data_2[,2:ncol(sil_trend_cor_data_2)], as.numeric)
sil_trend_cor_data_2 <- sil_trend_cor_data_2 %>% rename(name=sample)

write_csv(x = sil_trend_cor_data_2, 
          file = paste(project_dir, "/", Sys.Date(), "_signal_correction_results",  "/sil_trend_cor_data.csv", sep="")
          )


samPeno <- paste(project_dir, "/", Sys.Date(), "_signal_correction_results", "/sil_trend_cor_meta.csv", sep="")
samFile <- paste(project_dir, "/", Sys.Date(), "_signal_correction_results",  "/sil_trend_cor_data.csv", sep="")

signal_drift_method <- "blank"
if(workflow_choice == "default"){
  signal_drift_method <- "RF"
}
while(signal_drift_method != "loess"& signal_drift_method != "RF"){
  signal_drift_method <- dlgInput("What method do you wish to use? Loess or random forrest?", "loess/RF")$res
}

if(signal_drift_method == "RF"){
shiftCor(samPeno = samPeno,
         samFile =  samFile,
         Frule = 0.8,
         ntree = 500,
         MLmethod = 'QCRFSC',
         QCspan = 0,
         imputeM = "KNN",
         plot = FALSE,
         coCV = 1000
         )
}

if(signal_drift_method == "loess"){
  shiftCor(samPeno = samPeno,
           samFile =  samFile,
           Frule = 0.8,
           MLmethod = 'QCRLSC',
           QCspan = 0,
           imputeM = "KNN",
           plot = FALSE,
           coCV = 1000
  )
}

corrected_data <- read_csv(paste(project_dir, "/", Sys.Date(), "_signal_correction_results", "/statTarget/shiftCor/After_shiftCor/shift_all_cor.csv", sep=""))
corrected_data <- sil_trend_cor_meta %>% select(sampleID, sample) %>% right_join(corrected_data, by = 'sample') %>% select(-sample, -class) %>% arrange(sampleID)

#get column names back to correct format (replace _ with (:) e.g. CE_14_0_ to CE(14:0)
#first add bracket after lipid class
for (idx_name in lipid_class_list$value){
  colnames(corrected_data)[grep(paste0(idx_name, "_"), colnames(corrected_data))] <- gsub(idx_name, paste0(idx_name, "("), colnames(corrected_data)[grep(paste0(idx_name, "_"), colnames(corrected_data))])
}

colnames(corrected_data) <- paste0(colnames(corrected_data), ")")
colnames(corrected_data) <- sub("\\(_", "\\(", colnames(corrected_data))
colnames(corrected_data) <- sub("\\_)", "\\)", colnames(corrected_data))
colnames(corrected_data) <- sub("_",":",colnames(corrected_data))
colnames(corrected_data) <- sub("_","/",colnames(corrected_data))
colnames(corrected_data) <- sub("_",":",colnames(corrected_data))
colnames(corrected_data) <- sub("/","_",colnames(corrected_data))
colnames(corrected_data)[which(colnames(corrected_data) == "sampleID)")] <- "sampleID"

corrected_data <- lipid_exploreR_data$individual_lipid_data_sil_tic_intensity_filtered_ratio %>% 
  select(sampleID, plateID) %>% right_join(corrected_data, by = "sampleID")


corrected_lipid_list <- corrected_data %>% select(contains("(")) %>% colnames()

#because the correction changes the concentrations of the lipids, this next section re-scales the values based on the change (ratio) between pre and post corrected signal mean in the LTR QCs
signal_drift_corrected_data <- lapply(corrected_lipid_list, function(FUNC_LIPID_NORM){
  #browser()
  
  corrected_data_mean <- corrected_data %>% 
    filter(grepl("LTR", sampleID)) %>% 
    select(all_of(FUNC_LIPID_NORM)) %>% 
    as.matrix() %>% 
    mean()
  
  pre_corrected_data_mean <- lipid_exploreR_data$individual_lipid_data_sil_tic_intensity_filtered_ratio %>% 
    as_tibble() %>% 
    filter(grepl("LTR", sampleID)) %>% 
    select(all_of(FUNC_LIPID_NORM)) %>% 
    as.matrix() %>% 
    mean()
  
  normalization_ratio <- corrected_data_mean/pre_corrected_data_mean
  corrected_data_norm <- corrected_data %>% select(FUNC_LIPID_NORM)/normalization_ratio
  corrected_data_norm
}) %>% bind_cols %>% as_tibble()

signal_drift_corrected_data <- signal_drift_corrected_data %>% add_column(select(corrected_data, sampleID, plateID), .before = 1)
#

### correlate data with non-signal drift ratio data to prevent over correction
###only keep features where correlate r >0.85

corr_out = NULL

for(idx in signal_drift_corrected_data %>% select(contains("(")) %>% names()){
  
  #extract data from ratio data (before signal correction)
  loop_data <- lipid_exploreR_data$individual_lipid_data_sil_tic_intensity_filtered_ratio %>%
    select(sampleID, idx)
  
    processed_data <- signal_drift_corrected_data %>% 
      select(sampleID, plateID, all_of(idx))
    
    test_data <- processed_data %>%
      left_join(loop_data, by = "sampleID")
     
    corr_result <- cor(test_data[,3], test_data[,4]) %>% as_tibble(rownames = "lipid")
    colnames(corr_result)[2] <- "r_value"
    
    corr_out <- rbind(corr_out, 
                      corr_result)
    
  }

corr_out$lipid <- gsub(".x", "", corr_out$lipid)


corr_filter_list <- corr_out %>%
  filter(r_value > 0.75)
  
#select corrected data to those lipids that are correlated with > 0.85 to ratio data alone 
signal_drift_corrected_data <- signal_drift_corrected_data %>% 
  select(sampleID, plateID, all_of(corr_filter_list$lipid))

signal_drift_corrected_class_data <- create_lipid_class_data_summed(signal_drift_corrected_data)

}
