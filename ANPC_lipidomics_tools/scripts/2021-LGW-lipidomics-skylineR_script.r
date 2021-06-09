package_list <- c("plyr", "tidyverse", "janitor", "gridExtra", "ggpubr", "readxl", "cowplot", "scales", "devtools", "metabom8", "shiny", "plotly", "svDialogs", "DataEditR", "htmlwidgets", "httr")
loaded_packages <- lapply(package_list, require, character.only = TRUE)
rm(loaded_packages)

dlg_message("Welcome to skylineR! :-)", type = 'ok')

dlg_message("Please select your project folder", type = 'ok'); project_dir <- rstudioapi::selectDirectory(); setwd(project_dir) # save project directory root location and temporarily set the working directory

#pop-up user input here for project name and user initials
if(exists("project_name") != TRUE){project_name <- dlgInput("what is the name of the project?", "example_project")$res}
if(exists("user_name") != TRUE){user_name <- dlgInput("Insert your initials", "example_initials")$res}

# read in lipid MS target transition information file
dlg_message("please select your lipid MS target transition information file", type = 'ok'); transition_metadata <- read_csv(file.choose(.))
transition_metadata_headers <- colnames(transition_metadata)
transition_metadata <- clean_names(transition_metadata)
metabolite_target_list <- transition_metadata %>% select(precursor_name)

library(MSnbase)
dlg_message("Select the folder containing the mzML files", type = 'ok'); mzML_directory <- rstudioapi::selectDirectory()
#dlg_message("Select the mzML file created from the FIRST LTR from the run to read into R", type = 'ok'); test_spectra_1 <- MSnbase::readSRMData(file.choose(.))
#dlg_message("Select the mzML file created from the LAST LTR from the run to read into R", type = 'ok'); test_spectra_2 <- MSnbase::readSRMData(file.choose(.))

mzML_filelist <- list.files(mzML_directory, pattern = ".mzML") %>% as_tibble() %>% filter(grepl("LTR", value)) %>% filter(!grepl("conditioning", value))
mzML_filelist_idx <- c(seq(1, nrow(mzML_filelist), by = floor(nrow(mzML_filelist)/4)), nrow(mzML_filelist))
mzML_filelist_crop <- mzML_filelist[mzML_filelist_idx,]

pb <- progress_bar$new(total = 100)
for(mzML_idx in 1:nrow(mzML_filelist_crop)){
  pb$tick()
  test_spectra <- MSnbase::readSRMData(paste0(mzML_directory, "/", mzML_filelist_crop$value[mzML_idx]))
  
  rt_find <- NULL
  for (transition_idx in 1:nrow(test_spectra)){
    precursor_mz <- test_spectra[transition_idx,]@precursorMz[1]
    product_mz <- test_spectra[transition_idx,]@productMz[1]
    #if(precursor_mz == 614.6){print(transition_idx)}
    max_intensity_idx <- which(test_spectra[transition_idx,]@intensity == max(test_spectra[transition_idx,]@intensity)) %>% median()
    #max_intensity_idx_2 <- which(test_spectra_2[transition_idx,]@intensity == max(test_spectra_2[transition_idx,]@intensity)) %>% median()
    scans_in_window <- length(test_spectra[transition_idx,]@rtime)
    #scans_in_window_2 <- length(test_spectra_1[transition_idx,]@rtime)
    apex_rt <- test_spectra[transition_idx,]@rtime[max_intensity_idx]
    precursor_name <- results_1$precursor_name[which(results_1$precursor_mz == precursor_mz & results_1$product_mz == product_mz)] %>% unique()
    rt_find <- rbind(rt_find, c(precursor_name, round(apex_rt,2)))
  }
  
  if(mzML_idx == 1){rt_find_master <-  rt_find %>% as_tibble() %>% setNames(c("lipid", paste0("LTR_", mzML_idx)))}
  if(mzML_idx > 1){rt_find_master <-  rt_find %>% as_tibble() %>% setNames(c("lipid", paste0("LTR_", mzML_idx))) %>% left_join(rt_find_master, by = "lipid")}
  
  #rt_find_master <-  rt_find %>% as_tibble() %>% setNames(c("lipid", paste0("LTR_", mzML_idx))) 
}
    
rt_find_master[,grepl("LTR", colnames(rt_find_master))] <- sapply(rt_find_master[,grepl("LTR", colnames(rt_find_master))], as.numeric)
rt_find_master$median_rt <- NA
for(lipid_idx in 1:nrow(rt_find_master)){
  rt_find_master$median_rt[lipid_idx] <- rt_find_master[lipid_idx, grepl("LTR", colnames(rt_find_master))] %>% as.matrix() %>% median(na.rm = TRUE)
  transition_metadata$explicit_retention_time[which(transition_metadata$precursor_name == rt_find_master$lipid[lipid_idx])] <- rt_find_master$median_rt[lipid_idx]
}

  # if(max_intensity_idx_1 < 5){
  #   max_intensity_idx_1 <- max_intensity_idx_1 + floor(scans_in_window_1/100*10)
  # }
  # if(max_intensity_idx_1 > (scans_in_window_1-5)){
  #   max_intensity_idx_1 <- max_intensity_idx_1 - floor(scans_in_window_1/100*10)
  # }
  # 
  # if(max_intensity_idx_2 < 5){
  #   max_intensity_idx_2 <- max_intensity_idx_2 + floor(scans_in_window_1/100*10)
  # }
  # if(max_intensity_idx_2 > (scans_in_window_2-5)){
  #   max_intensity_idx_2 <- max_intensity_idx_2 - floor(scans_in_window_1/100*10)
  # }
  
  
 
  #apex_rt_2 <- test_spectra_2[transition_idx,]@rtime[max_intensity_idx_2]
  #apex_rt <- (apex_rt_1+apex_rt_2)/2
  

detach("package:MSnbase", unload = TRUE)

# new_rt <- lapply(metabolite_target_list$precursor_name, function(idx_m){
#   (filter(results_1, precursor_name == idx_m) %>% filter(!is.na(area)) %>% arrange(desc(area)))[1,] %>% select(precursor_name, retention_time)
# }) %>% bind_rows
#transition_metadata$explicit_retention_time <- new_rt$retention_time #update transition list with new RT

transition_metadata_2 <- transition_metadata %>% filter(!is.na(explicit_retention_time)) # drop features without a new RT 
metabolite_target_list <- transition_metadata_2 %>% select(precursor_name)
colnames(transition_metadata_2) <- c(transition_metadata_headers)

write_csv(transition_metadata_2, (file = paste(project_dir, "/", user_name, "_", project_name, "_skylineR_RT_update.csv", sep="")))


#create list of mzML files
#dlg_message("Please select the folder that contains your mzML files", type = 'ok')
#mzml_file_location <- rstudioapi::selectDirectory() # save project directory root location
#filenames <- list.files(pattern = ".mzML", paste(mzml_file_location)) %>% as_tibble() %>% rename(FileName = value)

dlg_message("1. Please open skylineMS software", type = 'ok')
dlg_message("2. Create new small molecule file", type = 'ok')
dlg_message("3. Import the project transition list (in the project folder - 'skylineR_RT_update') from csv file by navigating to File -> import -> transition list", type = 'ok')
dlg_message("4. Save project", type = 'ok')
dlg_message("5. Import mzml data files for processing by navigating to File -> import -> results", type = 'ok')
dlg_message("6. Let skyline process.  Export results once complete.  Export reports must have the following headings: Replicate, Protein, Peptide, Area, Retention Time, Start Time and End Time", type = 'ok')
dlg_message("7. Now return to R Studio", type = 'ok')

dlg_message("Next select the export file for importing into skylineR", type = 'ok')
# 
# results_1 <- read_csv(file = file.choose(.)) %>% 
#   clean_names() %>% 
#   rename(`molecule_list_name` = `protein`, `precursor_name` = `peptide`) %>% 
#   filter(grepl("LTR", replicate))
# results_1$area <- sapply(results_1$area, as.numeric) #ensure area column is numeric
# 
# 
# 
# dlg_message("1. Please return to skylineMS software", type = 'ok')
# dlg_message("2. Select all of the transtions in the left window pain and press delete", type = 'ok')
# dlg_message("3. Import the new skylineR_RT_update.csv transition list from csv file by navigating to File -> import -> transition list", type = 'ok')
# dlg_message("4. Save project", type = 'ok')
# dlg_message("5. Export results. Export reports must have the following headings: Replicate, Protein, Peptide, Area, Retention Time, Start Time and End Time", type = 'ok')
# dlg_message("7. Now return to R Studio", type = 'ok')
# 
# dlg_message("Next select the new export file for importing into skylineR in R Studio", type = 'ok')

detach("package:dplyr", unload = TRUE)
library("dplyr")
results_2 <- read_csv(file = file.choose(.))  %>% 
  clean_names() %>% 
  rename(molecule_list_name = protein, precursor_name = peptide) 

filenames <- results_2$replicate %>% unique()

results_2_ltr <- results_2 %>% filter(grepl("LTR", replicate))
results_2_ltr$area <- sapply(results_2_ltr$area, as.numeric) #ensure area column is numeric

rt_boundary_output <- lapply(metabolite_target_list$precursor_name, function(FUNC_LIPID){
  #browser()
  print(FUNC_LIPID)
  rt_boundary <-  filter(results_2_ltr, precursor_name == FUNC_LIPID) %>% filter(!is.na(area)) %>% filter(!grepl("conditioning", replicate)) %>% arrange(retention_time)
  if(nrow(rt_boundary) > 0) {
  start_time <- rt_boundary %>% select(start_time) %>% sapply(as.numeric) %>% min()
  end_time <- rt_boundary %>% select(end_time) %>% sapply(as.numeric) %>% max()
  #end_time <- max(c(rt_boundary$end_time, (rt_boundary$end_time %>% sapply(as.numeric) %>% summary())[5] %>% as.numeric()))
  rt_boundary_filelist <- filenames %>% as_tibble() %>% rename(FileName = value)
  rt_boundary_filelist$FullPeptideName <- FUNC_LIPID
  rt_boundary_filelist$MinStartTime <- round(start_time,2)
  rt_boundary_filelist$MaxEndTime <- round(end_time, 2)
  rt_boundary_filelist
  }
}) %>% bind_rows

rt_boundary_output <- rt_boundary_output %>% filter(!is.na(MinStartTime)) %>% filter(!is.na(MaxEndTime)) #remove any NAs that are in the peak boundary guide

write_csv(rt_boundary_output, (file = paste(project_dir, "/", user_name, "_", project_name, "_skylineR_boundary_update.csv", sep="")))

dlg_message("1. Please return to skylineMS software", type = 'ok')
dlg_message("2. Import the new skylineR_boundary_update.csv transition list from csv file by navigating to File -> import -> peak boundaries", type = 'ok')
dlg_message("4. Save project", type = 'ok')
dlg_message("5. Export results. Export reports must have the following headings: Replicate, Protein, Peptide, Area, Retention Time, Start Time and End Time", type = 'ok')
dlg_message("7. Now return to R Studio and run the lipid_exploreR to QC check data", type = 'ok')