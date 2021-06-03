package_list <- c("plyr", "tidyverse", "janitor", "gridExtra", "ggpubr", "readxl", "cowplot", "scales", "devtools", "metabom8", "shiny", "plotly", "svDialogs", "DataEditR", "htmlwidgets", "httr")
loaded_packages <- lapply(package_list, require, character.only = TRUE)
rm(loaded_packages)

dlg_message("Welcome to skylineR! :-)", type = 'ok')

dlg_message("Please select your project folder", type = 'ok')

project_dir <- rstudioapi::selectDirectory() # save project directory root location
setwd(project_dir) #temporarily set the working directory

#pop-up user input here for project name and user initials
if(exists("project_name") != TRUE){project_name <- dlgInput("what is the name of the project?", "example_project")$res}
if(exists("user_name") != TRUE){user_name <- dlgInput("Insert your initials", "example_initials")$res}

# read in lipid MS target transition information file
dlg_message("please select your lipid MS target transition information file", type = 'ok')
transition_metadata <- read_csv(file.choose(.))
transition_metadata_headers <- colnames(transition_metadata)
transition_metadata <- clean_names(transition_metadata)
metabolite_target_list <- transition_metadata %>% select(precursor_name)

#create list of mzML files
#dlg_message("Please select the folder that contains your mzML files", type = 'ok')
#mzml_file_location <- rstudioapi::selectDirectory() # save project directory root location
#filenames <- list.files(pattern = ".mzML", paste(mzml_file_location)) %>% as_tibble() %>% rename(FileName = value)

dlg_message("1. Please open skylineMS software", type = 'ok')
dlg_message("2. Create new small molecule file", type = 'ok')
dlg_message("3. Import transition list from csv file by navigating to File -> import -> transition list", type = 'ok')
dlg_message("4. Save project", type = 'ok')
dlg_message("5. Import mzml data files for processing by navigating to File -> import -> results", type = 'ok')
dlg_message("6. Let skyline process.  Export results once complete.  Export reports must have the following headings: Replicate, Protein, Peptide, Area, Retention Time, Start Time and End Time", type = 'ok')
dlg_message("7. Now return to R Studio", type = 'ok')

dlg_message("Next select the export file for importing into skylineR", type = 'ok')

results_1 <- read_csv(file = file.choose(.)) %>% 
  clean_names() %>% 
  rename(molecule_list_name = protein, precursor_name = peptide) %>% 
  filter(grepl("LTR", replicate))
results_1$area <- sapply(results_1$area, as.numeric) #ensure area column is numeric

new_rt <- lapply(metabolite_target_list$precursor_name, function(idx_m){
  (filter(results_1, precursor_name == idx_m) %>% filter(!is.na(area)) %>% arrange(desc(area)))[1,] %>% select(precursor_name, retention_time)
}) %>% bind_rows

transition_metadata$explicit_retention_time <- new_rt$retention_time #update transition list with new RT
transition_metadata_2 <- transition_metadata %>% filter(!is.na(explicit_retention_time)) # drop features without a new RT 
metabolite_target_list <- transition_metadata_2 %>% select(precursor_name)
colnames(transition_metadata_2) <- c(transition_metadata_headers)

write_csv(transition_metadata_2, (file = paste(project_dir, "/", user_name, "_", project_name, "_skylineR_RT_update.csv", sep="")))

dlg_message("1. Please return to skylineMS software", type = 'ok')
dlg_message("2. Select all of the transtions in the left window pain and press delete", type = 'ok')
dlg_message("3. Import the new skylineR_RT_update.csv transition list from csv file by navigating to File -> import -> transition list", type = 'ok')
dlg_message("4. Save project", type = 'ok')
dlg_message("5. Export results. Export reports must have the following headings: Replicate, Protein, Peptide, Area, Retention Time, Start Time and End Time", type = 'ok')
dlg_message("7. Now return to R Studio", type = 'ok')

dlg_message("Next select the new export file for importing into skylineR in R Studio", type = 'ok')

results_2 <- read_csv(file = file.choose(.))  %>% 
  clean_names() %>% 
  rename(molecule_list_name = protein, precursor_name = peptide) 

filenames <- results_2$replicate %>% unique()

results_2 <- results_2 %>% filter(grepl("LTR", replicate))
results_2$area <- sapply(results_2$area, as.numeric) #ensure area column is numeric

rt_boundary_output <- lapply(metabolite_target_list$precursor_name, function(idx_m){
  #browser()
  rt_boundary <-  filter(results_1, precursor_name == idx_m) %>% filter(!is.na(area))
  start_time <- rt_boundary %>% select(start_time) %>% sapply(as.numeric) %>% min()
  start_time <- rt_boundary %>% select(end_time) %>% sapply(as.numeric) %>% max()
  #end_time <- max(c(rt_boundary$end_time, (rt_boundary$end_time %>% sapply(as.numeric) %>% summary())[5] %>% as.numeric()))
  rt_boundary_filelist <- filenames
  rt_boundary_filelist$FullPeptideName <- idx_m
  rt_boundary_filelist$MinStartTime <- round(start_time,2)
  rt_boundary_filelist$MaxEndTime <- round(end_time, 2)
  rt_boundary_filelist
}) %>% bind_rows

rt_boundary_output <- rt_boundary_output %>% filter(!is.na(MinStartTime)) %>% filter(!is.na(MaxEndTime)) #remove any NAs that are in the peak boundary guide

write_csv(transition_metadata_2, (file = paste(project_dir, "/", user_name, "_", project_name, "_skylineR_boundary_update.csv", sep="")))

dlg_message("1. Please return to skylineMS software", type = 'ok')
dlg_message("2. Import the new skylineR_boundary_update.csv transition list from csv file by navigating to File -> import -> peak boundaries", type = 'ok')
dlg_message("4. Save project", type = 'ok')
dlg_message("5. Export results. Export reports must have the following headings: Replicate, Protein, Peptide, Area, Retention Time, Start Time and End Time", type = 'ok')
dlg_message("7. Now return to R Studio and run the lipid_exploreR to QC check data", type = 'ok')