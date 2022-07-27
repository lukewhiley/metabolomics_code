# Load packages and libraries

library(tidyverse)
library(janitor)


##################################	
###### USER DATA INPUT HERE#######
##################################

ms_instrument <- paste("pat02") # add instrument here
project <-  paste("paediatric_burns") # add project code here

csv_directory <- paste("/Users/lukewhiley/OneDrive - Murdoch University/projects/paediatric_burns/export_csv_files")
output_directory <- paste("/Users/lukewhiley/OneDrive - Murdoch University/projects/paediatric_burns/export_csv_files/r_output")
final_result_directory <- paste("/Users/lukewhiley/OneDrive - Murdoch University/projects/paediatric_burns/export_csv_files/final_results")



########################################	
## Import Waters TargetLynx MS data ###
#######################################

setwd(paste(csv_directory))


file_list <- as_tibble(list.files(pattern = ".csv")) #list targetlynx export files
file_list$matrix <- grepl("plasma", file_list$value)
file_list$matrix[grep(TRUE, file_list$matrix)] <- paste("plasma")
file_list$matrix[grep(FALSE, file_list$matrix)] <- paste("urine")
file_list <- filter(file_list, !grepl("target_compound",file_list$value))


data_list <- filter(file_list, grepl("data", file_list$value)) #list sample_lists
colnames(data_list) <- c("data_files", "matrix") # set column name

sample_list <- filter(file_list, grepl("sample", file_list$value)) #list sample_lists
colnames(sample_list) <- c("sample_list", "matrix") # set column name

quant_list <- read_csv("target_compounds.csv")
quant_list <- filter(quant_list, !is.na(quant_list$quant))

# combine export files from all plates
targetlynx_output <- NULL

for (x1 in 1:length(unique(file_list$matrix))){
  setwd(paste(csv_directory))

  targetlynx_output <- read_csv(paste(data_list$data_files[x1]))

  matrix <- paste(data_list$matrix[x1])

  targetlynx_output[6,1] <- paste("target") # create colname ---> target
  colnames(targetlynx_output) <- paste(targetlynx_output[6,]) #now use row 4 as header for overall table 
  targetlynx_output <- clean_names(targetlynx_output)

  #remove comma from target names to tidy data
  targetlynx_output$target <- gsub("[,]","",  targetlynx_output$target)


  #create list of all metaolite target transitions (e.g. metabolites and internal standards)
  transition_reference <- as_tibble(unique(targetlynx_output$target[grep("Compound", targetlynx_output$target)])) # unique list of transitions
  transition_reference <- as_tibble(transition_reference$value[!grepl("Summary", transition_reference$value)]) # removes unwanted line

  #create list of all target metabolites (remove internal standard data)
  target_metabolites <- transition_reference
  target_metabolites <- as_tibble(target_metabolites$value[!grepl("13C", target_metabolites$value)])

  #create list of targets names e.g. Compound 1: Picolinic acid converted to picolinic acid, this helps later when writing csv and png
  target_compounds <- quant_list$metabolite

  # loop to insert a metabolite name column for each meatbolite segment
  for (x2 in 1:(length(transition_reference$value))) {
      start_reference <- grep(paste(transition_reference$value[x2]), targetlynx_output$target)+1 # find start row of each target metabolite
  
      if (x2 < length(transition_reference$value)) {
      end_reference <- grep(paste(transition_reference$value[x2+1]), targetlynx_output$target)-1 # find end row of each target metabolite (finds next sergment and subtracts 1)
    
  } else {
      end_reference <- length(targetlynx_output$target)
    
  }
  
      for (x3 in 1:length(start_reference)){

      targetlynx_output$target[(start_reference[x3]+2):(end_reference[x3])] <- paste(transition_reference$value[x2])

  }
}

  #format the data

  setwd(paste(output_directory))
  dir.create(paste(matrix))

  setwd(paste(output_directory, "/", matrix, sep=""))
  
  # step_1 - extract data from samples
  for (x4 in 1:length(target_metabolites$value)){
    metabolite_data <- targetlynx_output[grep(target_metabolites$value[x4], targetlynx_output$target),]
    metabolite_data <- metabolite_data[-1,]
    metabolite_data$target <- paste(target_compounds[x4] )
    metabolite_data$ms_instrument <- paste(ms_instrument)
    metabolite_data$project <- paste(project)
    metabolite_data$quant <- quant_list$quant[which(target_compounds[x4]==quant_list$metabolite)]
    
    #metabolite_data$analysis_day <- paste(analysis_day)
  
  
    write_csv(metabolite_data, paste(matrix, "_", target_compounds[x4],  "_", project, ".csv", sep = ""))
  }

    meta_data <- select(metabolite_data, name, sample_text, type, std_conc,  ms_instrument)


    ##############################################################
    ######### re-import and create sub-sheets####################
    #############################################################

    
    setwd(paste(output_directory, "/", matrix, sep=""))
    
    csv_file_list <- as_tibble(list.files(pattern = ".csv")) #list targetlynx export files
    colnames(csv_file_list) <- paste("csv_files") # set column name


    combined_output <- meta_data # create meta data framework for final table

    for (x5 in 1:length(csv_file_list$csv_files)){
      read_csv(paste(csv_file_list$csv_files[x5])) %>%
        select(name, conc, area, type, quant) -> metabolite_conc
  
    metabolite <- sub(".csv", "", paste(csv_file_list$csv_files[x5])) # get metabolite name
  
    if (metabolite_conc$quant[1] == TRUE) {
    
      metabolite <- sub("_covid", "_conc_uM", metabolite)
      metabolite_conc$conc[grep("Analyte", metabolite_conc$type)] <-  metabolite_conc$conc[grep("Analyte", metabolite_conc$type)]*2 # adjust for dilution factor
      metabolite_conc <- select(metabolite_conc, name, conc) 
      metabolite_conc <- rename(metabolite_conc, !!metabolite := "conc")
    } else {
      metabolite <- sub("_covid", "_peak_area", metabolite)
      metabolite_conc <- select(metabolite_conc, name, area)
      metabolite_conc <- rename(metabolite_conc, !!metabolite := "area")
  } 
  
    
    
      combined_output <- left_join(combined_output, metabolite_conc, by = "name")
  
  }

setwd(paste(final_result_directory))

   #output analyte (study sample) data
   samples <- filter(combined_output, type == "Analyte")
    write_csv(samples, paste(Sys.Date(), project, matrix, "samples_data.csv", sep = "_"))

  #output standard data
  standards <- filter(combined_output, type == "Standards")
  write_csv(samples, paste(Sys.Date(), project, matrix, "qc_data.csv", sep = "_"))
  #output QC data
  qcs <- filter(combined_output, type == "QC")
  write_csv(combined_output, paste(Sys.Date(), project, matrix, "standard_data.csv", sep = "_"))

  }

