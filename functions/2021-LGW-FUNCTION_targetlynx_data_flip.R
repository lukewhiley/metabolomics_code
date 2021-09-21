#function to read in targetlynx txt file and then spit out a flipped tibble
# package requirements 
  # -> janitor
  # -> tidyverse

require(tidyverse)
require(janitor)

targetlynx_flippR <- function(FUNC_txt_file_path){

flippr_func_data <- read.delim(paste(FUNC_txt_file_path),
                               sep = "\t",
                               header=FALSE)

#find the row which contains the column numbers
for(idx_row in 1:10){
  temp_idx <- grep("Name", flippr_func_data[idx_row,])
  if(length(temp_idx==1)){
    colnames_row_idx <- c(idx_row)
  }
}
#set names from the found row idx
names(flippr_func_data) <- flippr_func_data[colnames_row_idx[1],]

#clean names - requires janitor package
flippr_func_data <- flippr_func_data %>% clean_names()

# create lists of targets -------------------------------------------------

flippr_func_compound_list <- unique(flippr_func_data$x[grep("Compound", flippr_func_data$x)])
flippr_func_compound_list <- flippr_func_compound_list[-grep("Report", flippr_func_compound_list)]

flippr_func_metabolite_list <- sub(".*:  ", "", flippr_func_compound_list)

# flip targetlynx_data --------------------------------------------------------

compound_start_idx <- grep(paste(flippr_func_compound_list,collapse="|"), 
                     flippr_func_data$x)

compound_end_idx <- c(compound_start_idx[-1]-1, nrow(flippr_func_data))


for(idx_analyte in 1:length(compound_start_idx)){
  range_idx <- seq(compound_start_idx[idx_analyte], 
                   compound_end_idx[idx_analyte], 
                   by = 1)
  if(idx_analyte == 1){
     temp_data <- flippr_func_data[range_idx,] %>%
       select(name, sample_text, type, conc)
  }
  
  if(idx_analyte > 1){
    temp_data <- flippr_func_data[range_idx,] %>%
      select(name, conc)
  }
  
  start_idx <- which(temp_data$name == "Name")
  temp_data <- temp_data[-(1:start_idx),]
  
  if(idx_analyte == 1){
    flippr_func_data_output <- temp_data
  }
  
  if(idx_analyte > 1){
  flippr_func_data_output <- left_join(flippr_func_data_output, temp_data, by = "name")
  }
}

names(flippr_func_data_output) <- c("sampleID", "sampleText", "sampleType", flippr_func_metabolite_list)

flippr_func_data_output

}



