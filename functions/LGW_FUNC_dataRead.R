# function to read in data from a single csv
#csv must contain sample annotation data and metabolite data
#annotation data must be tagged for filtering (e.g. with test string "sample")

#FUNC_datapath = STRING = datapath of csv file
#FUNC_annotation_tag = STRING = tag used to identify annotation data

LGW_FUNC_dataRead <- function(FUNC_datapath, FUNC_annotation_tag){
  
  function_temp_data <- read_csv(FUNC_datapath, show_col_types = FALSE) %>%
    clean_names()
  # function_temp_features <- function_temp_data %>% select(-contains("sample")) %>%
  #   colnames()
  
  list(annotation = function_temp_data %>% select(contains(FUNC_annotation_tag)),
       data = function_temp_data %>% select(!contains(FUNC_annotation_tag)),
       features = function_temp_data %>% select(!contains(FUNC_annotation_tag)) %>% colnames()
       )
}