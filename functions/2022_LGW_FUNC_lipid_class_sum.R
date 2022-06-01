# create summed lipid class data
#find all lipids within the same class and sums them together to create a summed class value

lgw_lipid_class_data_summed <- function(FUNC_data,
                                        FUNC_metabolite_list
                                       ){
#find unique lipid class
  lipid_class <- FUNC_metabolite_list
  lipid_class <- sub("\\(.*", "", lipid_class) %>% unique()

  #create tibble for storing data
  FUNC_class_data <- FUNC_data %>%
    select(!all_of(FUNC_metabolite_list))
  
  
  #loop to complete rowSums
    for(idx_lipid_class in lipid_class){
      #print(idx_lipid_class)
      #browser()
      FUNC_class_data <- bind_cols(
        FUNC_class_data,
       (rowSums(x = FUNC_data %>%
              select(contains(paste0(idx_lipid_class, "("))),
              na.rm = TRUE) %>% as_tibble() %>%
          setNames(idx_lipid_class))
     )
    }
  
#add column with total summed data (not class specific)

FUNC_class_data <-  bind_cols(
  FUNC_class_data, 
  (rowSums(x = FUNC_data %>%
             select(!contains("sample")),
           na.rm = TRUE) %>% as_tibble() %>%
     setNames("total_summed_signal"))
)

FUNC_class_data
    
}
