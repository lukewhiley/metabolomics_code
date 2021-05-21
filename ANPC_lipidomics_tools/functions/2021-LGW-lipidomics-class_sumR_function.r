# create summed lipid class data
#find all lipids within the same class and sums them together to create a summed class value

create_lipid_class_data_summed <- function(individual_lipid_data){
  #browser()
  lipid_class <- individual_lipid_data %>% select(contains("(")) %>% colnames() 
  lipid_class <- sub("\\(.*", "", lipid_class) %>% unique()
  lipid_class <- lipid_class[!grepl("sampleID", lipid_class)] %>% as_tibble()
  
  temp_class_data <- apply(lipid_class, 1, function(func_lipid_class){
    #browser()
    class_targets <- which(sub("\\(.*", "", colnames(individual_lipid_data)) == func_lipid_class) # find the columns in each lipid class
    temp_class_data <- individual_lipid_data %>% select(all_of(class_targets))  %>% mutate(rowsum = rowSums(., na.rm = TRUE)) %>% select(rowsum)
    colnames(temp_class_data) <- func_lipid_class
    temp_class_data
  }) %>% bind_cols()
  
  temp_class_data <- cbind(individual_lipid_data$sampleID, temp_class_data) %>% rename("sampleID" = "individual_lipid_data$sampleID")
  temp_class_data
}
