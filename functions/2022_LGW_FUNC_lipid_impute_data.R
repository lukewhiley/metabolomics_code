#impute function

lgw_lipid_impute <- function(FUNC_data,
                              FUNC_metabolite_list,
                             FUNC_option_impute_missing_data)
  {

#step 1 = replace all NA with 0 (missing value) prior to impute
  FUNC_data[is.na(FUNC_data)] <- 0

# #find missing value % in current data frame
# temp_frame <- master_list$impute_data[[idx_data]] %>% select(!contains("sample")) 

# total_values <- dim(temp_frame)[1] * dim(temp_frame)[2]

missing_values <- length(which(FUNC_data==0))
  
if(missing_values != 0 & FUNC_option_impute_missing_data == TRUE){
  #find idx of columns containing a 0 value:
  zero_value_column_idx <- which(
    colSums(
      FUNC_data == 0) > 0) %>% 
    as.matrix() %>% 
    c()
  metadata_idx <- grep("sample",colnames(FUNC_data))
  
  `%notin%` <- Negate(`%in%`) #create %notin% function
  
  zero_value_column_idx <- zero_value_column_idx[zero_value_column_idx %notin% metadata_idx] #remove any column index that are metadata
  
  #min value impute
  for(idx_num_impute in zero_value_column_idx){
    temp_impute_data <- FUNC_data[,idx_num_impute]
    non_zero_values_idx <- which(temp_impute_data > 0) 
    zero_values_idx <- which(temp_impute_data == 0)
    impute_value <- min(temp_impute_data[non_zero_values_idx,])/2
    FUNC_data[zero_values_idx,idx_num_impute] <- impute_value
  }
  
  }

FUNC_data
}