# create summed lipid class data
#find all lipids within the same class and sums them together to create a summed class value

# FUNC_data - data to be checked
# FUNC_metabolite_list -> metabolite list in the data
#FUNC_plot_title -> plot title (string)
# FUNC_IS_tag - are internal standards tagged with any string? e.g. SIL (stable isotope labelled) or IS  in their name

lgw_lipid_missing_value_filter <- function(FUNC_data,
                                        FUNC_metabolite_list,
                                        FUNC_IS_tag,
                                        FUNC_OPTION_missing_value_threshold,
                                        FUNC_OPTION_intensity_threshold
                                       ){
  
  FUNC_list <- list()
  
  #create tibble for storing data
  FUNC_summed_data <- FUNC_data %>%
    select(all_of(FUNC_metabolite_list)) %>%
    select(!contains(FUNC_IS_tag))
  
  zero_value_threshold <- (nrow(FUNC_summed_data)/100)*FUNC_OPTION_missing_value_threshold
  
  
  #step 1 - replace all NAs with 0
  
  FUNC_summed_data[is.na(FUNC_summed_data)] <- 0

  #step 2 - remove features with x amount of missing values
  
  FUNC_list$missing_value_fail <- names(FUNC_summed_data)[which(colSums(na.rm = TRUE,
                                                                        FUNC_summed_data < FUNC_OPTION_intensity_threshold) > zero_value_threshold) %>% 
                                                            as.numeric()]

  FUNC_list
  }

