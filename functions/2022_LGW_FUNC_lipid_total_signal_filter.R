# create summed lipid class data
#find all lipids within the same class and sums them together to create a summed class value

# FUNC_data - data to be checked
# FUNC_metabolite_list -> metabolite list in the data
#FUNC_plot_title -> plot title (string)
# FUNC_IS_tag - are internal standards tagged with any string? e.g. SIL (stable isotope labelled) or IS  in their name

lgw_total_summed_signal_filter <- function(FUNC_data,
                                        FUNC_metabolite_list,
                                        FUNC_plot_title,
                                        FUNC_IS_tag 
                                       ){
  
  FUNC_list <- list()

  #step 1 - intensity check for non-IS samples
  
  #create tibble for storing data
  FUNC_summed_data <- FUNC_data %>%
    select(all_of(FUNC_metabolite_list)) %>%
  select(!contains(FUNC_IS_tag))
  
  FUNC_row_sums <- FUNC_data %>% 
    select(!all_of(FUNC_metabolite_list)) %>%
    bind_cols(
  (rowSums(x = FUNC_summed_data,
          na.rm = TRUE) %>% as_tibble() %>% setNames("total_signal"))
    )
  
  FUNC_Q1 <- quantile(FUNC_row_sums$total_signal, .25)
  FUNC_Q3 <- quantile(FUNC_row_sums$total_signal, .75)
  FUNC_IQR <- IQR(FUNC_row_sums$total_signal)
  
  FUNC_list$pass_list <-  subset(FUNC_row_sums, 
         FUNC_row_sums$total_signal > (FUNC_Q1 - 1.5*FUNC_IQR) & 
           FUNC_row_sums$total_signal < (FUNC_Q3 + 1.5*FUNC_IQR))
  
  FUNC_list$fail_list <- subset(FUNC_row_sums, 
                      FUNC_row_sums$total_signal < (FUNC_Q1 - 1.5*FUNC_IQR) | 
                        FUNC_row_sums$total_signal > (FUNC_Q3 + 1.5*FUNC_IQR))

  par(mar=c(1,1,1,1))
  FUNC_list$pre_bp <- boxplot(FUNC_row_sums$total_signal, main = paste0(FUNC_plot_title, " - pre"))
  FUNC_list$post_bp <- boxplot(FUNC_list$pass_list$total_signal, main = paste0(FUNC_plot_title, " - post"))
  FUNC_list$summary <- (FUNC_row_sums$total_signal %>% length()) - (FUNC_list$pass_list$total_signal %>% length())

  #step 2 - repeat for SIL Internal Standards
  
  #create tibble for storing data
  FUNC_summed_data_SIL <- FUNC_data %>%
    select(all_of(FUNC_metabolite_list)) %>%
    select(contains(FUNC_IS_tag))
  
  FUNC_row_sums_SIL <- FUNC_data %>% 
    select(!all_of(FUNC_metabolite_list)) %>%
    bind_cols(
      (rowSums(x = FUNC_summed_data_SIL,
               na.rm = TRUE) %>% as_tibble() %>% setNames("total_signal"))
    )
  
  FUNC_Q1_SIL <- quantile(FUNC_row_sums_SIL$total_signal, .25)
  FUNC_Q3_SIL <- quantile(FUNC_row_sums_SIL$total_signal, .75)
  FUNC_IQR_SIL <- IQR(FUNC_row_sums_SIL$total_signal)
  
  FUNC_list$pass_list_SIL <-  subset(FUNC_row_sums_SIL, 
                                 FUNC_row_sums_SIL$total_signal > (FUNC_Q1_SIL - 1.5*FUNC_IQR_SIL) & 
                                   FUNC_row_sums_SIL$total_signal < (FUNC_Q3_SIL + 1.5*FUNC_IQR_SIL))
  
  FUNC_list$fail_list_SIL <- subset(FUNC_row_sums_SIL, 
                                FUNC_row_sums_SIL$total_signal < (FUNC_Q1_SIL - 1.5*FUNC_IQR_SIL) | 
                                  FUNC_row_sums_SIL$total_signal > (FUNC_Q3_SIL + 1.5*FUNC_IQR_SIL))
  
  par(mar=c(1,1,1,1))
  FUNC_list$pre_bp_SIL <- boxplot(FUNC_row_sums_SIL$total_signal, main = paste0(FUNC_plot_title, " - pre (SIL)"))
  FUNC_list$post_bp_SIL <- boxplot(FUNC_list$pass_list_SIL$total_signal, main = paste0(FUNC_plot_title, " - post (SIL)"))
  FUNC_list$summary_SIL <- (FUNC_row_sums_SIL$total_signal %>% length()) - (FUNC_list$pass_list_SIL$total_signal %>% length())
  
  FUNC_list$unique_fail_list <- unique(c(FUNC_list$fail_list$sample_name, FUNC_list$fail_list_SIL$sample_name))
  FUNC_list
  }
  
