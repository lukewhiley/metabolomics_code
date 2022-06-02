# create summed lipid class data
#find all lipids within the same class and sums them together to create a summed class value

# FUNC_data - data to be checked
# FUNC_metabolite_list -> metabolite list in the data
#FUNC_plot_title -> plot title (string)
# FUNC_IS_tag - are internal standards tagged with any string? e.g. SIL (stable isotope labelled) or IS  in their name

lgw_summed_class_signal_filter <- function(FUNC_data,
                                           FUNC_metabolite_list,
                                           FUNC_plot_title,
                                           FUNC_IS_tag,
                                           FUNC_OPTION_NON_SIL_filter,
                                           FUNC_OPTION_SIL_filter,
                                           FUNC_OPTION_show_boxplots
){
  
  FUNC_list <- list()
  FUNC_list$pass_list <- list()
  FUNC_list$fail_list <- list()
  FUNC_list$summary <- list()
  FUNC_list$pass_list_SIL <- list()
  FUNC_list$fail_list_SIL <- list()
  FUNC_list$summary_SIL <- list()
  FUNC_list$full_list <- list()
  
  #step 1 - intensity check for non-IS samples
  
  lipid_class <- FUNC_metabolite_list
  lipid_class <- sub("\\(.*", "", lipid_class) %>% unique()

  for(idx_lipid in lipid_class){
  
    if(!grepl(FUNC_IS_tag, idx_lipid)){
  #create tibble for storing data
  FUNC_summed_data <- FUNC_data %>%
    select(contains(paste0(idx_lipid, "("))) %>%
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
  
  #set pass filter to keep all features that are lower than [quartile 1 - 2* the interquartile range] - looking for samples that have a clear overall low signal (i.e. miss injection)
  FUNC_list$pass_list[[idx_lipid]] <-  subset(FUNC_row_sums, 
         FUNC_row_sums$total_signal > (FUNC_Q1 - FUNC_OPTION_NON_SIL_filter*FUNC_IQR)) #& FUNC_row_sums$total_signal < (FUNC_Q3 + FUNC_OPTION_NON_SIL_filter*FUNC_IQR)) 
  
  FUNC_list$pass_list[[idx_lipid]] <- FUNC_list$pass_list[[idx_lipid]] %>%
    add_column("class" = rep(idx_lipid, nrow(FUNC_list$pass_list[[idx_lipid]]))) %>%
    add_column("filter" = rep("total signal", nrow(FUNC_list$pass_list[[idx_lipid]])))
  
  #create fail list
  FUNC_list$fail_list[[idx_lipid]] <- subset(FUNC_row_sums, 
                      FUNC_row_sums$total_signal < (FUNC_Q1 - FUNC_OPTION_NON_SIL_filter*FUNC_IQR)) #| FUNC_row_sums$total_signal > (FUNC_Q3 + FUNC_OPTION_NON_SIL_filter*FUNC_IQR) 
  
  FUNC_list$fail_list[[idx_lipid]] <- FUNC_list$fail_list[[idx_lipid]] %>%
    add_column("class" = rep(idx_lipid, nrow(FUNC_list$fail_list[[idx_lipid]]))) %>%
    add_column("filter" = rep("total signal", nrow(FUNC_list$fail_list[[idx_lipid]])))

 
  if(FUNC_OPTION_show_boxplots == TRUE){
    par(mar=c(1,1,1,1))
    boxplot(FUNC_row_sums$total_signal, main = paste0(FUNC_plot_title, " - pre"))
    boxplot(FUNC_list$pass_list$total_signal, main = paste0(FUNC_plot_title, " - post"))
  }
  
  FUNC_list$summary[[idx_lipid]] <- (FUNC_row_sums$total_signal %>% length()) - (FUNC_list$pass_list[[idx_lipid]]$total_signal %>% length())
    }
 
     #step 2 - repeat for SIL Internal Standards
 
    if(grepl(FUNC_IS_tag, idx_lipid)){ 
  #create tibble for storing data
  FUNC_summed_data_SIL <- FUNC_data %>%
    select(contains(paste0(idx_lipid, "("))) %>%
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
  
  FUNC_list$pass_list_SIL[[idx_lipid]] <-  subset(FUNC_row_sums_SIL, 
                                 FUNC_row_sums_SIL$total_signal > (FUNC_Q1_SIL - FUNC_OPTION_SIL_filter*FUNC_IQR_SIL)) #& FUNC_row_sums_SIL$total_signal < (FUNC_Q3_SIL + FUNC_OPTION_SIL_filter*FUNC_IQR_SIL)) #wider threshold for IS - should be more consistent so smaller IQR
  
  FUNC_list$pass_list_SIL[[idx_lipid]] <- FUNC_list$pass_list_SIL[[idx_lipid]] %>%
    add_column("class" = rep(idx_lipid, nrow(FUNC_list$pass_list_SIL[[idx_lipid]]))) %>%
    add_column("filter" = rep("SIL signal", nrow(FUNC_list$pass_list_SIL[[idx_lipid]])))
  
  FUNC_list$fail_list_SIL[[idx_lipid]] <- subset(FUNC_row_sums_SIL, 
                                FUNC_row_sums_SIL$total_signal < (FUNC_Q1_SIL - FUNC_OPTION_SIL_filter*FUNC_IQR_SIL)) #| FUNC_row_sums_SIL$total_signal > (FUNC_Q3_SIL + FUNC_OPTION_SIL_filter*FUNC_IQR_SIL)) #wider threshold for IS - should be more consistent so smaller IQR
  
  FUNC_list$fail_list_SIL[[idx_lipid]] <- FUNC_list$fail_list_SIL[[idx_lipid]] %>%
    add_column("class" = rep(idx_lipid, nrow(FUNC_list$fail_list_SIL[[idx_lipid]]))) %>%
    add_column("filter" = rep("SIL signal", nrow(FUNC_list$fail_list_SIL[[idx_lipid]])))
  
  if(FUNC_OPTION_show_boxplots == TRUE){
    par(mar=c(1,1,1,1))
  FUNC_list$pre_bp_SIL <- boxplot(FUNC_row_sums_SIL$total_signal, main = paste0(FUNC_plot_title, " - pre (SIL)"))
  FUNC_list$post_bp_SIL <- boxplot(FUNC_list$pass_list_SIL$total_signal, main = paste0(FUNC_plot_title, " - post (SIL)"))
  }
  
  FUNC_list$summary_SIL[[idx_lipid]] <- (FUNC_row_sums_SIL$total_signal %>% length()) - (FUNC_list$pass_list_SIL[[idx_lipid]]$total_signal %>% length())
    }
  }
  
  FUNC_list$fail_list$all <- bind_rows(FUNC_list$fail_list)
  FUNC_list$fail_list_SIL$all <- bind_rows(FUNC_list$fail_list_SIL)
  
  FUNC_list$full_list <- bind_rows(FUNC_list$fail_list$all, 
                                   FUNC_list$fail_list_SIL$all)
  
  
  FUNC_list$unique_fail_list <- unique(FUNC_list$full_list$sample_name)
  FUNC_list
  }
  
