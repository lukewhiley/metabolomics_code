#impute function

lgw_lipid_conc_calc <- function(FUNC_data,
                              FUNC_metabolite_list,
                             FUNC_SIL_guide_path,
                             FUNC_conc_guide_path)
  {

 FUNC_list <- list()
  
  FUNC_SIL_guide <- read_csv(
   file = FUNC_SIL_guide_path, 
    show_col_types = FALSE) %>% 
    clean_names()
  

  FUNC_conc_guide <- read_csv(file = FUNC_conc_guide_path, show_col_types = FALSE) %>% 
    clean_names()
  
lipid_features <- FUNC_metabolite_list[!grepl("SIL",FUNC_metabolite_list)] 
    
    #set empty list
FUNC_list$out <- FUNC_data %>%
      select(sample_name)
    
    for(idx_lipid in lipid_features){
      #find appropriate SIL IS from SIL_guide
      SIL_used <- FUNC_SIL_guide %>%
        filter(precursor_name == idx_lipid) %>%
        select(note) %>% paste0()
      
      #find SIL concentration factor from template
      SIL_concentration_factor <- FUNC_conc_guide %>%
        filter(sil_name == all_of(SIL_used)) %>%
        select(concentration_factor) %>% as.numeric()
      
      FUNC_list$out <- FUNC_list$out %>% 
        #left join to keep samples and rows aligned
        left_join(FUNC_data %>%
                    #select lipid and its SIL IS columns
                    select(sample_name, all_of(idx_lipid), all_of(SIL_used)) %>% 
                    #create ratio between lipid and SIL IS
                    add_column("ratio" = (FUNC_data[[idx_lipid]]/
                                            FUNC_data[[SIL_used]])*
                                 #multiply by concentration factor (from template)
                                 SIL_concentration_factor) %>%
                    #drop lipid and ratio columns
                    select(-all_of(idx_lipid), -all_of(SIL_used)) %>% 
                    #rename column
                    rename_with(~all_of(idx_lipid), "ratio"),
                  by= "sample_name"
        )
    }
  }
  
FUNC_data
}