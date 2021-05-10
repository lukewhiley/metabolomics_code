# intensity threshold filter to remove lipids that are not present in the dataset

lipid_intensity_check <- individual_lipid_data_sil_tic_filtered %>% select(!contains("SIL"))

lipid_intensity_list <- lipid_intensity_check %>% select(contains("(")) %>% colnames()

#create a tibble listing the max intensity for each feature
lipid_intensity_max <- lapply(lipid_intensity_list, function(FUNC_INTENSITY){
  #browser()
  lipid_intensity_check %>% select(all_of(FUNC_INTENSITY)) %>% as.matrix() %>% max()
}) %>% c() %>% unlist() %>% as_tibble() %>% add_column(lipid_intensity_list, .before = 1) %>% arrange(value)

write_csv(lipid_intensity_max, "intensity.csv")


