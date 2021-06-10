###################################
###### summed SIL checkR #########
##################################

#this script fits into the lipidomics lipid exploreR. 
# the script will look for samples that have been got preparation errors
# the script sums the total intensity from SIL internal standard
# Samples are removed if the summed intensity is either > or < x standard deviations (SD) above or below the mean summed SIL internal standard signal of the dataset. 
# users define the SD cuttoff (x)
# Sample SIL intensity > x SD indicates high concentration of SIL, so likely as a result of excess SIL IS added
# Sample SIL intensity < x SD below mean indicates too little volume added of SIL internal standard added
# users can select if failed samples are removed

dlg_message("Internal standard check. This next step will assess the internal standards accross all of the samples. If internal standards have been incorrectly added the summed signal intensity will be too low/high.", type = 'ok')

total_summed_sil <- apply(individual_lipid_data %>% select(sampleID), 1, function(summedSIL){
  temp_data <- individual_lipid_data %>% filter(sampleID == summedSIL) %>% select(-sampleID) %>% select(contains("SIL")) %>% rowSums(na.rm = TRUE)
}) %>% c() %>% as_tibble() %>%  add_column(individual_lipid_data$sampleID, .before = 1) %>% 
  rename(SIL_TIC = value, sampleID = "individual_lipid_data$sampleID")


total_summed_sil <- new_project_run_order %>% left_join(total_summed_sil, by = "sampleID") %>% arrange(injection_order)
total_summed_sil$sample_idx <- c(1:nrow(total_summed_sil))
total_summed_sil$LOG_SIL_TIC <- log(total_summed_sil$SIL_TIC)

# while loop here
sil_check_status <- "change"
while(sil_check_status == "change"){

# #flag samples with SIL x standard deviations below mean
# temp_answer <- dlgInput("What do you wish to set for the fail cut off filter.  x number of median absolute deviations from the median", "e.g. recommended default x = 4")$res
# while(is.na(as.numeric(temp_answer))){
#   temp_answer <- dlgInput("You did not enter a numeric value.  What do you wish to set for the fail cut off filter.  x number of median absolute deviations from the median", "e.g. recommended default x = 4")$res
# }
# 
#  median_sil_tic <- median(total_summed_sil$SIL_TIC)
# mad_sil_tic <- mad(total_summed_sil$SIL_TIC)
# sil_cut_off_lower <- median_sil_tic - (as.numeric(temp_answer)*mad_sil_tic)
# sil_cut_off_upper <- median_sil_tic + (as.numeric(temp_answer)*mad_sil_tic)

 
 temp_answer <- dlgInput("What do you wish to set for the fail cut off filter.  x % from the median", "e.g. recommended default x = 50")$res
 while(is.na(as.numeric(temp_answer))){
   temp_answer <- dlgInput("You did not enter a numeric value.  What do you wish to set for the fail cut off filter.  x % from the median", "e.g. recommended default x = 50")$res
 }
 
  median_sil_tic <- median(total_summed_sil$SIL_TIC)
  #inter_quantile_range <- as.numeric(quantile(total_summed_sil$SIL_TIC, 0.75)) - as.numeric(quantile(total_summed_sil$SIL_TIC, 0.25))
  # sil_cut_off_lower <- median_sil_tic - (as.numeric(temp_answer)*inter_quantile_range)
  # sil_cut_off_upper <- median_sil_tic + (as.numeric(temp_answer)*inter_quantile_range)
  
  sil_cut_off_lower <- median_sil_tic - (median_sil_tic*as.numeric(temp_answer)/100)
  sil_cut_off_upper <- median_sil_tic + (median_sil_tic*as.numeric(temp_answer)/100)

#create lists of which samples have failed the SIL internal standard check
sil_qc_fail <- total_summed_sil$sampleID[which(total_summed_sil$SIL_TIC < sil_cut_off_lower | total_summed_sil$SIL_TIC > sil_cut_off_upper)] %>% as_tibble %>% rename(sampleID = value)
sil_qc_fail$fail_point <- "sil"
sil_qc_fail_ltr <- sil_qc_fail %>% filter(grepl("LTR", sampleID))
sil_qc_fail_samples <- sil_qc_fail %>% filter(!grepl("LTR", sampleID))

# visualise for reports
total_summed_sil$removed <- "pass_qc"
total_summed_sil$removed[total_summed_sil$sampleID %in% sil_qc_fail$sampleID] <- "removed"

total_summed_sil_removed <- total_summed_sil %>% filter(grepl("removed", removed))
total_summed_sil_pass <- total_summed_sil %>% filter(grepl("pass_qc", removed))

# create a plate list ID
plate_number <- unique(plateID) %>% substr(14,14) %>% unique()
plate_number <- seq(1, length(plate_number))

plateIDx <- lapply(unique(plateID), function(FUNC_plateID){
  #browser()
  grep(FUNC_plateID, total_summed_sil$plateID)[1]}) %>% unlist()


#set y axis limits
if(sil_cut_off_lower < min(total_summed_sil$SIL_TIC)){
  y_limit_lower <- log(sil_cut_off_lower-(sil_cut_off_lower/100*25))
}
if(sil_cut_off_lower > min(total_summed_sil$SIL_TIC)){
  y_limit_lower <- log(min(total_summed_sil$SIL_TIC)-(min(total_summed_sil$SIL_TIC)/100*25))
}
if(sil_cut_off_upper > max(total_summed_sil$SIL_TIC)){
  y_limit_upper <- log(max(total_summed_sil$SIL_TIC)+(max(total_summed_sil$SIL_TIC)/100*25))
}
if(sil_cut_off_upper < max(total_summed_sil$SIL_TIC)){
  y_limit_upper <- log(max(total_summed_sil$SIL_TIC)+(max(total_summed_sil$SIL_TIC)/100*25))
}

# create a layout list of extra lines to add
p_threshold_lines <- list(list(type='line', x0= min(total_summed_sil$sample_idx), x1= (max(total_summed_sil$sample_idx)+10), y0=log(sil_cut_off_lower), y1=log(sil_cut_off_lower),
                          line=list(dash='dot', width=3, color = '#FF0000')),
                     list(type='line', x0= min(total_summed_sil$sample_idx), x1= (max(total_summed_sil$sample_idx)+10), y0=log(sil_cut_off_upper), y1=log(sil_cut_off_upper),
                          line=list(dash='dot', width=3, color = '#FF0000')),
                     list(type='line', x0= min(total_summed_sil$sample_idx), x1= (max(total_summed_sil$sample_idx)+10), y0=log(median_sil_tic), y1=log(median_sil_tic),
                          line=list(dash='dot', width=3, color = '#000000'))
)
p_plate_list <- lapply(plateIDx[2:length(plateIDx)], function(FUNC_P_PLATE_LIST){
   list(type='line', x0 = FUNC_P_PLATE_LIST, x1= FUNC_P_PLATE_LIST, y0=y_limit_lower-log(median_sil_tic)*.05, y1=y_limit_upper+log(median_sil_tic)*.05,
            line=list(dash='dot', width=2, color = '#808080'))
})

#only add plate lines if multiple plates exist
# if(is.na(plateIDx)){
#   p_plot_lines <- p_threshold_lines
# }

if(length(plateIDx) == 1){
p_plot_lines <- p_threshold_lines
}

if(length(plateIDx) > 1){
  p_plot_lines <- c(p_threshold_lines, p_plate_list)
}

#create a list of axis settings for plot_ly
x_axis_settings <- list(
  zeroline = FALSE,
  showline = TRUE,
  linecolor = toRGB("black"),
  linewidth = 2,
  showgrid = FALSE,
  range = c(0, max(total_summed_sil$sample_idx)+1),
  title = "Sample index"
)

y_axis_settings <- list(
  zeroline = FALSE,
  showline = TRUE,
  linecolor = toRGB("black"),
  linewidth = 2,
  showgrid = TRUE,
  title = "Lipid total ion count (Log)",
  range = c(y_limit_lower-log(median_sil_tic)*.05, 
            y_limit_upper+log(median_sil_tic)*.05)
)

p <- plot_ly(
  type = "scatter", mode = "markers", data = total_summed_sil_pass, x = ~sample_idx, y = ~LOG_SIL_TIC, text = ~sampleID, color = ~removed, colors = c('#1E90FF', '#FF0000'), 
  marker = list(size = 7, color = '#1E90FF', opacity = 0.5,
                line = list(color = '#000000',width = 1))
 ) %>% 
  add_trace(type = "scatter", data = total_summed_sil_removed, x = ~sample_idx, y = ~LOG_SIL_TIC, text = ~sampleID, color = ~removed, 
            marker = list(size = 8, color = '#FF0000')
            ) %>%
  layout(xaxis = x_axis_settings,
         yaxis = y_axis_settings
         ) %>%
  layout(shapes=p_plot_lines)

#create html widget and display it in the users internet browser
sil_check_p <- p

saveWidget(sil_check_p, file = paste(project_dir_html, "/", project_name, "_", user_name, "_SIL_check_plot.html", sep=""))# save plotly widget
browseURL(paste(project_dir_html, "/", project_name, "_", user_name, "_SIL_check_plot.html", sep="")) #open plotly widget in internet browser

#sil_qc_fail - ask the user if they wish to continue or change the threshold
sil_check_status <- "blank"
while(sil_check_status != "continue" & sil_check_status != "change"){
  sil_check_status <- dlgInput(paste(nrow(sil_qc_fail), "samples FAILED the SIL QC check.  continue or change the exclusion threshold?"), "continue/change")$res
}
}

#sil_qc_fail - ask the user if they wish to remove all/none/samples/LTR which failed the QC check
temp_answer <- "blank"
while(temp_answer != "all" & temp_answer != "none" & temp_answer != "samples" & temp_answer != "LTR"){
  temp_answer <- dlgInput(paste("of the ", nrow(sil_qc_fail), "FAILED samples.  ",  nrow(sil_qc_fail_ltr),"were LTRs.  Do you want to remove failed samples?"), "all/none/samples/LTR")$res
  if(temp_answer == "all"){individual_lipid_data_sil_filtered <- individual_lipid_data %>% filter(!sampleID %in% sil_qc_fail$sampleID)}
  if(temp_answer == "samples"){individual_lipid_data_sil_filtered <- individual_lipid_data %>% filter(!sampleID %in% sil_qc_fail_samples$sampleID)}
  if(temp_answer == "LTR"){individual_lipid_data_sil_filtered <- individual_lipid_data %>% filter(!sampleID %in% sil_qc_fail_ltr$sampleID)}
  if(temp_answer == "none"){individual_lipid_data_sil_filtered <- individual_lipid_data}
}

#tidy up environment

 # remove_list <- c("mad_sil_tic", "median_sil_tic", "plate_number", "sil_check_status", "sil_cut_off_lower", "sil_cut_off_upper", "remove_list",
 #                 "p", "total_summed_sil")
 # rm(list = remove_list)
      

