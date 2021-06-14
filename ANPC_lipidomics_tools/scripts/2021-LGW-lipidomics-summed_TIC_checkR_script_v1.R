###################################
###### summed TIC checkR #########
##################################

#this script fits into the lipidomics lipid exploreR. 
# the script will look for samples that have been got preparation errors
# the script sums the total intensity from all endogenous targets (e.g. not internal standards)
dlg_message("Summed TIC check. This next step will assess the summed TIC accross all of the samples. If samples have been incorrectly prepared the summed TIC intensity will be too low/high. Used to remove obvious outliers only.", type = 'ok')

total_summed_tic <- apply(individual_lipid_data_sil_filtered %>% select(sampleID), 1, function(summedTIC){
  #browser()
  temp_data <- individual_lipid_data_sil_filtered %>% filter(sampleID == summedTIC) %>% select(-sampleID, -plateID) %>% select(!contains("SIL")) %>% rowSums(na.rm = TRUE)
}) %>% c() %>% as_tibble() %>%  add_column(individual_lipid_data_sil_filtered$sampleID, .before = 1) %>% 
  rename(summed_TIC = value, sampleID = "individual_lipid_data_sil_filtered$sampleID")

total_summed_tic <- new_project_run_order %>% left_join(total_summed_tic, by = "sampleID") %>% arrange(injection_order) %>% filter(!is.na(summed_TIC))
total_summed_tic$sample_idx <- c(1:nrow(total_summed_tic))
total_summed_tic$LOG_summed_TIC <- log(total_summed_tic$summed_TIC)

#while loop here
tic_check_status <- "change"
while(tic_check_status == "change"){

#flag samples with SIL x number of standard deviations below mean
# temp_answer <- dlgInput("What do you wish to set for the fail cut off filter.  x number of median absolute deviations from the median", "e.g. recommended default x = 4")$res
# while(is.na(as.numeric(temp_answer))){
#   temp_answer <- dlgInput("You did not enter a numeric value.  What do you wish to set for the fail cut off filter.  x number of median absolute deviations from the median", "e.g. recommended default x = 4")$res
# }
# 
# temp_answer <- as.numeric(temp_answer)
# median_summed_tic <- median(total_summed_tic$summed_TIC)
# mad_summed_tic <- mad(total_summed_tic$summed_TIC)
# tic_cut_off_lower <- median_summed_tic - (as.numeric(temp_answer)*mad_summed_tic)
  
  
  temp_answer <- dlgInput("What do you wish to set for the fail cut off filter.  x % from the median", "e.g. recommended default x = 50")$res
  while(is.na(as.numeric(temp_answer))){
    temp_answer <- dlgInput("You did not enter a numeric value.  What do you wish to set for the fail cut off filter.  x % from the median", "e.g. recommended default x = 50")$res
  }
  
  median_summed_tic <- median(total_summed_tic$summed_TIC)
  # inter_quantile_range <- as.numeric(quantile(total_summed_tic$summed_TIC, 0.50)) - as.numeric(quantile(total_summed_tic$summed_TIC, 0.25))
  # tic_cut_off_lower <- median_summed_tic - (as.numeric(temp_answer)*inter_quantile_range)
  
  median_summed_tic <- median(total_summed_tic$summed_TIC)
  tic_cut_off_lower <- median_summed_tic - (median_summed_tic*as.numeric(temp_answer)/100)


tic_qc_fail <- total_summed_tic$sampleID[which(total_summed_tic$summed_TIC < tic_cut_off_lower)] %>% as_tibble %>% rename(sampleID = value)
tic_qc_fail$fail_point <- "tic"
tic_qc_fail_ltr <- tic_qc_fail %>% filter(grepl("LTR", sampleID))# create tibble of failed LTRs
tic_qc_fail_samples <- tic_qc_fail %>% filter(!grepl("LTR", sampleID)) # create tibble of failed samples (not LTRS)


#visualise for reports
total_summed_tic$removed <- "pass_qc"
total_summed_tic$removed[total_summed_tic$sampleID %in% tic_qc_fail$sampleID] <- "removed"

total_summed_tic_removed <- total_summed_tic %>% filter(grepl("removed", removed))
total_summed_tic_pass <- total_summed_tic %>% filter(grepl("pass_qc", removed))

# create a plate list ID
plate_number <- unique(plateID) %>% substr(14,14) %>% unique()
plateIDx <- plateIDx <- lapply(unique(plateID), function(FUNC_plateID){
  #browser()
  grep(FUNC_plateID, total_summed_tic$plateID)[1]}) %>% unlist()

#set y axis limits
if(tic_cut_off_lower < min(total_summed_tic$summed_TIC)){
  y_limit_lower <- log(tic_cut_off_lower-(tic_cut_off_lower/100*25))
}
if(tic_cut_off_lower > min(total_summed_tic$summed_TIC)){
  y_limit_lower <- log(min(total_summed_tic$summed_TIC)-(min(total_summed_tic$summed_TIC)/100*25))
}

y_limit_upper <- log(max(total_summed_tic$summed_TIC)+(max(total_summed_tic$summed_TIC)/100*25))



# create a layout list of extra lines to add
p_threshold_lines <- list(list(type='line', x0= min(total_summed_tic$sample_idx), x1= (max(total_summed_tic$sample_idx)+10), y0=log(tic_cut_off_lower), y1=log(tic_cut_off_lower),
                               line=list(dash='dot', width=3, color = '#FF0000')),
                          list(type='line', x0= min(total_summed_tic$sample_idx), x1= (max(total_summed_tic$sample_idx)+10), y0=log(median_summed_tic), y1=log(median_summed_tic),
                               line=list(dash='dot', width=3, color = '#000000'))
)

p_plate_list <- lapply(plateIDx[2:length(plateIDx)], function(FUNC_P_PLATE_LIST){
  list(type='line', x0 = FUNC_P_PLATE_LIST, x1= FUNC_P_PLATE_LIST, y0=y_limit_lower-log(median_summed_tic)*0.1, y1=y_limit_upper,
       line=list(dash='dot', width=2, color = '#808080'))
})

#only add plate lines if multiple plates exist
if(is.na(plateIDx)){
  p_plot_lines <- p_threshold_lines
}

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
  range = c(0, max(total_summed_tic$sample_idx)+1),
  title = "Sample index"
)

y_axis_settings <- list(
  zeroline = FALSE,
  showline = TRUE,
  linecolor = toRGB("black"),
  linewidth = 2,
  showgrid = TRUE,
  range = c(y_limit_lower-log(median_summed_tic)*0.1, 
            y_limit_upper),
  title = "Lipid total ion count (Log)"
)

p <- plot_ly(
  type = "scatter", mode = "markers", data = total_summed_tic_pass, x = ~sample_idx, y = ~LOG_summed_TIC, text = ~sampleID, color = ~removed, colors = c('#1E90FF', '#FF0000'), 
  marker = list(size = 7, color = '#1E90FF', opacity = 0.5,
                line = list(color = '#000000',width = 1))
) %>% 
  add_trace(type = "scatter", data = total_summed_tic_removed, x = ~sample_idx, y = ~LOG_summed_TIC, text = ~sampleID, color = ~removed, 
            marker = list(size = 8, color = '#FF0000')
  ) %>%
  layout(xaxis = x_axis_settings,
         yaxis = y_axis_settings
         ) %>%
  layout(shapes=p_plot_lines)

#create html widget and display it in the users internet browser
tic_check_p <- p

saveWidget(tic_check_p, file = paste(project_dir_html, "/", project_name, "_", user_name, "_TIC_check_plot.html", sep=""))# save plotly widget
browseURL(paste(project_dir_html, "/", project_name, "_", user_name, "_TIC_check_plot.html", sep="")) #open plotly widget in internet browser

#tic_qc_fail - ask the user if they wish to continue or change the threshold
tic_check_status <- "blank"
while(tic_check_status != "continue" & tic_check_status != "change"){
  tic_check_status <- dlgInput(paste(nrow(tic_qc_fail),  "samples FAILED the SIL QC check.  continue or change the exclusion threshold?"), "continue/change")$res
}
}

#tic_qc_fail - ask the user if they wish to remove all/none/samples/LTR which failed the QC check
temp_answer <- "blank"
while(temp_answer != "all" & temp_answer != "none" & temp_answer != "samples" & temp_answer != "LTR"){
  temp_answer <- dlgInput(paste("of the ", nrow(tic_qc_fail), "FAILED samples.  ",  nrow(tic_qc_fail_ltr),"  were LTRs.  Do you want to remove failed samples?"), "all/none/samples/LTR")$res
  if(temp_answer == "all"){individual_lipid_data_sil_tic_filtered <- individual_lipid_data_sil_filtered %>% filter(!sampleID %in% tic_qc_fail$sampleID)}
  if(temp_answer == "samples"){individual_lipid_data_sil_tic_filtered <- individual_lipid_data_sil_filtered %>% filter(!sampleID %in% tic_qc_fail_samples$sampleID)}
  if(temp_answer == "LTR"){individual_lipid_data_sil_tic_filtered <- individual_lipid_data_sil_filtered %>% filter(!sampleID %in% tic_qc_fail_ltr$sampleID)}
  if(temp_answer == "none"){individual_lipid_data_sil_tic_filtered <- individual_lipid_data_sil_filtered}
}

#tidy up environment
 # remove_list <- c("idx_line", "mad_summed_tic", "median_summed_tic", "plate_number", "tic_check_status", "tic_cut_off_lower","p", "total_summed_tic", "remove_list")
 # rm(list = remove_list)
