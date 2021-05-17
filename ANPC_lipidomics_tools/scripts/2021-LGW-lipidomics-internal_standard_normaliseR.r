################################################
###### Internal standard normalization #########
#################################################

# this section normalizes each SRM lipid target with the appropriate internal standard
# Requires a template guide with internal standard transition for each target lipid SRM transition

dlg_message("Time for normalization using the internal standards :-)", type = 'ok')
dlg_message("REQUIRES - a reference csv file listing each lipid target and the assigned internal standard. Column headings required are 'Precursor Name' containing the lipid target  and 'Note' containing the IS name", type = 'ok')

#import transition report 3
filtered_data <- individual_lipid_data_sil_tic_intensity_filtered %>% filter(!grepl("conditioning", sampleID))
filtered_data[is.na(filtered_data)] <- 0

dlg_message("Please select this template file now.", type = 'ok')

temp_answer <- "change"
while(temp_answer == "change"){
sil_target_list <- read_csv(file = file.choose(.)) %>% clean_names
lipid_data <- filtered_data %>% select(-sampleID, - plate_id) %>% select(!contains("SIL"))
sil_data <- filtered_data %>% select(-sampleID, - plate_id) %>% select(contains("SIL"))

# this section checks each of the SIL IS used in the target list template in the LTRS. It evaluates if:
##  a: is the internal standard present in the LTR samples? Some batches of IS do not contain every IS availible. This alos prevents user error if the IS batch has not been made correctly.
##  b: the RSD of the internal standard signal intensity in LTRs. If the SIL IS is suitable for use in the dataset its signal should be stable in the LTRs. (raises a warning for >30%)

dlg_message("Checks to see if all internal standards are present in the SIL internal standard mix", type = 'ok')

sil_data_check <- individual_lipid_data_sil_tic_filtered %>% select(sampleID, plate_id, contains("SIL")) %>% filter(grepl("LTR", sampleID))

sil_list <- sil_target_list %>% filter(grepl("SIL", note)) %>% select(note) %>% unique() 

sil_sum <- lapply(sil_list$note, function(FUNC_SIL){
  #browser()
  temp_func_data_sum <- sil_data_check %>% select(all_of(FUNC_SIL)) %>% as.matrix() %>% sum() %>% log()
}) %>% c() %>% unlist() %>% as_tibble() %>% rename(SIL_SUM = value) %>% add_column(sil_list, .before = 1) %>% arrange(SIL_SUM)

sil_sum_q1 <- quantile(sil_sum$SIL_SUM, 0.25) %>% as.numeric()
inter_quantile_range <- as.numeric(quantile(sil_sum$SIL_SUM, 0.75)) - as.numeric(quantile(sil_sum$SIL_SUM, 0.25))
sil_sum_lower_threshold <- sil_sum_q1 - inter_quantile_range

#create a list of IS that fail the test
sil_list_warning <- sil_sum$note[which(sil_sum$SIL_SUM < sil_sum_lower_threshold)]

sil_rsd <- lapply(sil_list$note, function(FUNC_SIL){
  #browser()
  temp_func_data <- sil_data_check %>% select(all_of(FUNC_SIL))
  temp_func_data_mean <- temp_func_data %>% as.matrix() %>% mean()
  temp_func_data_sd <- temp_func_data %>% as.matrix() %>% sd()
  temp_func_data_rsd <- (temp_func_data_sd/temp_func_data_mean)*100
  temp_func_data_rsd
}) %>% c() %>% unlist() %>% as_tibble() %>% rename(SIL_RSD = value) %>% add_column(sil_list, .before = 1)

sil_list_warning_2 <- sil_rsd %>% filter(SIL_RSD > 30) 

dlg_message(paste( "############################################", 
                   "Warning 1! These internal standards have a very low signal and may not be present in the mixture:",
                   paste(sil_list_warning, collapse = ";     ."),
                   "############################################", 
                   "Warning 2! These internal standards have a %RSD >30%: ",
                   paste(sil_list_warning_2$note, collapse = ";     ."),
                   "############################################",
                   "double check skyline!"), 
            type = 'ok')


sil_list_warning <- c(sil_list_warning, unlist(sil_list_warning_2$note))

temp_answer <- dlgInput("Do you wish to continue or use different internal standards?", "continue/change")$res

#add in a check in case the user enters the incorrect entry. It must be "continue" or "change" to continue
while(temp_answer != "continue" & temp_answer!="change"){
  temp_answer <- dlgInput("Do you wish to continue or use different internal standards?", "continue/change")$res
}

# if change has been selected the user now has an opportunity to change the internal standard import file template
if(temp_answer == "change"){
  dlg_message("OK - please edit and select a new internal standard reference file now")
}
}

# now multiply by the IS concentration to create a concentration factor
dlg_message("We also need to calculate lipid concentrations from the internal standard. Please select the concentration template csv - NOTE: use the correct lot for your analysis", type = 'ok')

#read in SIL-internal standard concentration file. aim to put on github to centralise
sil_concentrations <- read_csv(file = file.choose(.)) %>% clean_names

sil_batch <- "blank"
while(is.na(as.numeric(sil_batch))){
  sil_batch <- dlgInput("What batch did you use?", "101/102/103")$res
}


# this apply function creates:
# 1. a response ratio by dividing the signal area for each target lipid by the peak area from the appropriate SIL IS metabolite. 
# 2. a final estimated concentration using the pre-defined internal standard as a single point calibration


ratio_data <- apply(as_tibble(colnames(lipid_data)), 1, function(FUNC_IS_RATIO){
  browser()
  # step 1 - create a ratio between lipid target and the appropriate internal standard as pre-defined in the reference file
  func_data <- lipid_data %>% select(all_of(FUNC_IS_RATIO))
  sil_to_use <- sil_target_list$note[which(sil_target_list$precursor_name==FUNC_IS_RATIO)]
  func_data_sil <- sil_data %>% select(all_of(sil_to_use))
  normalised_data <- func_data/func_data_sil
  
  #step 2 - select concentration factor from csv template and multiply normalised data by concentration factor
  func_concentration_factor <- sil_concentrations %>% filter(sil_name == sil_to_use) %>% select(concentration_factor) %>% as.numeric()
  concentration_data <- normalised_data*func_concentration_factor
  concentration_data
  
}) %>% bind_cols() %>% add_column(filtered_data$sampleID, filtered_data$plate_id, .before = 1)

colnames(ratio_data) <- c("sampleID", "plateID", colnames(lipid_data))

ltr_rsd <- apply(as_tibble(colnames(lipid_data)), 1, function(LTR_RSD){
  #browser()
  func_data <- ratio_data %>% filter(grepl("LTR", sampleID)) %>% select(all_of(LTR_RSD))
  (sd(func_data$value)*100)/mean(func_data$value)
}) %>% as_tibble() %>% add_column(colnames(lipid_data), .before = 1)

colnames(ltr_rsd) <- c("lipid", "RSD")
dlg_message(paste(nrow(ltr_rsd)-length(which(ltr_rsd$RSD < 30)), " lipid targets had a LTR RSD of > 30% so were removed from the dataset", sep=""))
dlg_message(paste("number of feature ratios with with an LTR RSD of <30% =", length(which(ltr_rsd$RSD < 30))), type = 'ok')
dlg_message(paste("number of feature ratios with with an LTR RSD of <20% =", length(which(ltr_rsd$RSD < 20))), type = 'ok')
dlg_message(paste("number of feature ratios with with an LTR RSD of <15% =", length(which(ltr_rsd$RSD < 15))), type = 'ok')
dlg_message(paste("number of feature ratios with with an LTR RSD of <10% =", length(which(ltr_rsd$RSD < 10))), type = 'ok')

lipid_keep_list <- ltr_rsd %>% filter(RSD < 30)

final_dataset <- ratio_data %>% select(sampleID, plateID, all_of(lipid_keep_list$lipid))
final_dataset[is.na(final_dataset)] <- 0






# visualisation of normalized data
# first - produce a plot of all normalized features to see if there are any overall trends in the data

total_summed_ratio <- apply(final_dataset %>% select(sampleID), 1, function(summedTIC){
  #browser()
  temp_data <- final_dataset %>% filter(sampleID == summedTIC) %>% select(-sampleID, -plateID) %>% rowSums(na.rm = TRUE)
}) %>% c() %>% as_tibble() %>%  add_column(final_dataset$sampleID, .before = 1) %>% 
  rename(summed_TIC = value, sampleID = "final_dataset$sampleID")
total_summed_ratio$sample_idx <- c(1:nrow(total_summed_ratio))

sd(total_summed_ratio$summed_TIC*100)/mean(total_summed_ratio$summed_TIC)

total_summed_ratio$sample <- "sample"
total_summed_ratio$sample[grep("LTR", total_summed_ratio$sampleID)] <- "LTR"

# create a plotly plot to visualize
total_summed_ratio$log_summed_TIC <- log(total_summed_ratio$summed_TIC+1)

total_summed_ratio_samples <- total_summed_ratio %>% filter(!grepl("LTR", sampleID))
total_summed_ratio_LTR <- total_summed_ratio %>% filter(grepl("LTR", sampleID))

# create a plate list ID
plate_number <- unique(plate_id) %>% substr(14,14) %>% unique()
plate_idx <- lapply(unique(plateid), function(plateID){grep(plateID, total_summed_ratio$sampleID)[1]}) %>% unlist()

# create a layout list of extra lines to add
p_plate_list <- lapply(plate_idx[2:length(plate_idx)], function(FUNC_P_PLATE_LIST){
  list(type='line', x0 = FUNC_P_PLATE_LIST, x1= FUNC_P_PLATE_LIST, 
       y0=log(min(total_summed_ratio_samples$summed_TIC)-(min(total_summed_ratio_samples$summed_TIC)/100*50)), 
       y1=log(max(total_summed_ratio_samples$summed_TIC)+(max(total_summed_ratio_samples$summed_TIC)/100*25)),
       line=list(dash='dot', width=2, color = '#808080'))
})

#only add plate lines if multiple plates exist
if(is.na(plate_idx)){
  p_plot_lines <- NULL
}

if(length(plate_idx) == 1){
  p_plot_lines <- NULL
}

if(length(plate_idx) > 1){
  p_plot_lines <- c(p_plate_list)
} 

#create a list of axis settings for plot_ly
x_axis_settings <- list(
  zeroline = FALSE,
  showline = TRUE,
  linecolor = toRGB("black"),
  linewidth = 2,
  showgrid = FALSE,
  range = c(0, max(total_summed_ratio_samples$sample_idx)),
  title = "Sample index"
)

y_axis_settings <- list(
  zeroline = FALSE,
  showline = TRUE,
  linecolor = toRGB("black"),
  linewidth = 2,
  showgrid = TRUE,
  range = c(log(min(total_summed_ratio_samples$summed_TIC)-(min(total_summed_ratio_samples$summed_TIC)/100*50)), 
            log(max(total_summed_ratio_samples$summed_TIC)+(max(total_summed_ratio_samples$summed_TIC)/100*25))),
  title = "Summed lipid target response ratio (Log)"
)

p <- plot_ly(
  type = "scatter", mode = "markers", data = total_summed_ratio_samples, x = ~sample_idx, y = ~log_summed_TIC, text = ~sampleID, color = ~sample, colors = c('#1E90FF', '#FF0000'), 
  marker = list(size = 7, color = '#1E90FF', opacity = 0.5,
                line = list(color = '#000000',width = 1))
) %>% 
  add_trace(type = "scatter", data = total_summed_ratio_LTR, x = ~sample_idx, y = ~log_summed_TIC, text = ~sampleID, color = ~sample, 
            marker = list(size = 8, color = '#FF0000')
  ) %>%
  layout(xaxis = x_axis_settings,
         yaxis = y_axis_settings
  ) %>%
  layout(shapes=p_plot_lines)

normalized_check_p <- p

#create html widget and display it in the users internet browser
if(!dir.exists(paste(project_dir, "/html_files", sep=""))){dir.create(paste(project_dir, "/html_files", sep=""))} # create a new directory to store html widgets
saveWidget(normalized_check_p, file = paste(project_dir, "/html_files/",project_name, "_", user_name, "_normalized_check_plot.html", sep=""))# save plotly widget
browseURL(paste(project_dir, "/html_files/",project_name, "_", user_name, "_normalized_check_plot.html", sep="")) #open plotly widget in internet browser

dlg_message("Check plot for summed all normilized features. Press OK to continue", type = 'ok')

# second - produce a plot of normalized features, summed by class to see if there are any overall trends in the lipid class data

dlg_message("Now we are going to look at the data summed by lipid class", type = 'ok')

final_individual_lipid_data <- final_dataset

final_class_lipid_data <- create_lipid_class_data_summed(final_individual_lipid_data)

lipid_class_list <- final_individual_lipid_data %>% select(contains("(")) %>% colnames() 
lipid_class_list <- sub("\\(.*", "", lipid_class_list) %>% unique()
lipid_class_list <- lipid_class_list[!grepl("sampleID", lipid_class_list)] %>% as_tibble()

#add ltr TRUE/FALSE column

final_class_lipid_data$is_ltr <- "sample"
final_class_lipid_data$is_ltr[grep("LTR", final_class_lipid_data$sampleID)] <- "LTR"

plotlist <- apply(lipid_class_list %>% select(value), 1, function(lipidClass){
  #browser()
  plot_data <- final_class_lipid_data %>% select("sampleID", "is_ltr", all_of(lipidClass)) %>% rename(ms_response = value) 
  plate_id <- str_extract(plot_data$sampleID, "PLIP.*")
  plate_id <- substr(plate_id, 0,15)
  plot_data$sample_index <- paste(plate_id, sub(".*\\_", "", plot_data$sampleID), sep="_")
  plot_data <- plot_data %>% arrange(sample_index)
  plot_data$sample_idx <- 1:nrow(plot_data)
  
  plot_data_ltr <- plot_data %>% filter(is_ltr == "LTR")
  plot_data <- plot_data %>% filter(is_ltr == "sample")
  
  plate_idx <- lapply(unique(plate_id), function(plateID){
    grep(plateID, plot_data$sampleID)[1]
  }) %>% unlist
  
  # create a plate list ID
  plate_number <- unique(plate_id) %>% substr(14,14) %>% unique()
  plate_idx <- lapply(unique(plateid), function(plateID){grep(plateID, total_summed_ratio$sampleID)[1]}) %>% unlist()
  
  # create a layout list of extra lines to add
  p_plate_list <- lapply(plate_idx[2:length(plate_idx)], function(FUNC_P_PLATE_LIST){
    list(type='line', x0 = FUNC_P_PLATE_LIST, x1= FUNC_P_PLATE_LIST, 
         y0=(log(min(plot_data$ms_response)+1)-(log(min(plot_data$ms_response)+1)/100*50)), 
         y1=(log(max(plot_data$ms_response)+1)+(log(max(plot_data$ms_response)+1)/100*25)),
         line=list(dash='dot', width=2, color = '#808080'))
  })
  
  #only add plate lines if multiple plates exist
  if(is.na(plate_idx)){
    p_plot_lines <- NULL
  }
  
  if(length(plate_idx) == 1){
    p_plot_lines <- NULL
  }
  
  if(length(plate_idx) > 1){
    p_plot_lines <- c(p_plate_list)
  } 
  
  
  
  #create a list of axis settings for plot_ly
  x_axis_settings <- list(
    zeroline = FALSE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = FALSE,
    range = c(0, max(plot_data$sample_idx)),
    title = paste(lipidClass)
    )
  
  y_axis_settings <- list(
    zeroline = FALSE,
    showline = TRUE,
    linecolor = toRGB("black"),
    linewidth = 2,
    showgrid = TRUE,
    range = c(log(min(plot_data$ms_response)+1)-(log(min(plot_data$ms_response)+1)/100*50), 
              log(max(plot_data$ms_response)+1)+(log(max(plot_data$ms_response)+1)/100*25)
              ),
    title = "Summed lipid target/internal standard ratio (Log)"
  )

  
  p <- plot_ly(
    type = "scatter", mode = "markers",  colors = c('#1E90FF', '#FF0000'), data = plot_data, x = ~sample_idx, y = ~log(ms_response+1), text = ~sampleID, color = ~is_ltr, 
    marker = list(size = 5, color = '#1E90FF', opacity = 0.5,
                  line = list(color = '#000000',width = 1)),
    showlegend = FALSE
  ) %>% 
    add_trace(type = "scatter", data = plot_data_ltr, x = ~sample_idx, y = ~log(ms_response+1), text = ~sampleID, color = ~is_ltr, 
              marker = list(size = 6, color = '#FF0000'),
              showlegend = FALSE
    ) %>%
    layout(xaxis = x_axis_settings,
           yaxis = y_axis_settings
           ) %>%
    layout(shapes=p_plot_lines)
  p
})

normalized_check_class_p <- subplot(plotlist, nrows = 4, titleX = TRUE, margin = c(0.01,0.01,0.05,0.05))

if(!dir.exists(paste(project_dir, "/html_files", sep=""))){dir.create(paste(project_dir, "/html_files", sep=""))} # create a new directory to store html widgets
saveWidget(normalized_check_class_p, file = paste(project_dir, "/html_files/",project_name, "_", user_name, "_normalized_check_class_plot.html", sep=""))# save plotly widget
browseURL(paste(project_dir, "/html_files/",project_name, "_", user_name, "_normalized_check_class_plot.html", sep="")) #open plot_ly widget in internet browser

dlg_message("Check plot for summed lipid class normilized features. Press OK to continue", type = 'ok')

