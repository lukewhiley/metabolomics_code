###################################
###### LTR checks #########
##################################

# now switching to LTR comparisons across the dataset

dlg_message("Let's look at LTRs! :-)", type = 'ok')

qc_data <- individual_lipid_data_sil_tic_filtered %>% filter(!grepl("conditioning", sampleID)) %>% filter(grepl("LTR", sampleID))

area_sum <- apply(qc_data, 1, function(peakArea){
  #browser()
  func_data <- peakArea %>% as_tibble(rownames = "lipid") %>% filter(!grepl("sampleID", lipid)) %>% filter(grepl("SIL", lipid)) %>% select(value) %>% sapply(as.numeric) %>% colSums(na.rm = TRUE)
}) %>% as_tibble() %>% add_column(qc_data$sampleID, qc_data$plate_id, .before = 1) %>% 
  rename("sampleID" = "qc_data$sampleID", "plate_id" = "qc_data$plate_id", "summed_area" = "value")

temp_answer <- dlgInput("What do you wish to set for the fail cut off filter.  x number of standard deviations from the mean", "e.g.   x = 2")$res
while(is.na(as.numeric(temp_answer))){
  temp_answer <- dlgInput("You did not enter a numeric value.  What do you wish to set for the fail cut off filter.  x number of standard deviations from the mean", "e.g.   x = 2")$res
}

area_mean <- mean(area_sum$summed_area)
(sd(area_sum$summed_area)*100)/area_mean
upper_cutoff <- area_mean + temp_answer*sd(area_sum$summed_area)
lower_cutoff <- area_mean - temp_answer*sd(area_sum$summed_area)

area_sum$normalisation_factor <- area_mean/area_sum$summed_area
area_sum$sample_idx <- c(1:nrow(area_sum))
#area_sum$summed_area <- area_sum$summed_area * area_sum$normalisation_factor
area_sum$LOG_summed_TIC <- log(area_sum$summed_area)

#visualise
plot_limits <- c(area_mean/2, area_mean*2) %>% log()

p <- plot_ly(type = "scatter", mode = "markers", area_sum, x = ~sample_idx, y = ~LOG_summed_TIC, text = ~sampleID, colors = c("lightblue3")) %>%
  layout(yaxis = list(range = c(plot_limits)))

#add vertical lines to represent plate boundaries - only runs if multiple plates exist using ANPC LIMS fine name style
plateid <- str_extract(area_sum$sampleID, "PLIP.*")
if(100/length(plateid)*length(which(is.na(plateid))) < 50){
  plateid <- substr(plateid, 0,15)
  plate_id <- paste(plateid, sub(".*\\_", "", area_sum$sampleID), sep="_")
  plate_number <- unique(plate_id) %>% substr(14,14) %>% unique()
  plate_idx <- lapply(unique(plateid), function(plateID){grep(plateID, area_sum$sampleID)[1]}) %>% unlist()
  
  for (idx_line in 2:length(plate_idx)){
    p <- p %>% add_segments(x = plate_idx[idx_line], xend = plate_idx[idx_line], y = plot_limits[1], yend = plot_limits[2], line = list(color = "grey", dash = "dash"), showlegend = FALSE)
  }
}

p <- p %>% add_segments(y = log(area_mean), yend = log(area_mean), x = 0, xend = nrow(area_sum), line = list(color = "black", dash = "dash"), showlegend = FALSE)
p <- p %>% add_segments(y = log(upper_cutoff), yend = log(upper_cutoff), x = 0, xend = nrow(area_sum), line = list(color = "red", dash = "dot"), showlegend = FALSE)
p <- p %>% add_segments(y = log(lower_cutoff), yend = log(lower_cutoff), x = 0, xend = nrow(area_sum), line = list(color = "red", dash = "dot"), showlegend = FALSE)


LTR_sil_check_p <- p
if(!dir.exists(paste(project_dir, "/html_files", sep=""))){dir.create(paste(project_dir, "/html_files", sep=""))} # create a new directory to store html widgets
saveWidget(LTR_sil_check_p, file = paste(project_dir, "/html_files/",project_name, "_", user_name, "_LTR_SIL_check_plot.html", sep=""))# save plotly widget
browseURL(paste(project_dir, "/html_files/",project_name, "_", user_name, "_LTR_SIL_check_plot.html", sep="")) #open plotly widget in internet browser


sil_area_sum <- area_sum