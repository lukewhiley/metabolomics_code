#' ---
#' title: "Target Lipid ExploreR QC Report"
#' author: ANPC
#' output: html_document
#' 
#' ---

#' Thank you for using LGW SkylineR and Lipid ExploreR. Your lipidomics quality control evaluation report will now be produced.
#'
#'
#' ### 1. Stable isotope label (SIL) internal standard (IS) quality control check
#'
#' The first QC check summed the peak areas from all stable isotope labeled (SIL) internal standards (IS) to identify those samples that either:  
#'
#' * contained low/no internal standard. Suggesting a preparation error where low volume of the SIL IS had been added to the well 
#' * contained excessive internal standard. Suggesting a preparation error where a droplet may have entered the well increasing the volume of the SIL IS
#' 
#' 
#' 
#' #### Internal standard QC plot
#' 
#' This plot is the result of summing all the internal standards into a single signal for each sample. All samples are plotted including LTR samples.
#' 
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=4
sil_check_p
#'
#'
##+ echo=FALSE, message=FALSE, fig.width=10, fig.height=4
print(paste(nrow(sil_qc_fail_samples), " study samples and ", nrow(sil_qc_fail_ltr), " LTR samples were removed from the dataset during this check", sep = ""))
#'
#' The following samples were removed from the dataset at this check point
#' 
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=4 
knitr::kable(sil_qc_fail)
#'
#'
#'
#' ### 2. Total ion count (TIC) quality control check
#' 
#' The second QC check summed the peak area of every lipid target to produce a total ion count signal for each sample. The check is used to identify those samples that:  
#' 
#' * contained low/no signal. Suggesting a preparation error where low volume of the sample had been added to the well  
#' 
#' 
#' 
#' #### Total ion count QC plot
#' 
#' This plot is the result of summing all the lipid target peak areas to form a single total ion count signal. All samples are plotted including LTR samples.
#' 
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=4
tic_check_p 
#'
#'
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=4
print(paste(nrow(tic_qc_fail_samples), " study samples and ", nrow(tic_qc_fail_ltr), " LTR samples were removed from the dataset during this check", sep = ""))
#'
#' The following samples were removed from the dataset at this check point
#'  
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=4 
knitr::kable(tic_qc_fail)
#'
#'
#'
#'
#' ### 3. Target lipid feature quality control check
#' The QC checks now switched to the lipid mass spectrometry targets. Lipids were filtered out due to a user defined signal intensity and appearance frequency.
#' For example the default threshold for inclusion is to only keep features with a signal intensity (height) greater than 2500 counts in over 50% of LTR samples
#'
#'
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=4

if(intensity_threshold_ltr == "LTR"){
  print(paste("A lipid feature was kept in the dataset if it was present in ", intensity_threshold_percentage, " % of ",  intensity_threshold_ltr, " samples with a signal intensity (peak height) greater than ", intensity_threshold, " counts", sep = ""))
}

if(intensity_threshold_ltr == "samples"){
  print(paste("A lipid feature was kept in the dataset if it was present in ", intensity_threshold_percentage, " % of study",  intensity_threshold_ltr, " (not LTR) with a signal intensity (peak height) greater than ", intensity_threshold, " counts", sep = ""))
}

if(intensity_threshold_ltr == "both"){
  print(paste("A lipid feature was kept in the dataset if it was present in ", intensity_threshold_percentage, " % of both study samples and LTR samples with a signal intensity (peak height) greater than ", intensity_threshold, " counts", sep = ""))
}

print(paste(length(lipid_intensity_filter_fail),"lipid targets were removed from the dataset as they did not meet these criteria"))
#'
#'
#' ### 4. Creation of response values and response ratio QC check
#' Following the creation of a response ratio with an appropriate internal standard (as pre-defined by the user) features are removed at this QC checkpoint if the percentage relative standard deviation (%RSD) was greater than 30% in the LTR QC samples
#'
#'  
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=4

print(paste(nrow(ltr_rsd_1)-length(which(ltr_rsd_1$RSD < 30)), " lipid targets had a LTR RSD of > 30% so were removed from the dataset", sep=""))
print(paste("Total number of lipid target response ratios with with an LTR RSD of <30% =", length(which(ltr_rsd_1$RSD < 30))))
print(paste("Total number of lipid target response ratios with with an LTR RSD of <20% =", length(which(ltr_rsd_1$RSD < 20))))
print(paste("Total number of lipid target response ratios with with an LTR RSD of <15% =", length(which(ltr_rsd_1$RSD < 15))))
print(paste("Total number of lipid target response ratios with with an LTR RSD of <10% =", length(which(ltr_rsd_1$RSD < 10))))

#' #### Response ratio plots
#' 
#' This plot is created by summing the response ratio from every lipid target to create a single value for each sample. The plot includes both samples and LTRs
#' 
#' 
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=4

normalized_check_p_1

#' 
#' 
#' This plot is created by summing the response ratio from every lipid target to create a single value PER LIPID CLASS for each sample. The plot includes both samples and LTRs
#' 
#' 
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=8

normalized_check_class_p_1

#'
#'
#'
#' ### 6. PCA plot to visualize final dataset variance (LTR and samples)
#' 
#' 
#' 
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=4

print(paste("For creating the QC PCA plots",  pca_scale_used_1, "scaling was used"))
print(paste("PCA plot created using the",  length(which(ltr_rsd_1$RSD < 30)), "lipids in the final dataset" ))
pca_p_1[[1]][[1]]

#print(paste("PCA plot using lipid class created by summing individual lipid response ratios grouped by lipid class"))
#pca_p_1[[2]][[1]]

#'
#'
#' ### 7. Data correction for signal drift
#' 
#' Using the LTR pooled QC samples the signals was corrected for signal drift across the run
#'  
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=4

if(signal_drift_method == "loess"){
print(paste("For correction of the signal drift a loess method (QC-RLSC) was employed using the statTarget package (bioconductor)"))
}

if(signal_drift_method == "RF"){
  print(paste("For correction of the signal drift a Random Forest method (QC-RFSC) was employed using the statTarget package (bioconductor"))
}

print(paste(nrow(ltr_rsd_2)-length(which(ltr_rsd_2$RSD < 30)), " lipid targets had a LTR RSD of > 30% so were removed from the dataset", sep=""))
print(paste("Total number of lipid target response ratios with with an LTR RSD of <30% =", length(which(ltr_rsd_2$RSD < 30))))
print(paste("Total number of lipid target response ratios with with an LTR RSD of <20% =", length(which(ltr_rsd_2$RSD < 20))))
print(paste("Total number of lipid target response ratios with with an LTR RSD of <15% =", length(which(ltr_rsd_2$RSD < 15))))
print(paste("Total number of lipid target response ratios with with an LTR RSD of <10% =", length(which(ltr_rsd_2$RSD < 10))))


#' #### Response ratio plots
#' 
#' This plot is created by summing the response ratio from every lipid target to create a single value for each sample. The plot includes both samples and LTRs
#' This plot is created post-correction for sample drift
#' 
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=4

normalized_check_p_2

#' 
#' 
#' This plot is created by summing the response ratio from every lipid target to create a single value PER LIPID CLASS for each sample. The plot includes both samples and LTRs.
#' This plot is created post-correction for sample drift
#' 
#' 
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=8

normalized_check_class_p_2

#'
#' ### 8. PCA plot to visualize final corrected dataset variance (LTR and samples)
#' 
#' 
#' 
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=4

print(paste("For creating the QC PCA plots",  pca_scale_used_2, "scaling was used"))
print(paste("PCA plot created using the",  length(which(ltr_rsd_2$RSD < 30)), "lipids in the final dataset" ))
pca_p_2[[1]][[1]]

#print(paste("PCA plot using lipid class created by summing individual lipid response ratios grouped by lipid class"))
#pca_p_2[[2]][[1]]

