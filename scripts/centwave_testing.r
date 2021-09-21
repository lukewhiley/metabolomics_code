cwp <- CentWaveParam(ppm = 1000,
                     peakwidth = c(2, 20),
                     snthresh = 2, 
                     prefilter = c(3, 10), 
                     mzCenterFun = "wMean",
                     integrate=1, 
                     mzdiff = -0.001,
                     fitgauss = FALSE,
                     noise = 0, 
                     verboseColumns =FALSE,
                     #roiList = 0,
                     firstBaselineCheck = TRUE,
                     #roiScales=0
                     )
  
res <- findChromPeaks(test_spectra, param = cwp)



mfp <- MatchedFilterParam(binSize = 0.5,
                          impute = "none",
                          baseValue = NULL,
                          distance = NULL,
                          fwhm = 10,
                          sigma = 12.7399,
                          max =5,
                          )