# code for rOPLS

# create rOPLS model

# for PCA

pca_model <- ropls::opls(
  x = , # data matrix (cols = features, rows = samples)
  y = NULL, #do not need y for PCA
  predI = 2, #how many components
  orthoI = NA,  # not required for PCA
  algoC = "nipals", #choice of algorithm
  crossvalI = 3, #cross validation segments
  permI = 20, #default is 20
  log10L = TRUE, #log data? TRUE/FALSE
  scaleC = "pareto", #scale data
  fig.pdfC = "none", # display defalt figures none/interactive/ or file xport path
  plotSubC = NA, #subtitle
  subset = NULL, #indices to creat training/test data
  info.txtC = "none",  #file name to export results or "none"
)

# extract data for plotting in ggplot2
# data frame is in same order as x in pca model, may need to add sample info columns for colour control etc
scores_data <- pca_model@scoreMN 
#then perform ggplot

#loadings data
loadings_data <- pca_model@loadingMN


#for OPLS-DA
opls_model <- ropls::opls(
  x = , # data matrix (cols = features, rows = samples)
  y = , #vector of y values (e.g. control or PD) 
  predI = 1, #how many components #1 here for OPLS-DA
  orthoI = 1,  # orthogonal components
  algoC = "nipals", #choice of algorithm
  crossvalI = 3, #cross validation segments
  permI = 20, #default is 20
  log10L = TRUE, #log data? TRUE/FALSE
  scaleC = "pareto", #scale data
  fig.pdfC = "none", # display defalt figures none/interactive/ or file xport path
  plotSubC = NA, #subtitle
  subset = NULL, #indices to creat training/test data
  info.txtC = "none",  #file name to export results or "none"
)

# extract data for plotting in ggplot2
# data frame is in same order as x in pca model, may need to add sample info columns for colour control etc
scores_data <- opls_model@scoreMN 
#then perform ggplot

#loadings data
loadings_data <- opls_model@loadingMN