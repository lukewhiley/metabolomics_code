library(FELLA)

g <- buildGraphFromKEGGREST(
  organism = "hsa", 
  filter.path = "hsa01100"
)
buildDataFromGraph(g)
