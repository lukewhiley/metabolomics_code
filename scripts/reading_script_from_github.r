# IN USE AND WORKING

library(httr)

# Source R script from Github
script <-
  GET(
    url = "https://raw.githubusercontent.com/lukewhiley/metabolomics_code/main/lgw_boxplots.r",
    authenticate("luke.whiley@murdoch.edu.au", "cda3533bb0a8248333824bd91029b48c367164d4")  ,   # Instead of PAT, could use password
    accept("application/vnd.github.v3.raw")
  ) %>% content(as = "text")

eval(parse(text = script), envir = .GlobalEnv)

test_plot_data <-  as_tibble(cbind(runif(100)*100, factor(rep(seq(1:2),50)))) 
colnames(test_plot_data) <- c("concentration", "class")

lgw_boxplot(test_plot_data)
   