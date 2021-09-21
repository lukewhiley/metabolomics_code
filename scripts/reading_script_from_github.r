# IN USE AND WORKING

library(httr)

# Source R script from Github
script <-
  GET(
    url = "https://raw.githubusercontent.com/lukewhiley/metabolomics_code/main/lgw_boxplots.r",
    accept("application/vnd.github.v3.raw")
  ) %>% content(as = "text")

eval(parse(text = script), envir = .GlobalEnv)

test_plot_data <-  as_tibble(cbind(runif(100)*100, factor(rep(seq(1:2),50)))) 
colnames(test_plot_data) <- c("concentration", "class")

lgw_boxplot(test_plot_data)
   