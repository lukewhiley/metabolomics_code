# 
library(tidyverse)
library(janitor)
data_a <- read_csv(file = file.choose()) %>% clean_names()

metabolite_list <- names(data_a)[2:ncol(data_a)]

MW_resultList <- lapply(metabolite_list, function(FUNC_STATS){
  #browser()
  func_data <- data_a %>% select(treatment_allocation_blind, all_of(FUNC_STATS))
  colnames(func_data) <- c("treatment", "metabolite")
  mw_res <- wilcox.test(metabolite ~ as.factor(treatment),
                        data = func_data,
                        exact=FALSE
                        )
  mw_res$p.value
}) %>% unlist() %>% c() %>% as_tibble() %>% add_column(metabolite_list, .before = 1) %>% arrange(value)

colnames(MW_resultList) <- c("metabolite", "pvalue")

MW_resultList

MW_boxplotList <- lapply(metabolite_list, function(FUNC_STATS){
  browser()
  func_data <- data_a %>% select(treatment_allocation_blind, all_of(FUNC_STATS))
  colnames(func_data) <- c("treatment", "metabolite")
  boxplot(metabolite ~ as.factor(treatment),
         data = func_data
         )
  title(paste(FUNC_STATS))
  
})



