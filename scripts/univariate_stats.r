#!/usr/bin/env Rscript


# Required libraries
# -----------------------------------------------------------------------------


library(tidyverse)
library(janitor)


# Function definitions
# -----------------------------------------------------------------------------


resultListFunc <- function(func_stats, data) {
  # FIXME:
  # - UNDEFVAR treatment_allocation_blind
  func_data <- data %>% select(treatment_allocation_blind, all_of(func_stats))
  colnames(func_data) <- c("treatment", "metabolite")
  mw_res <- wilcox.test(
    metabolite ~ as.factor(treatment),
    data = func_data,
    exact = FALSE)
  return(mw_res$p.value)
}


boxplotListFunc <- function(func_stats, data) {
  # FIXME:
  # - UNDEFVAR treatment_allocation_blind
  func_data <- data %>% select(treatment_allocation_blond, all_of(func_stats))
  colnames(func_data) <- c("treatment", "metabolite")
  boxplot <- boxplot(
    metabolite ~ as.factor(treatment),
    data = func_data,
    main = paste(func_stats))
  return(boxplot)
}


# Operations
# -----------------------------------------------------------------------------


# Choose a file to parse as CSV, and clean the names
data_a <- read_csv(file = file.choose()) %>% clean_names()

# Subset the data to the metabolite data
metabolite_list <- names(data_a)[2:ncol(data_a)]

# Process data and reformat as Tibble
MW_resultList <- lapply(metabolite_list, resultListFunc, data = data_a) %>%
  unlist() %>%
  c() %>%
  as_tibble() %>%
  add_column(metabolite_list, .before = 1) %>%
  arrange(value)

# Replace the column names of the result list Tibble
colnames(MW_resultList) <- c("metabolite", "pvalue")

# Produce boxplots of each group in the metabolite list
MW_boxplotList <- lapply(metabolite_list, boxplotListFunc, data = data_a)
