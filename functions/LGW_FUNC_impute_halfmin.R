# Function to impute missing values and NA values. Assumes data missing due to <LLOQ or <LLOD so imputation method = halfmin (value[metabolite]/2)

# x = tibble with numeric metabolite data, make sure that there is no character vectors or annotation data

#impute function (method = min/2)
lgw_impute <- function(x){
    map(.x = x, .f = ~ (min(.x[.x > 0], na.rm = TRUE))/2) %>%
    #use replace_na to replace NAs with min/2 value
    replace_na(
      data = x %>% na_if(x = ., y = 0), #note - replace zeros with NA to make compatible with replace_na()
      replace = .) #note - replace with list of min/2 values generated from map function in pipe (.)
}
