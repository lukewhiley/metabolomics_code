# kruskal walis and dunn test function

lgw_kw_dunn <- function(FUNC_data,
                        FUNC_metabolite_list,
                        FUNC_HEADER_class){
  #require(rstatix)
  
  #browser()
  
  kw_result <- FUNC_metabolite_list %>%
    as_tibble() %>%
    rename(lipid = value)
  kw_result$pvalue <- NA
  
  dunn_out <- NULL
  
  for(idx_lipid in FUNC_metabolite_list){
    temp_data <- FUNC_data %>% select(all_of(FUNC_HEADER_class), all_of(idx_lipid))
    colnames(temp_data) <- c("class", "data")
    temp_kw_result <- kruskal.test(data ~ class, data = temp_data)
    temp_dunn_result <- rstatix::dunn_test(data = temp_data, 
                                           formula = data ~ class,
                                           p.adjust.method = "BH")
    dunn_out <- rbind(dunn_out, temp_dunn_result$p.adj)
    
    kw_result$pvalue[which(kw_result$lipid == idx_lipid)] <- temp_kw_result$p.value
  }
  
  #browser()
  kw_result$qvalue <- kw_result$pvalue %>% p.adjust(method = "BH")
  
  #browser()
  
  dunn_out <- as_tibble(dunn_out)
  colnames(dunn_out) <- paste(temp_dunn_result$group1,
                              "_",
                              temp_dunn_result$group2)
  
  kw_result_neg <- cbind(kw_result,
                         dunn_out) %>%
    arrange(pvalue) #%>%
    #filter(pvalue < 0.05)
  
}
  
  