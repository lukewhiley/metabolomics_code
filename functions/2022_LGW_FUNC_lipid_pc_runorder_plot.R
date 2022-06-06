#ANPC PCA quality control visualisation

# FUNC_data = data containing individual lipid data - MUST CONTAIN sampleID column 
# FUNC_colour_by = how to colour the plot (e.g. sample class, is_ltr or cohort)
# FUNC_plot label = what to label the scores plot with (e.g. sampleID)
# FUNC_scaling = UV or Pareto

lgw_pc_run_plot <- function(FUNC_data, 
                    FUNC_metabolite_list, 
                    FUNC_colour_by, 
                    FUNC_plot_label, 
                    FUNC_scaling,
                    FUNC_title,
                    FUNC_project_colours,
                    #FUNC_option_invert_y,
                    #FUNC_option_invert_x,
                    FUNC_option_plot_qc
                    #FUNC_option_iqr_filter_samples,
                    #FUNC_option_iqr_filter_qc
                    ){
  require(metabom8)
  require(RColorBrewer)
  require(tidyverse)
  require(plotly)
  
  pca_output <- list()
  pca_output$plots <- list()
  
  
  #browser()
  
  title_text <- FUNC_title
  
  qc_idx <- which(FUNC_data[["sample_type"]] == "qc")
  
  if(FUNC_option_plot_qc == FALSE){
    FUNC_data <- FUNC_data %>% 
      filter(sample_type != "qc")
  }
  
  
  
  #create data matrix for PCA
  pca_x <- FUNC_data %>%  select(all_of(FUNC_metabolite_list)) %>% as.matrix()+1 
  pca_x[pca_x == 1] <- NA #remove all 0 values (above adds 1 to all values therefore anything that = 1 was a 0)
  pca_x[is.infinite(pca_x)] <- NA #remove all infinite values
  min_value <- min(pca_x, na.rm = TRUE) # find the lowest value in the matrix
  pca_x[is.na(pca_x)] <- min_value/100 # replace all NA, Inf, and 0 values with the lowest value in the matrix/100 to represent a value below limit of detection
  
  pca_x <- log(pca_x+1) #log values for plotting
  
  #create PCA model
  pca_output$pca_model <- pca(pca_x, pc=3, scale = paste(FUNC_scaling), center = TRUE)
  
  # extract score values for plotting in plot_ly
  PC1 <- as.numeric(as.matrix(pca_output$pca_model@t[,1]))
  PC2 <- as.numeric(as.matrix(pca_output$pca_model@t[,2]))
  PC3 <- as.numeric(as.matrix(pca_output$pca_model@t[,3]))

  plot_Val <- as_tibble(cbind(PC1, PC2, PC3, FUNC_data$sample_idx, FUNC_data$sample_name,  FUNC_data$sample_type, FUNC_data$sample_batch)) %>% 
    setNames(c("PC1", "PC2", "PC3", "sample_idx", "sample_name", "sample_type", "sample_batch"))
  plot_Val$sample_idx <- c(1:nrow(plot_Val)) %>% as.numeric()
  plot_Val$PC1 <- plot_Val$PC1 %>% as.numeric()
  plot_Val$PC2 <- plot_Val$PC2 %>% as.numeric()
  plot_Val$PC3 <- plot_Val$PC3 %>% as.numeric()

  #set object to store
  pca_output$PC1 <- list()
  #pca_output$PC1$sample_fail_idx <- list()
  #pca_output$PC1$qc_fail_idx <- list()
  
  #set object to store
  pca_output$PC2 <- list()
  #pca_output$PC2$sample_fail_idx <- list()
  #pca_output$PC2$qc_fail_idx <- list()
  
  #set object to store
  pca_output$PC3 <- list()
  #pca_output$PC3$sample_fail_idx <- list()
  #pca_output$PC3$qc_fail_idx <- list()
  
  #loop for each component
  # for(idx_PC in c("PC1", "PC2", "PC3")){
  # 
  # all_data <- plot_Val# %>% select(all_of(idx_PC))
  # FUNC_all_data_median <- median(all_data[[idx_PC]])
  # FUNC_all_data_sd <- sd(all_data[[idx_PC]])
  # FUNC_all_data_q1 <- quantile(all_data[[idx_PC]], 0.25) %>% as.numeric()
  # FUNC_all_data_q3 <- quantile(all_data[[idx_PC]], 0.75) %>% as.numeric()
  # FUNC_all_data_iqr <- IQR(all_data[[idx_PC]])
  # 
  # 
  # PC_threshold_low <- FUNC_all_data_q1 - (FUNC_all_data_iqr*FUNC_option_iqr_filter_samples)
  # PC_threshold_high <- FUNC_all_data_q3 + (FUNC_all_data_iqr*FUNC_option_iqr_filter_samples)
  # 
  # pca_output[[idx_PC]]$sample_fail_idx <- which(plot_Val[[idx_PC]] < PC_threshold_low | 
  #                                         plot_Val[[idx_PC]] > PC_threshold_high)
  #                                 #which(plot_Val$sample_type == "sample"))
  # 
  # 
  # 
  # #for qc data
  # 
  # FUNC_qc_data <- plot_Val %>% filter(sample_type == "qc") #%>% select(all_of(idx_PC))
  # FUNC_qc_data_median <- median(FUNC_qc_data[[idx_PC]])
  # FUNC_qc_data_sd <- sd(FUNC_qc_data[[idx_PC]])
  # FUNC_qc_data_q1 <- quantile(FUNC_qc_data[[idx_PC]], 0.25) %>% as.numeric()
  # FUNC_qc_data_q3 <- quantile(FUNC_qc_data[[idx_PC]], 0.75) %>% as.numeric()
  # FUNC_qc_data_iqr <- IQR(FUNC_qc_data[[idx_PC]])
  # 
  # PC_threshold_low_qc <- FUNC_qc_data_q1 - (FUNC_qc_data_iqr*FUNC_option_iqr_filter_qc)
  # PC_threshold_high_qc <- FUNC_qc_data_q3 + (FUNC_qc_data_iqr*FUNC_option_iqr_filter_qc)
  # 
  # pca_output[[idx_PC]]$qc_fail_idx <- intersect(which(plot_Val[[idx_PC]] < PC_threshold_low_qc | 
  #                                     plot_Val[[idx_PC]] > PC_threshold_high_qc),
  #                            which(plot_Val$sample_type == "qc"))
  
                     
  #pca_output[[idx_PC]]$fail_idx_unique <- c(pca_output[[idx_PC]]$sample_fail_idx, pca_output[[idx_PC]]$qc_fail_idx) %>% unique()

  #pca_output[[idx_PC]]$fail_sample <- FUNC_data$sample_name[pca_output[[idx_PC]]$fail_idx_unique]
  
  
  # set plot attributes (controlled by FUNC_colour_by and FUNC_plot_label)
  pca_colour <- list()
  pca_colour <- FUNC_data %>% select(all_of(FUNC_colour_by)) #%>% as.matrix()
  colnames(pca_colour) <- "pca_colour" 
  #pca_colour[[idx_PC]]$pca_colour <- factor(pca_colour[[idx_PC]]$pca_colour, levels = c(levels(pca_colour[[idx_PC]]$pca_colour), "sample fail", "qc fail"), ordered = TRUE)
  #pca_colour[[idx_PC]]$pca_colour[pca_output[[idx_PC]]$sample_fail_idx %>% unique()] <- "sample fail"
  #pca_colour[[idx_PC]]$pca_colour[pca_output[[idx_PC]]$qc_fail_idx %>% unique()] <- "qc fail"
  pca_plot_colour <- pca_colour$pca_colour
  pca_plot_colour[is.na(pca_plot_colour)] <- "none"
  
  #set colours
  plot_colours <- c(FUNC_project_colours)


 
  ##################  ##################  ##################  ##################
  #produce run order plot
  ##################  ##################  ##################  ##################
  
  for(idx_PC in plot_Val %>% select(-contains("sample")) %>% names()){
  
  bp <- ggplot(data=plot_Val,
               aes(x=sample_idx,
                   y=get(idx_PC))
  )
  
  bp <- bp + geom_point(aes(#text = sample_name,
                            fill = pca_plot_colour#,
                            #color = sample_type_factor
  ),
  shape = 21
  )
  
  bp <- bp + scale_fill_manual(values = c(plot_colours))
  #bp <- bp + scale_color_manual("black")
  
  bp <- bp + labs(x = paste("Sample order"),
                  y = paste0(idx_PC))
  bp <- bp + ggtitle(paste0(FUNC_title, " - ", idx_PC))
  bp <- bp + theme_cowplot() 
  bp <- bp + theme(
    plot.title = element_text(hjust = 0.5, size=14),
    axis.text.y = element_text(size = 12, margin = margin(t = 0, r = 0, b = 0, l = 2)),
    #axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
    axis.text.x = element_blank(),
    axis.title = element_text(size = 14),
    #legend.text=element_text(size=12),
    legend.title=element_blank(),
    axis.ticks.x = element_blank(),
    
  )
  
  #create vertical lines to separate classes on plot
  FUNC_data_batches <- plot_Val$sample_batch %>% unique()
  batch_idx <- NULL
   for(idx_batch in FUNC_data_batches[2:length(FUNC_data_batches)]){
     batch_idx <- c(batch_idx, min(which(plot_Val$sample_batch == idx_batch)))
   }

  bp <- bp + geom_vline(xintercept=c(batch_idx),color="grey")
  
  pca_output$plots[[idx_PC]]$plotly <- bp %>% ggplotly() %>% layout(legend = list(orientation = "h",   # show entries horizontally
                                                                                      xanchor = "center",  # use center of legend as anchor
                                                                                      x = 0.5,
                                                                                      y = -0.2,
                                                                                      title=list(text="Group comparison:"))
  )
     
  
 
  
  
  
  
  
  
   }
 
  pca_output$plots
}
