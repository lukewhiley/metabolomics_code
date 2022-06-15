#ANPC PCA quality control visualisation

# FUNC_data = data containing individual lipid data - MUST CONTAIN sampleID column 
# FUNC_colour_by = how to colour the plot (e.g. sample class, is_ltr or cohort)
# FUNC_plot label = what to label the scores plot with (e.g. sampleID)
# FUNC_scaling = UV or Pareto

lgw_lipid_3D_plot <- function(FUNC_data,
                           FUNC_plot_comparisons,
                           FUNC_plot_colour_low,
                           FUNC_plot_colour_high
                    ){
  require(metabom8)
  require(RColorBrewer)
  require(tidyverse)
  require(plotly)
  
  lipid_plot_output <- list()
  
  #browser()
  
  FUNC_sidechains <- list()
  #extract chain length data from inside brackets
  FUNC_sidechains$all <- gsub("[\\(\\)]", "", regmatches(FUNC_data$feature, gregexpr("\\(.*?\\)", FUNC_data$feature)))
  
  #extract data either side of "/" for sidechain 1 and sidechain 2
  FUNC_sidechains$sidechain_1 <- sub('/.*', '', FUNC_sidechains$all)
  FUNC_sidechains$sidechain_2 <- sub('.*/', '', FUNC_sidechains$all)
  
  #delete TAG from sidechain 1 because they do not contain specific sidechain data
  FUNC_sidechains$sidechain_1[which(master_list$lipid_plot$data$subclass == "TAG")] <- NA
  FUNC_sidechains$sidechain_1[!grepl("/", master_list$lipid_plot$data$feature)] <- NA
  
  
  #drop the extra P- and O- and d, allows for focus on just chain length
  FUNC_sidechains$sidechain_1 <- sub('P-', '', FUNC_sidechains$sidechain_1)
  FUNC_sidechains$sidechain_1 <- sub('O-', '', FUNC_sidechains$sidechain_1)
  FUNC_sidechains$sidechain_1 <- sub('d', '', FUNC_sidechains$sidechain_1)
  FUNC_sidechains$sidechain_2 <- sub('FA', '', FUNC_sidechains$sidechain_2)
  
  #melt into longer list and stack side chains into single column
  FUNC_sidechains$plot <- FUNC_data %>% 
    add_column(sidechain = FUNC_sidechains$sidechain_1) %>%
    bind_rows(FUNC_data %>% 
                add_column(sidechain = FUNC_sidechains$sidechain_2)) %>%
    filter(!is.na(sidechain))
  
  plot_Val <- FUNC_sidechains$plot %>% 
    select(feature_idx, sidechain, subclass, feature, p, all_of(FUNC_plot_comparisons)) 
  
  #find max -log10 (p)
 max_log10_p <- FUNC_data %>% 
    select(all_of(FUNC_plot_comparisons)) %>%
    log() %>%
    abs() %>%
    max()
 

  for(idx_str_data in FUNC_plot_comparisons){
  temp_column_idx <- which(colnames(plot_Val) == idx_str_data)
    plot_Val_2 <- plot_Val %>%
      add_column(plot_Val[,temp_column_idx] %>%
                   log() %>%
                   abs() %>%
                   setNames(paste0("-Log10 (p)")))
 

  #create factor for subclass
  plot_Val_2$lipid_class_factor <- plot_Val_2$subclass %>% factor(levels = unique(plot_Val_2$subclass %>% sort()), ordered = TRUE)
  
  plot_Val_2$lipid_sidechain_factor <- plot_Val_2$sidechain %>% factor(levels = unique(plot_Val_2$sidechain %>% sort()), ordered = TRUE)

  
  #create values for point size
  plot_Val_2 <- plot_Val_2 %>%
    add_column(point_size = rank(plot_Val_2$`-Log10 (p)`))
  
  #reset scale so that all features p<0.05 are set to size 1
  plot_Val_2$point_size <- plot_Val_2$point_size - which(plot_Val_2$`-Log10 (p)` < (log(0.05) %>% abs())) %>% length()
  
  #ser scale so that all features p<0.05 are set to size 1
  plot_Val_2$point_size[which(plot_Val_2$`-Log10 (p)` < (log(0.05) %>% abs()))] <- 1
 
  #arrange by point_size
  plot_Val_2 <- plot_Val_2 %>% arrange(`-Log10 (p)`)
  
  plot_Val_2$point_size_rank <- 1
  
  #split off all features that are significant
  temp_plot_Val <- plot_Val_2 %>%
    filter(point_size >1)
  
  #set breaks
  temp_plot_Val$point_size_rank <- as.integer(cut(temp_plot_Val$point_size, breaks = 100))
  #define size breaks
  temp_plot_Val$point_size <- 2.5
  temp_plot_Val$point_size[which(temp_plot_Val$point_size_rank>25)] <- 5
  temp_plot_Val$point_size[which(temp_plot_Val$point_size_rank>50)] <- 10
  temp_plot_Val$point_size[which(temp_plot_Val$point_size_rank>80)] <- 20
  temp_plot_Val$point_size[which(temp_plot_Val$point_size_rank>95)] <- 40
  
  #re-build plot_Val_2 table
  plot_Val_2 <- plot_Val_2 %>%
    filter(point_size == 1) %>%
    bind_rows(temp_plot_Val) %>% 
    #arrange(feature_idx)
    arrange(`-Log10 (p)`)
    #arrange(desc(sidechain))
  
  #browser()
  
  #produce plot
  bp <- ggplot(data=plot_Val_2,
               aes(x=lipid_class_factor,
                   y=lipid_sidechain_factor)
               )
  
 

    #add points
  bp <- bp + geom_quasirandom(#cex = 2,
    #position = "beeswarm",
                        #dodge.width = 0.1,
                        #size = 0.5,
                        colour = "black",
                        shape = 21,
                        aes(fill = `-Log10 (p)`,
                            #color = `-Log10 (p)`,
                            size = point_size
                            ))
  
  bp <- bp + labs(x = paste("Lipid Class"),
                  y = paste("Sidechain"))
  bp <- bp + ggtitle(paste0(str_split(idx_str_data, "_")[[1]][1], " vs ", str_split(idx_str_data, "_")[[1]][2]))
  bp <- bp + theme_cowplot() 
  bp <- bp + theme(plot.title = element_text(hjust = 0.5)) 
  bp <- bp + theme(plot.title = element_text(size=14)) 
  bp <- bp + theme(axis.text.y = element_text(size = 12, margin = margin(t = 0, r = 0, b = 0, l = 2)))
  bp <- bp + theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1))
  bp <- bp + theme(axis.title = element_text(size = 14)) 
  bp <- bp + theme(legend.title=element_text(size=12), 
                   legend.text=element_text(size=12))
  bp <- bp + scale_size(limits = c(0,100))
  bp <- bp + scale_fill_gradient(low = FUNC_plot_colour_low, high = FUNC_plot_colour_high)
  bp <- bp + guides(size="none")
  
  #create vertical lines to seprate classes on plot
  x_lipid_sequence <- seq(1:(length(plot_Val_2$lipid_class_factor %>% unique())-1))+0.5  
  y_lipid_sequence <- seq(1:(length(plot_Val_2$lipid_sidechain_factor %>% unique())-1))+0.5 
  
  bp <- bp + geom_vline(xintercept=c(x_lipid_sequence),color="grey")
  bp <- bp + geom_hline(yintercept=c(y_lipid_sequence),color="grey")
  
  #bp <- bp + scale_color_gradient(low = FUNC_plot_colour_low, high = FUNC_plot_colour_high)

  
  
  #bp <- bp + scale_fill_viridis_c(option = "magma", limits = c(0, max_log10_p))
  #bp <- bp + scale_color_viridis_c(option = "magma", limits = c(0, max_log10_p))
  #bp$labels$fill <- paste0(FUNC_HEADER_temp_colour) %>% str_to_title()
  lipid_plot_output[[idx_str_data]] <- list()
  lipid_plot_output[[idx_str_data]]$ggplot <- bp 
  
  # #create plot_ly
  lipid_plot_output[[idx_str_data]]$plotly <- bp %>% ggplotly()
  
  }
  lipid_plot_output
  
}
  

