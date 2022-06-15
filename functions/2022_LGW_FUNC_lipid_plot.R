#ANPC PCA quality control visualisation

# FUNC_data = data containing individual lipid data - MUST CONTAIN sampleID column 
# FUNC_colour_by = how to colour the plot (e.g. sample class, is_ltr or cohort)
# FUNC_plot label = what to label the scores plot with (e.g. sampleID)
# FUNC_scaling = UV or Pareto

lgw_lipid_plot <- function(FUNC_data,
                           FUNC_plot_comparisons,
                    #FUNC_colour_by, 
                    #FUNC_plot_label, 
                    #FUNC_title,
                    FUNC_project_colour_reference,
                    FUNC_project_colours,
                    FUNC_plot_class_or_sidechain){
  #require(metabom8)
  require(RColorBrewer)
  require(tidyverse)
  require(plotly)
  
  lipid_plot_output <- list()
  
  #browser()
  
  #title_text <- FUNC_title
  
  # create plot values
  if(FUNC_plot_class_or_sidechain == "class"){
  plot_Val <- FUNC_data %>% 
    select(feature_idx, subclass, feature, p, all_of(FUNC_plot_comparisons)) %>%
    rename(lipid_class = subclass)
}
if(FUNC_plot_class_or_sidechain == "sidechain"){

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
  select(feature_idx, sidechain, subclass, feature, p, all_of(FUNC_plot_comparisons)) %>%
  rename(lipid_class = subclass)
}
  
  
  
  temp_list <- list()
  for(idx_plot_comparisons in 1:length(FUNC_plot_comparisons)){
    temp_column_idx <- which(colnames(plot_Val) %in% FUNC_plot_comparisons[idx_plot_comparisons])
    temp_list[[idx_plot_comparisons]] <- plot_Val %>%
      add_column(plot_Val[,temp_column_idx] %>% 
                   log() %>% 
                   abs() %>% 
                   setNames(paste0("-Log10 (p)"))) %>%
      add_column(`Dunn's test comparison` = factor(FUNC_plot_comparisons[idx_plot_comparisons]))
 }
  
  plot_Val_2 <- bind_rows(temp_list)
  
  #browser()
  
  #create factor for subclass
  if(FUNC_plot_class_or_sidechain=="class"){
    plot_Val_2$lipid_class_factor <- plot_Val_2$lipid_class %>% factor(levels = unique(plot_Val_2$lipid_class %>% sort()), ordered = TRUE)
  }
  
  if(FUNC_plot_class_or_sidechain=="sidechain"){
  plot_Val_2$lipid_sidechain_factor <- plot_Val_2$sidechain %>% factor(levels = unique(plot_Val_2$sidechain %>% sort()), ordered = TRUE)
  }

  #create values for point size
  plot_Val_2 <- plot_Val_2 %>%
    add_column(point_size = rank(plot_Val_2$`-Log10 (p)`))
  
  #reset scale so that all features p<0.05 are set to size 1
  plot_Val_2$point_size <- plot_Val_2$point_size - which(plot_Val_2$`-Log10 (p)` < (log(0.05) %>% abs())) %>% length()
  
  #reset scale so that all features p<0.05 are set to size 1
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
    arrange(feature_idx)
  
  #produce plot
  if(FUNC_plot_class_or_sidechain == "class"){
  bp <- ggplot(data=plot_Val_2,
               aes(x=lipid_class_factor,
                   y=`-Log10 (p)`),
               group =`Dunn's test comparison`)
  }
  
  if(FUNC_plot_class_or_sidechain == "sidechain"){
    bp <- ggplot(data=plot_Val_2,
                 aes(x=lipid_sidechain_factor,
                     y=`-Log10 (p)`),
                 group =`Dunn's test comparison`)
  }
  
  #browser()
  #create color reference
  unique_comparisons <- unique(plot_Val_2$`Dunn's test comparison`)
  unique_comparison_test <- sub('.*_', '', unique_comparisons) %>% str_to_lower()
  colour_idx <- which((FUNC_project_colour_reference%>% str_to_lower()) %in% unique_comparison_test)

    
  # add scatter points to plot
  bp <- bp + geom_point(aes(#text = feature,
                            fill = `Dunn's test comparison`,
                            group =`Dunn's test comparison`,
                            size = point_size
                            ),
                          position=position_dodge(0.7),
                          #width = 0.01,
                          shape = 16
                         )
  bp <- bp + scale_fill_manual(values = FUNC_project_colours[colour_idx])
  
  
  #edit titles
if(FUNC_plot_class_or_sidechain=="class"){
  plot_x_title <- "lipid class"
}
  if(FUNC_plot_class_or_sidechain=="sidechain"){
    plot_x_title <- "lipid sidechain"
  }
  bp <- bp + labs(x = paste(plot_x_title),
                  y = paste("Kruskal Wallis post-hoc Dunn's Test [-Log10 (p)]"))
  bp <- bp + ggtitle("Kruskal Wallis post-hoc Dunn's Test Results")
  
  #edit theme and plot appearance
  bp <- bp + theme_cowplot() 
  bp <- bp + theme(plot.title = element_text(hjust = 0.5)) 
  bp <- bp + theme(plot.title = element_text(size=10)) 
  bp <- bp + theme(axis.text.y = element_text(size = 10, margin = margin(t = 0, r = 0, b = 0, l = 2)))
  bp <- bp + theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust = 1))
  bp <- bp + theme(axis.title = element_text(size = 10)) 
  bp <- bp + theme(legend.title=element_text(size=10), 
                   legend.text=element_text(size=10))
  
  #add line guides to plot
  #create vertical lines to seprate classes on plot
  if(FUNC_plot_class_or_sidechain == "class"){lipid_sequence <- seq(1:(length(plot_Val_2$lipid_class %>% unique())-1))+0.5}
  if(FUNC_plot_class_or_sidechain == "sidechain"){lipid_sequence <- seq(1:(length(plot_Val_2$lipid_sidechain_factor %>% unique())-1))+0.5}
  #add the lines to the plot
  bp <- bp + geom_vline(xintercept=c(lipid_sequence),color="grey")
  
  #add horizontal line to represent p=0.05 threshold
  h_line <- c(0.05 %>% log() %>% abs())
  bp <- bp + geom_hline(yintercept=h_line, color="red")
  
  #label horizontal line
  # bp <- bp + geom_text(aes(x = 21, 
  #                          y= h_line, 
  #                          label = "p=0.05", 
  #                          vjust = - 1, 
  #                          hjust = 2))
  
  #remove legends
  bp <- bp + guides(size= "none")
  
  
 #add scaleing factor that controls points on the plot
  bp <- bp + scale_size(limits = c(0,100))
  
 
  
  lipid_plot_output <- bp %>%
    ggplotly() %>%
    layout(legend = list(orientation = "h",   # show entries horizontally
                         xanchor = "center",  # use center of legend as anchor
                         x = 0.5,
                         y = -0.2,
                         title=list(text="Group comparison:"))
    )
  lipid_plot_output
  
}


  

