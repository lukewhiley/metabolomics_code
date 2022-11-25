#ANPC PCA quality control visualisation

# FUNC_compare_means_table = data containing individual lipid data - MUST CONTAIN sampleID column 
# FUNC_colour_by = how to colour the plot (e.g. sample class, is_ltr or cohort)
# FUNC_plot label = what to label the scores plot with (e.g. sampleID)
# FUNC_scaling = UV or Pareto

lgw_lipid_3D_plot <- function(FUNC_compare_means_table,
                              FUNC_plot_comparisons,
                              FUNC_plot_colour_low,
                              FUNC_plot_colour_high,
                              FUNC_sidechain_sep
                    ){

  lipid_plot_output <- list()
  
  #Add lipid class column to results table
  FUNC_compare_means_table <- FUNC_compare_means_table %>%
    add_column("FUNC_class" =
                 sub("\\(.*", "", FUNC_compare_means_table$feature),
               .after = "feature") %>%
    add_column("FUNC_feature_idx" = 
                 c(1:nrow(FUNC_compare_means_table)),
               .before = "feature")
  
  FUNC_sidechains <- list()
  #extract chain length data from inside brackets
  FUNC_sidechains$all <- regmatches(FUNC_compare_means_table$feature, gregexpr("\\(.*?\\)", FUNC_compare_means_table$feature)) %>% unlist()
  FUNC_sidechains$all <- sub("\\(", "", FUNC_sidechains$all)
  FUNC_sidechains$all <- sub("\\)", "", FUNC_sidechains$all)
  
  #extract data either side of "/" for sidechain 1 and sidechain 2
  #FUNC_sidechains$sidechain_1 <- sub('/.*', '', FUNC_sidechains$all)
  FUNC_sidechains$sidechain_1 <- sub(paste0(FUNC_sidechain_sep, ".*"), '', FUNC_sidechains$all)
  #FUNC_sidechains$sidechain_2 <- sub('.*/', '', FUNC_sidechains$all)
  FUNC_sidechains$sidechain_2 <- sub(paste0(".*", FUNC_sidechain_sep), '', FUNC_sidechains$all)
  
  #delete TAG from sidechain 1 because they do not contain specific 1x sidechain data
  FUNC_sidechains$sidechain_1[which(FUNC_compare_means_table$FUNC_class == "TAG")] <- NA
  
  #delete single chain lipids from sidechain 1 because they only have 1x chain and it will be stored in sidechain 2 in this table
  FUNC_sidechains$sidechain_1[!grepl(FUNC_sidechain_sep, FUNC_compare_means_table$feature)] <- NA

  
  #drop the extra P- and O- and d, allows for focus on just chain length
  FUNC_sidechains$sidechain_1 <- sub('P-', '', FUNC_sidechains$sidechain_1)
  FUNC_sidechains$sidechain_1 <- sub('O-', '', FUNC_sidechains$sidechain_1)
  FUNC_sidechains$sidechain_1 <- sub('d', '', FUNC_sidechains$sidechain_1)
  FUNC_sidechains$sidechain_2 <- sub('FA', '', FUNC_sidechains$sidechain_2)
  
  #melt into longer list and stack side chains into single column
  FUNC_sidechains$plot <- FUNC_compare_means_table %>% 
    add_column(sidechain = FUNC_sidechains$sidechain_1) %>%
    bind_rows(FUNC_compare_means_table %>% 
                add_column(sidechain = FUNC_sidechains$sidechain_2)) %>%
    filter(!is.na(sidechain))
  
  plot_Val <- FUNC_sidechains$plot %>% 
    select(FUNC_feature_idx, sidechain, FUNC_class, feature, p, all_of(FUNC_plot_comparisons))  
  
  ###########
  #pivot to a long table to control plotting scales
  plot_Val_long <- plot_Val %>% 
    select(-p) %>%
    pivot_longer(
      cols = all_of(FUNC_plot_comparisons),
      names_to = "comparison",
      values_to = "p")
  plot_Val_long <- plot_Val_long %>%
    #add column for -Log10 p
    add_column("-Log10 (p)" = 
                 plot_Val_long$p %>% log() %>% abs()) %>%
    #add column for lipid class as unique ordered factors
    add_column("lipid_class_factor" =
                 plot_Val_long$FUNC_class %>% 
                 factor(levels = unique(plot_Val_long$FUNC_class %>% sort()), 
                        ordered = TRUE)) %>%
    #add column for lipid sidechains as unique ordered factors
    add_column("lipid_sidechain_factor" =
                 plot_Val_long$sidechain %>% 
                 factor(levels = unique(plot_Val_long$sidechain %>% sort()), 
                        ordered = TRUE))
    #add column for plot point size as rank ordered -log10 p 
  plot_Val_long <-  plot_Val_long %>%
    add_column("point_size" = 
                 rank(plot_Val_long$`-Log10 (p)`))
  
  #reset plot point size scale so that all features p<0.05 are set to size 1
  plot_Val_long$point_size[which(plot_Val_long$p > 0.05)] <- 1
  #create plot point size integer breaks
  plot_Val_long$point_size[which(plot_Val_long$p < 0.05)] <- as.integer(
    cut(plot_Val_long$point_size[which(plot_Val_long$p < 0.05)], 
        breaks = 100))
  # define plot point size integer break
  plot_Val_long$point_size[intersect(which(plot_Val_long$point_size > 1),
                                     which(plot_Val_long$point_size < 26))] <- 2.5
  plot_Val_long$point_size[intersect(which(plot_Val_long$point_size > 2.5),
                                     which(plot_Val_long$point_size < 51))] <- 5
  plot_Val_long$point_size[intersect(which(plot_Val_long$point_size > 5),
                                     which(plot_Val_long$point_size < 81))] <- 10
  plot_Val_long$point_size[intersect(which(plot_Val_long$point_size > 10),
                                     which(plot_Val_long$point_size < 81))] <- 20
  plot_Val_long$point_size[which(plot_Val_long$point_size > 20)] <- 40
 
  #find max value (resulting from smallest p) for plot axis
  max_log10_p <- plot_Val_long$`-Log10 (p)` %>% max()
  
  #run loop to create boxplot
for(idx_str_data in FUNC_plot_comparisons){
  plot_Val_2 <- plot_Val_long %>%
    filter(comparison == idx_str_data)

  #produce ggboxplot
  
  bp <- ggplot(data=plot_Val_2,
               aes(x=lipid_class_factor,
                   y=lipid_sidechain_factor))
    #add points
  bp <- bp + geom_quasirandom(
    colour = "black",
    shape = 21,
    aes(fill = `-Log10 (p)`,
    size = point_size))
  
  bp <- bp + labs(x = paste("Lipid Class"),
                  y = paste("Sidechain"))
  bp <- bp + ggtitle(idx_str_data)
  #bp <- bp + ggtitle(paste0(str_split(idx_str_data, " - ")[[1]][1], " vs ", str_split(idx_str_data, " - ")[[1]][2]))
  bp <- bp + theme_cowplot() 
  bp <- bp + theme(plot.title = element_text(hjust = 0.5)) 
  bp <- bp + theme(plot.title = element_text(size=14)) 
  bp <- bp + theme(axis.text.y = element_text(size = 12, margin = margin(t = 0, r = 0, b = 0, l = 2)))
  bp <- bp + theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1))
  bp <- bp + theme(axis.title = element_text(size = 14)) 
  bp <- bp + theme(legend.title=element_text(size=12), 
                   legend.text=element_text(size=12))
  bp <- bp + scale_size(limits = c(0,100))
  bp <- bp + scale_fill_gradient(low = FUNC_plot_colour_low, high = FUNC_plot_colour_high,
                                 limits = c(0, max_log10_p))
  bp <- bp + guides(size="none")
  
  #create vertical lines to seprate classes on plot
  x_lipid_sequence <- seq(1:(length(plot_Val_2$lipid_class_factor %>% unique())-1))+0.5  
  y_lipid_sequence <- seq(1:(length(plot_Val_2$lipid_sidechain_factor %>% unique())-1))+0.5 
  
  bp <- bp + geom_vline(xintercept=c(x_lipid_sequence),color="grey")
  bp <- bp + geom_hline(yintercept=c(y_lipid_sequence),color="grey")

 
  lipid_plot_output[[idx_str_data]] <- bp
  
  #lipid_plot_output[[idx_str_data]] <- list()
  #lipid_plot_output[[idx_str_data]]$ggplot <- bp 
  # #create plot_ly
  #lipid_plot_output[[idx_str_data]]$plotly <- bp %>% ggplotly()
  
  }
  lipid_plot_output
  
}


#Archive backup code

#bp <- bp + scale_color_gradient(low = FUNC_plot_colour_low, high = FUNC_plot_colour_high)
#bp <- bp + scale_fill_viridis_c(option = "magma", limits = c(0, max_log10_p))
#bp <- bp + scale_color_viridis_c(option = "magma", limits = c(0, max_log10_p))
#bp$labels$fill <- paste0(FUNC_HEADER_temp_colour) %>% str_to_title()
  

