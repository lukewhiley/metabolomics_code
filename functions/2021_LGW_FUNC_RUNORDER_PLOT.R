#ANPC plot signal over the course of a run plot



lgw_runorder_plot <- function(FUNC_data, 
                        FUNC_metabolite_list, 
                        FUNC_HEADER_colour_by, 
                        FUNC_HEADER_plot_label,
                        FUNC_HEADER_batch,
                        FUNC_title,
                        FUNC_project_colours
                    #FUNC_option_invert_y,
                    #FUNC_option_invert_x,
                    #FUNC_option_plot_qc
                    ){
  
  require(tidyverse)
  require(ggplot2)
  
  run_plot_output <- list()
  
  #browser()
  
  title_text <- FUNC_title
  
  for(idx_feature in FUNC_metabolite_list){
   # browser()
    # create plot values
  plot_Val <- FUNC_data %>% 
    select(
      all_of(FUNC_HEADER_colour_by), 
      all_of(FUNC_HEADER_plot_label),
      all_of(FUNC_HEADER_batch),
      all_of(idx_feature)
      ) %>%
    rename("concentration" = all_of(idx_feature))
  #log data
  plot_Val$concentration <- (plot_Val$concentration + 1) %>% log()
  
  #create sample idx factor for plotting
  plot_Val$sample_idx <- seq(1:nrow(plot_Val)) %>% factor()


    #produce plot
  bp <- ggplot(data=plot_Val,
               aes(x=sample_idx,
                   y=concentration)
  )
  
  bp <- bp + geom_point(aes(text = sample_name,
                            fill = sample_type_factor#,
                            #color = sample_type_factor
  ),
  shape = 21
  )
  
  bp <- bp + scale_fill_manual(values = c(FUNC_project_colours))
  #bp <- bp + scale_color_manual("black")
  
  bp <- bp + labs(x = paste("Sample order"),
                  y = paste("Log[Summed Class Concentration]"))
  bp <- bp + ggtitle(paste0(FUNC_title, " - ", idx_feature))
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

  run_plot_output[[idx_feature]]$ggplot <- bp
  run_plot_output[[idx_feature]]$plotly <- bp %>% ggplotly() %>% layout(legend = list(orientation = 'h'))
  }

  run_plot_output
  }
