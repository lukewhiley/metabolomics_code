lgw_boxplot <- function(x){
  ggplot(data=x, 
         aes(x=as.factor(class), 
             y=as.numeric(concentration))) + 		
    geom_boxplot(aes(color = class, fill = class), outlier.shape = NA, lwd = 1, alpha=0.05)  +
    geom_jitter(aes(color = class), width = 0.05) +
    labs(x = paste(""), y = paste("Concentration µM")) +
    ggtitle(idx_bp) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size=50)) +
    theme(axis.text.y=element_text(size = 40, margin = margin(t = 0, r = 0, b = 0, l = 2))) +
    theme(axis.title = element_text(size = 45)) +
    theme(axis.text.x = element_blank()) +
    theme(legend.position = "right",  
          legend.title = element_text(color = "black", size = 40),
          legend.text = element_text(color = "black", size = 40),
          legend.key.size = unit(2, "cm"))
}
