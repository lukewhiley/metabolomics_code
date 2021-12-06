lgw_colour_pie_select <- function(FUNC_OPTION_colour_choice){
  
  colour_pie <- c(
    "dodgerblue2", 
    "#E31A1C", 
    "green4",
    "#6A3D9A", 
    "#FF7F00", 
    "black", 
    "white",
    "gold1",
    "skyblue2", 
    "#FB9A99", 
    "palegreen2",
    "#CAB2D6", 
    "#FDBF6F", 
    "gray70", 
    "khaki2",
    "maroon", 
    "orchid1", 
    "deeppink1", 
    "blue1", 
    "steelblue4",
    "darkturquoise", 
    "green1", 
    "yellow4", 
    "yellow3",
    "darkorange4", 
    "brown"
  )

  pie(rep(1, length(FUNC_OPTION_colour_choice)), col = colour_pie[FUNC_OPTION_colour_choice])
  
  paste0(colour_pie[FUNC_OPTION_colour_choice])
  
}