lgw_colour_pie <- function(){
  
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

pie(rep(1, 26), col = colour_pie)

# print(paste0("what ", FUNC_number_of_groups, " colours do you want to use? Seperate each one with a comma"))
# user_input <- readline()
# 
# colour_choice <- colour_pie[user_input]
# 
# pie(rep(1, FUNC_number_of_groups), col = colour_choice)

}