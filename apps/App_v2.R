library(shiny)
library(janitor)
library(tidyverse)
library(png)
library(vroom)

# theme for ggplot 2

#plotting theme for ggplot2
.theme<- theme(
  axis.line = element_line(colour = 'gray', size = .75),
  panel.background = element_blank(),
  plot.background = element_blank()
)


# ui object
ui_title <- titlePanel(title=div(img(src="ANPC_logo.png", width = "240px", height = "140px"), 
                                   "LGW - Metabolomics Workbook", style = "color:#3474A7"))

#create a user interface for uplopading the metabolite data values
ui_data_upload <- sidebarLayout(
  
  #side panel
  sidebarPanel(
    p("select a .csv or .tsv containing metabolite concentration values. Column headers should contain the metabolite name matching the metabolite list below. sample/filenames should be in a column called: sample name"),
    fileInput("metabolite_data_file",
              label = "choose file",
              multiple = TRUE,
              accept = c(".csv", ".tsv")),
    p("does the file have a header?"),
    checkboxInput("header", "Header", TRUE),
    radioButtons("sep", "Separator",
                 choices = c(Semicolon = ";",
                             Comma = ",",
                             Tab = "\t"),
                 selected = ","),
    numericInput("r1", "Rows to preview", 5, min = 1),
    numericInput("c1", "Columns to preview", 5, min = 1)
  ),
  
mainPanel(
    h3("Raw data"),
    tableOutput("metabolite_data_file")
  )
)

#create a user interface for uploading a list of target metabolites
ui_metabolite_target_upload <- sidebarLayout(
  sidebarPanel(
    p("select a .csv or .tsv contatining a list of metabolite targets in a column named: metabolite"),
    fileInput("metabolite_target_file", NULL, accept = c(".csv", ".tsv")),
    numericInput("r2", "Rows to preview", 5, min = 1),
  ),
  mainPanel(
    h3("Metabolite targets"),
    tableOutput("metabolite_target_file")
  )
)

#create a user interface for uploading corresponding metadata
ui_metadata_upload <- sidebarLayout(
  sidebarPanel(
    p("select a .csv or .tsv contatining any project metadata. Column headers should contain the metadata headers. Sample/filenames should be in a column called: sample name. This should match the raw data file uploaded above"),
    fileInput("metabolite_metadata_file", NULL, accept = c(".csv", ".tsv")),
    numericInput("r3", "Rows to preview", 5, min = 1),
    numericInput("c3", "Columns to preview", 5, min = 1)
  ),
  mainPanel(
    h3("Metadata"),
    tableOutput("metabolite_metadata_file")
  )
)


#create a user interface for making boxplots
ui_boxplot <- pageWithSidebar(
  # title
  headerPanel("Select Options"),
  
  #input
  sidebarPanel
  (
    # Input: Select a file ----
    
    fileInput("file1", "Choose CSV File",
              multiple = TRUE,
              accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")),
    # Input: Checkbox if file has header ----
    checkboxInput("header", "Header", TRUE),
    
    # Input: Select separator ----
    radioButtons("sep", "Separator",
                 choices = c(Semicolon = ";",
                             Comma = ",",
                             Tab = "\t"),
                 selected = ","),
    # Horizontal line ----
    tags$hr(),
    
    
    # Input: Select what to display
    selectInput("dataset","Data:",
                choices =list(iris = "iris",
                              uploaded_file = "inFile"), selected=NULL),
    selectInput("variable","Variable:", choices = NULL),
    selectInput("group","Group:", choices = NULL),
    selectInput("plot.type","Plot Type:",
                list(boxplot = "boxplot", histogram = "histogram", density = "density", bar = "bar")
    ),
    checkboxInput("show.points", "show points", TRUE)
  ),
  
  # output
  mainPanel(
    h3(textOutput("caption")),
    #h3(htmlOutput("caption")),
    uiOutput("plot") # depends on input
  )
)
  



ui <- fluidPage(
  ui_title,
  ui_data_upload,
  ui_metabolite_target_upload,
  ui_metadata_upload,
  ui_boxplot
)


# server()
server<-(function(input, output, session){
    #update group and
    #variables based on the data
    observe({
      #browser()
      if(!exists(input$dataset)) return() #make sure upload exists
      var.opts<-colnames(get(input$dataset))
      updateSelectInput(session, "variable", choices = var.opts)
      updateSelectInput(session, "group", choices = var.opts)
    })
    
    output$caption<-renderText({
      switch(input$plot.type,
             "boxplot" 	= 	"Boxplot",
             "histogram" =	"Histogram",
             "density" 	=	"Density plot",
             "bar" 		=	"Bar graph")
    })
    
    output$plot <- renderUI({
      plotOutput("p")
    })
    
    #get data object
    get_data<-reactive({
      
      if(!exists(input$dataset)) return() # if no upload
      
      check<-function(x){is.null(x) || x==""}
      if(check(input$dataset)) return()
      
      obj<-list(data=get(input$dataset),
                variable=input$variable,
                group=input$group
      )
      
      #require all to be set to proceed
      if(any(sapply(obj,check))) return()
      #make sure choices had a chance to update
      check<-function(obj){
        !all(c(obj$variable,obj$group) %in% colnames(obj$data))
      }
      
      if(check(obj)) return()
      
      
      obj
      
    })
    
    #plotting function using ggplot2
    output$p <- renderPlot({
      
      plot.obj<-get_data()
      
      #conditions for plotting
      if(is.null(plot.obj)) return()
      
      #make sure variable and group have loaded
      if(plot.obj$variable == "" | plot.obj$group =="") return()
      
      #plot types
      plot.type<-switch(input$plot.type,
                        "boxplot" 	= geom_boxplot(),
                        "histogram" =	geom_histogram(alpha=0.5,position="identity"),
                        "density" 	=	geom_density(alpha=.75),
                        "bar" 		=	geom_bar(position="dodge")
      )
      
      
      if(input$plot.type=="boxplot")	{		#control for 1D or 2D graphs
        p<-ggplot(plot.obj$data,
                  aes_string(
                    x 		= plot.obj$group,
                    y 		= plot.obj$variable,
                    fill 	= plot.obj$group # let type determine plotting
                  )
        ) + plot.type
        
        if(input$show.points==TRUE)
        {
          p<-p+ geom_point(color='black',alpha=0.5, position = 'jitter')
        }
        
      } else {
        
        p<-ggplot(plot.obj$data,
                  aes_string(
                    x 		= plot.obj$variable,
                    fill 	= plot.obj$group,
                    group 	= plot.obj$group
                    #color 	= as.factor(plot.obj$group)
                  )
        ) + plot.type
      }
      
      p<-p+labs(
        fill 	= input$group,
        x 		= "",
        y 		= input$variable
      )  +
        .theme
      print(p)
    })
    
    # set uploaded file
    upload_data<-reactive({
      
      inFile <- input$metabolite_data_file
      
      if (is.null(inFile))
        return(NULL)
      
      #could also store in a reactiveValues
      read.csv(inFile$datapath,
               header = input$header,
               sep = input$sep)
    })
    
    observeEvent(input$fmetabolite_data_file,{
      inFile<<-upload_data()
    })
    
    
  })
  
  
#   #server on importing raw data
#   metabolite_data_file <- reactive({
#     req(input$metabolite_data_file)
#     
#     ext <- tools::file_ext(input$metabolite_data_file$name)
#     switch(ext,
#            csv = vroom::vroom(input$metabolite_data_file$datapath, delim = ","),
#            tsv = vroom::vroom(input$metabolite_data_file$datapath, delim = "\t"),
#            validate("Invalid file; Please upload a .csv or .tsv file")
#     )
#   })
#   
#   output$metabolite_data_file <- renderTable({
#     head(metabolite_data_file(), c(input$r1, input$c1))
#     
#   })
#   
#   #server on importing metabolite target list
#   metabolite_target_file <- reactive({
#     req(input$metabolite_target_file)
#     
#     ext <- tools::file_ext(input$metabolite_target_file$name)
#     switch(ext,
#            csv = vroom::vroom(input$metabolite_target_file$datapath, delim = ","),
#            tsv = vroom::vroom(input$metabolite_target_file$datapath, delim = "\t"),
#            validate("Invalid file; Please upload a .csv or .tsv file")
#     )
#   })
#   
#   output$metabolite_target_file <- renderTable({
#     head(metabolite_target_file(), input$r2)
#     
#   })
# 
# #server on importing metadata
# metabolite_metadata_file <- reactive({
#   req(input$metabolite_metadata_file)
#   
#   ext <- tools::file_ext(input$metabolite_metadata_file$name)
#   switch(ext,
#          csv = vroom::vroom(input$metabolite_metadata_file$datapath, delim = ","),
#          tsv = vroom::vroom(input$metabolite_metadata_file$datapath, delim = "\t"),
#          validate("Invalid file; Please upload a .csv or .tsv file")
#   )
# })
# 
# output$metabolite_metadata_file <- renderTable({
#   head(metabolite_metadata_file(), c(input$r3, input$c3))
#   
# })
# 
# 
# # server for boxplot ################################################
# 
# # data <- reactive({
# #   inFile <- input$metabolite_data_file
# #   if(!is.null(inFile)){
# #     read.csv(inFile$datapath, header = input$header, stringsAsFactors = FALSE)    
# #   }
# # })
# # 
# # output$plot1 <- renderPlot({
# #   req(data())
# #   print("help") 
# # })
# 
# plot_data <- reactive({
#   if (is.null(input$metabolite_data_file)) return(NULL)
#   read.csv(input$metabolite_data_file$datapath, header = input$header)
# })
# 
# output$contents <- renderTable({
#   req(plot_data()) # if user_data() is null than the reactive will not run
#   head(plot_data(), c(input$r1, input$c1))
# })
# 
# output$plot1 <- renderPlot({
#   #req(user_data())
#   #plot(user_data$speed, user_data$dist)
# })
# 
# 
# 
# 
# })#final bracket for server

# shinyApp()
shinyApp(ui = ui, server = server)

