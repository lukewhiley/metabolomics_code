library(shiny)
library(janitor)
library(tidyverse)
library(png)
library(vroom)

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
ui_boxplot <- sidebarLayout(
  sidebarPanel(
    p("test boxplot"),
    selectInput("dataset","Data:",
                choices =list(iris = "iris", mtcars = "mtcars",
                              uploaded_file = "inFile"), selected=NULL)
  ),
    mainPanel(
      mainPanel(
        h3("boxplot"),
        plotOutput(outputId = "plot1"), # depends on input
        tableOutput(outputId = "content")
      )
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
server <- shinyServer(function(input, output, session) {
  
  #server on importing raw data
  metabolite_data_file <- reactive({
    req(input$metabolite_data_file)
    
    ext <- tools::file_ext(input$metabolite_data_file$name)
    switch(ext,
           csv = vroom::vroom(input$metabolite_data_file$datapath, delim = ","),
           tsv = vroom::vroom(input$metabolite_data_file$datapath, delim = "\t"),
           validate("Invalid file; Please upload a .csv or .tsv file")
    )
  })
  
  output$metabolite_data_file <- renderTable({
    head(metabolite_data_file(), c(input$r1, input$c1))
    
  })
  
  #server on importing metabolite target list
  metabolite_target_file <- reactive({
    req(input$metabolite_target_file)
    
    ext <- tools::file_ext(input$metabolite_target_file$name)
    switch(ext,
           csv = vroom::vroom(input$metabolite_target_file$datapath, delim = ","),
           tsv = vroom::vroom(input$metabolite_target_file$datapath, delim = "\t"),
           validate("Invalid file; Please upload a .csv or .tsv file")
    )
  })
  
  output$metabolite_target_file <- renderTable({
    head(metabolite_target_file(), input$r2)
    
  })

#server on importing metadata
metabolite_metadata_file <- reactive({
  req(input$metabolite_metadata_file)
  
  ext <- tools::file_ext(input$metabolite_metadata_file$name)
  switch(ext,
         csv = vroom::vroom(input$metabolite_metadata_file$datapath, delim = ","),
         tsv = vroom::vroom(input$metabolite_metadata_file$datapath, delim = "\t"),
         validate("Invalid file; Please upload a .csv or .tsv file")
  )
})

output$metabolite_metadata_file <- renderTable({
  head(metabolite_metadata_file(), c(input$r3, input$c3))
  
})


# server for boxplot ################################################

# data <- reactive({
#   inFile <- input$metabolite_data_file
#   if(!is.null(inFile)){
#     read.csv(inFile$datapath, header = input$header, stringsAsFactors = FALSE)    
#   }
# })
# 
# output$plot1 <- renderPlot({
#   req(data())
#   print("help") 
# })

plot_data <- reactive({
  if (is.null(input$metabolite_data_file)) return(NULL)
  read.csv(input$metabolite_data_file$datapath, header = input$header)
})

output$contents <- renderTable({
  req(plot_data()) # if user_data() is null than the reactive will not run
  head(plot_data(), c(input$r1, input$c1))
})

output$plot1 <- renderPlot({
  #req(user_data())
  #plot(user_data$speed, user_data$dist)
})




})#final bracket for server

# shinyApp()
shinyApp(ui = ui, server = server)


#branch_test
