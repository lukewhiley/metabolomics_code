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
  sidebarPanel(
    p("select a .csv or .tsv containing metabolite concentration values. Column headers should contain the metabolite name matching the metabolite list below. sample/filenames should be in a column called: sample name"),
    fileInput("metabolite_data_file", NULL, accept = c(".csv", ".tsv")),
    numericInput("r", "Rows to preview", 10, min = 1),
    numericInput("c", "Columns to preview", 5, min = 1)
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
    numericInput("r", "Rows to preview", 5, min = 1),
  ),
  mainPanel(
    h3("Metabolite targets"),
    tableOutput("metabolite_target_file")
  )
)

ui <- fluidPage(
  ui_title,
  ui_data_upload,
  ui_metabolite_target_upload
)


# server()
server <- function(input, output, session) {
  
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
    head(metabolite_data_file(), c(input$r, input$c))
    
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
    head(metabolite_target_file(), input$r)
    
  })
}

# shinyApp()
shinyApp(ui = ui, server = server)

