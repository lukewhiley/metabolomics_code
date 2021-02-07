library(shiny)
library(png)
library(vroom)

# ui object
ui_title <- titlePanel(title=div(img(src="ANPC_logo.png", width = "240px", height = "140px"), 
                                   "LGW - Metabolomics Workbook", style = "color:#3474A7"))
  
ui_upload <- sidebarLayout(
  sidebarPanel(
    fileInput("file", NULL, accept = c(".csv", ".tsv")),
    numericInput("r", "Rows to preview", 10, min = 1),
    numericInput("c", "Columns to preview", 5, min = 1)
  ),
  mainPanel(
    h3("Raw data"),
    tableOutput("head")
  )
)

ui <- fluidPage(
  ui_title,
  ui_upload
)


# server()
server <- function(input, output, session) {
  data <- reactive({
    req(input$file)
    
    ext <- tools::file_ext(input$file$name)
    switch(ext,
           csv = vroom::vroom(input$file$datapath, delim = ","),
           tsv = vroom::vroom(input$file$datapath, delim = "\t"),
           validate("Invalid file; Please upload a .csv or .tsv file")
    )
  })
  
  output$head <- renderTable({
    head(data(), c(input$r, input$c))
    
  })
}

# shinyApp()
shinyApp(ui = ui, server = server)

