library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyFiles)
#library(data.table)
library(DT)
library(stringr)
library(dplyr)

# Define UI for application that draws a histogram
ui <- navbarPage("FACkit Analysis",
                 mainPanel(), #includeMarkdown(path="front_page_notes.md") ## insert into main panel for documentation
                 tabPanel("Data Import",
                          mainPanel(
                            h5("Select Data Folder"),
                            fileInput("file1", "Choose CSV File",multiple = TRUE,
                                      accept = c(
                                        "text/csv",
                                        "text/comma-separated-values,text/plain",
                                        ".csv")
                            ),
                            numericInput(inputId="n.cond.cols", label="Number of Condition Columns", value=3, min=0, max=Inf, step=1, width="15%"),

                            h5("Check column names:"),
                            h4("output$dir"),
                            verbatimTextOutput("dir"), br(),
                            h4("Files in that dir"),
                            verbatimTextOutput("files"),
                            h4("table"),
                            DTOutput("table")
                          )
                 )
)

server <- function(input, output) {
  data.folder <- reactiveValues()



  observe({
    files <- input$file1$name
    if(is.null(files)){return(NULL)}
    n.cols <- input$n.cond.cols
    if(n.cols > 0)
    {
      data.files <- matrix(nrow=length(files), ncol=1+n.cols,
                           dimnames=list(c(1:length(files)), c("Data.Files",paste("Cond",1:n.cols, sep="."))))
    }
    else{
      data.files <- matrix(nrow=length(files), ncol=1+n.cols,
                           dimnames=list(c(1:length(files)), c("Data.Files")))
    }
    data.files <- as.data.frame(data.files, stringsAsFactors=FALSE)
    data.files$Data.Files <- files
    data.files[] <- lapply(data.files, as.character)
    data.folder$files <- data.files
  })

  output$table <- renderDT(data.folder[["files"]], editable = TRUE, selection = "none", server = TRUE)

  proxy = dataTableProxy("table")
  observeEvent(input$table_cell_edit, {
    info = input$table_cell_edit
    str(info)
    i = info$row
    j = info$col
    v = info$value
    data.folder[["files"]][i, j] <<- DT::coerceValue(v, data.folder[["files"]][i, j])
    replaceData(proxy, data.folder[["files"]],resetPaging = FALSE)  # important
    str(data.folder[["files"]])
    data.folder[["files"]] %>% as_data_frame %>% print
  })

  output$files <- renderPrint(data.folder[["files"]])

}

# Run the application
shinyApp(ui = ui, server = server)

