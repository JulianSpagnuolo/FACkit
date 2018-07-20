library(shiny)
library(shinydashboard)
library(shinyjs)
library(shinyFiles)
library(data.table)
library(DT)
library(stringr)
library(dplyr)

library(ggplot2)
library(ggalt)
library(ggbeeswarm)
library(gridExtra)
library(pals)
library(ggthemes)

library(plotly)

library(FACkit)


## TODO Make the layout pretty with shinydashboards.
## TODO Add documentation guides in the boxes
## TODO implement download of workspace
## TODO implement download of figures
## TODO implement dynamic report creation and download
ui <- dashboardPage(dashboardHeader(title="FACkit Analysis"),
                    dashboardSidebar(sidebarMenu(menuItem("Home", tabName = "home", icon=icon("home", lib = "font-awesome"), selected=TRUE),
                                                 menuItem("Data Import", tabName = "data_import", icon=icon("import", lib = "glyphicon")),
                                                 menuItem("Transformation", tabName = "transformation", icon=icon("equalizer", lib = "glyphicon")))),
                 dashboardBody(
                   tabItems(
                     tabItem(tabName = "home",
                           fluidRow(
                             #includeMarkdown(path="front_page_notes.md") ## insert into main panel for documentation
                             )
                           ),
                   tabItem(tabName = "data_import",
                           h2("Data Import"),
                           fluidRow(h4("Select Data Folder"),
                                    column(fileInput("file1", "Choose CSV File",multiple = TRUE,
                                                     accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
                                           width=4),
                                    column(numericInput(inputId="n.cond.cols", label="Number of Condition Columns", value=3, min=0, max=Inf, step=1, width="25%"),
                                           width=8)
                                    ),
                           fluidRow(h4("Add Conditional Data:"),
                                    DTOutput("table")),
                           fluidRow(h4("Check Marker Names:"),
                                    DTOutput("column.names")),
                           fluidRow(actionButton(inputId = "upload", label = "Upload Data", icon = icon("upload", lib="font-awesome")))
                   ),
                   tabItem(tabName = "transformation",
                           h2("Data Transformation"),
                           ## TODO Add method for tabulating the conditional data columns
                           ## TODO add downsampling methods
                           h4("Define Cutoffs"),
                           ## TODO add option to define cutoffs from scatter plot
                           fluidRow(box(plotOutput("raw.dist", click = "cutoff.click")),
                                    box(selectInput("raw.dist.marker", label = h5("Select Marker"),
                                                    choices = list(),
                                                    selected = 1),
                                        sliderInput("raw.dist.range", label = h5("Plot Range"),
                                                    min = -100, max = 100, value = c(30,40)))),
                           h4("Cutoff Values"),
                           fluidRow(DTOutput("cutoffs")),
                           h4("Run Transformation:"),
                           numericInput(inputId = "asincofac", label = h5("Arc Sin Cofactor"), value = 25, min=0, max = Inf,width = "25%"),
                           actionButton("transform", label = "Apply Transformation", icon = icon("ok",lib="glyphicon")),
                           h3("Check Transformed Data:"),
                           fluidRow(box(plotOutput("norm.dist", click = "bin.click")),
                                    box(selectInput("norm.dist.marker", label = h5("Select Marker"),
                                                    choices = list(), selected = 1)
                                                 )),
                           ## TODO Implement bin defs for scatter plot??
                           fluidRow(box(plotlyOutput("norm.scat")),
                                    box(selectInput("norm.scat.x", label=h5("Select X-axis Marker"),
                                                    choices = list(), selected = 1),
                                        selectInput("norm.scat.y", label=h5("Select Y-axis Marker"),
                                                    choices = list(), selected = 1))),
                           br(),
                           h3("Marker Enrichment Modelling"),
                           h4("Choose Markers to Include In Analysis:"),
                           fluidRow(uiOutput("marker.select")),
                           ## TODO implement the MEM step
                           h4("MEM"),
                           fluidRow(box())
                           )
                 )
                 )
)


server <- function(input, output, session) {
  # TODO Modularise the shiny app.
  data.folder <- reactiveValues()
  expdata <- reactiveValues()
  # input$file1 will be NULL initially. After the user selects
  # and uploads a file, it will be a data frame with 'name',
  # 'size', 'type', and 'datapath' columns. The 'datapath'
  # column will contain the local filenames where the data can
  # be found.

  # Add Conditional Columns
  observe({
    ## TODO add method to give custom Cond column names
    ## TODO add method to store cond column names (if present)
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
  })


  # Check Column Names
  observe({
    path <- input$file1$datapath
    if(is.null(path)){return(NULL)}
    header <- vector(length=length(path), mode="list")
    str(path)
    for(i in 1:length(path))
    {
      header[[i]] <- read.csv(file = path[i],nrows = 1, header = FALSE, stringsAsFactors = FALSE)
    }
    if(length(unique(unlist(lapply(header, length)))) == 1)
    {
      col.names  <- as.data.frame(matrix(nrow=length(path), ncol=length(header[[1]]), data=unlist(header), byrow=TRUE), stringsAsFactors=FALSE)
    }
    else{
          col.names <- matrix(nrow=length(path), ncol=max(unlist(lapply(header, length))))
          for(i in 1:length(header))
          {
            if(length(header[[i]]) == ncol(col.names))
            {
              col.names[i,] <- header[[i]]
            }
            else{
                  col.names[i,] <- c(header[[i]], rep(NA, ncol(col.names)-length(header[[i]])))
                }
          }
          col.names <- as.data.frame(col.names, stringsAsFactors = FALSE)
    }
    # check for identical column names, return single row matrix if TRUE.
    if (all(apply(col.names, MARGIN = 2, FUN = function(x){length(unique(x))}) == 1))
    {
      col.names <- col.names[1,]
    }
    data.folder$col.names <- col.names

    rm(col.names);rm(header)
  })

  output$column.names <- renderDT(data.folder[["col.names"]], editable = TRUE, selection = "none", server = TRUE)
  proxy.cols = dataTableProxy("column.names")
  observeEvent(input$column.names_cell_edit, {
    info = input$column.names_cell_edit
    i = info$row
    j = info$col
    v = info$value
    data.folder[["col.names"]][i, j] <<- DT::coerceValue(v, data.folder[["col.names"]][i, j])
    replaceData(proxy.cols, data.folder[["col.names"]],resetPaging = FALSE)  # important
    str(data.folder[["col.names"]])

  })

  ## Data Upload
  observeEvent(input$upload,{
    paths <- input$file1$datapath

    # Create Matrix of Column Names if it was reduced to one row (b/c they were all the same)
    if(nrow(data.folder[["col.names"]]) == 1)
    {
      column.names <- vector(mode="list", length=length(paths))
      for(i in 1:length(column.names))
      {
        column.names[[i]] <- data.folder[["col.names"]]
      }
      column.names <- do.call("rbind", column.names)
      column.names <- as.data.frame(column.names, stringsAsFactors=FALSE)
    }
    else{
      column.names <- as.data.frame(data.folder[["col.names"]], stringsAsFactors=FALSE)
    }
    str(column.names[1,])
    n.cols <- input$n.cond.cols
    raw.data.list <- vector(length=length(paths), mode="list")
    metadata.list <- vector(length=length(paths), mode="list")
    for(i in 1:length(paths))
    {
       raw.data <- fread(file = paths[i], col.names = as.character(column.names[i,]))
       if(n.cols > 0)
       {
         ## TODO make this function handle custom column names if present
         metadata <- matrix(nrow=nrow(raw.data), ncol=n.cols,
                            dimnames=list(c(1:nrow(raw.data)), c(paste("Cond",1:n.cols, sep="."))))
         for(n in colnames(metadata))
         {
           metadata[,n] <- as.character(data.folder[["files"]][i,n])
         }
         metadata.list[[data.folder[["files"]][i,"Data.Files"]]] <- as.data.frame(metadata, stringsAsFactors=FALSE)
       }
       raw.data.list[[data.folder[["files"]][i,"Data.Files"]]] <- as.data.frame(raw.data)
    }
    raw.data <- as.data.frame(do.call("rbind", raw.data.list), stringsAsFactors=FALSE)
    if(n.cols > 0)
    {
      metadata <- as.data.frame(do.call("rbind", metadata.list), stringsAsFactors=FALSE)
      expdata$markers.raw <- colnames(raw.data)
      expdata$metadata <- colnames(metadata)

      raw.data <- cbind(raw.data, metadata)
      expdata$raw.data <- raw.data
    }
    else{
      expdata$raw.data <- raw.data
      expdata$markers.raw <- colnames(raw.data)
    }

    updateSelectInput(session, "raw.dist.marker",
                      choices = as.vector(expdata[["markers.raw"]]),
                      selected = as.vector(expdata[["markers.raw"]])[1])

    updateSelectInput(session, "norm.dist.marker",
                      choices = as.vector(expdata[["markers.raw"]]),
                      selected = as.vector(expdata[["markers.raw"]])[1])

    updateSelectInput(session, "norm.scat.x",
                      choices = as.vector(expdata[["markers.raw"]]),
                      selected = as.vector(expdata[["markers.raw"]])[1])
    updateSelectInput(session, "norm.scat.y",
                      choices = as.vector(expdata[["markers.raw"]]),
                      selected = as.vector(expdata[["markers.raw"]])[1])

    expdata$cutoffs <- matrix(nrow=1, ncol = length(as.vector(expdata[["markers.raw"]])), dimnames=list(c(1),c(as.vector(expdata[["markers.raw"]]))))
  })


  # Data Transformation

  ## Choose Marker to Plot Density Distribution of Raw Expression Values
  observeEvent(input$raw.dist.marker,{
    marker <- input$raw.dist.marker
    updateSliderInput(session, "raw.dist.range",
                      min = min(expdata[["raw.data"]][,marker]),
                      max = max(expdata[["raw.data"]][,marker]), step = 1,
                      value = c(min(expdata[["raw.data"]][,marker]),
                                max(expdata[["raw.data"]][,marker])))
    data.folder[["cut"]] <- 0
  })
  output$raw.dist <- renderPlot({
    ggplot(expdata[["raw.data"]], aes_string(x=input$raw.dist.marker)) +geom_bkde() +geom_vline(xintercept=as.vector(data.folder[["cut"]]), colour="red") +scale_x_continuous(limits=input$raw.dist.range)
  })

  ## Define Cutoff Value from Plot
  observeEvent(input$cutoff.click,{
    data.folder$cut <- input$cutoff.click$x
    expdata[["cutoffs"]][,input$raw.dist.marker] <- input$cutoff.click$x
  })

  ## Define Cutoff Value by Manually Entering in DataTable
  output$cutoffs <- renderDT(expdata[["cutoffs"]], editable = TRUE, selection = "none", server = TRUE)
  proxy.cols = dataTableProxy("cutoffs")
  observeEvent(input$cutoffs_cell_edit, {
    info = input$cutoffs_cell_edit
    i = info$row
    j = info$col
    v = as.numeric(info$value)
    expdata[["cutoffs"]][i, j] <<- DT::coerceValue(v, expdata[["cutoffs"]][i, j])
    replaceData(proxy.cols, expdata[["cutoffs"]],resetPaging = FALSE)  # important
    str(expdata[["cutoffs"]])
  })

  ## Run Transformation
  observeEvent(input$transform,{
    norm.data <- facsnorm(x=expdata[["raw.data"]][,colnames(expdata[["cutoffs"]])], cutoffs = as.numeric(expdata[["cutoffs"]][1,]), asinCofac = input$asincofac, method = "arcsin")

    norm.data <- cbind(norm.data, expdata[["raw.data"]][,expdata[["metadata"]]])
    expdata$norm.data <- norm.data
    expdata$bin.defs <- matrix(nrow=2, ncol=length(expdata[["markers.raw"]]), dimnames = list(c("pos","neg"),c(expdata[["markers.raw"]])))
  })

  ## Check Transformed Data
  output$norm.dist <- renderPlot({
    if(!is.numeric(expdata[["bin.defs"]][,input$norm.dist.marker])){
      ggplot(expdata[["norm.data"]], aes_string(x=input$norm.dist.marker)) +geom_bkde()
    }else{
      ggplot(expdata[["norm.data"]], aes_string(x=input$norm.dist.marker)) +geom_bkde() +geom_vline(xintercept = expdata[["bin.defs"]][,input$norm.dist.marker], colour=c("red","red"))
    }

  })

  output$norm.scat <- renderPlotly({
    plot_ly(x = expdata[["norm.data"]][,input$norm.scat.x], y = expdata[["norm.data"]][,input$norm.scat.y], type="scattergl", mode = "markers",
            hoverinfo="none", marker = list(size = 3, color = 'rgba(0, 0, 0, .5)')) %>%
      layout(yaxis = list(title=input$norm.scat.y),
             xaxis = list(title=input$norm.scat.x))

  })

  observeEvent(input$bin.click,{
    if(sign(input$bin.click$x) == 1){
      expdata[["bin.defs"]]["pos",input$norm.dist.marker] <- input$bin.click$x
    }
    if(sign(input$bin.click$x) == -1){
      expdata[["bin.defs"]]["neg",input$norm.dist.marker] <- input$bin.click$x
    }
  })

  output$marker.select <- renderUI({
    checkboxGroupInput(inputId = "marker.select",
                       label = "Select Markers",
                       choices = expdata[["markers.raw"]], inline = TRUE)
  })
  observe({
    expdata$markers.norm <- input$marker.select ## input$marker.select is a character vector of markers selected in the checkbox.
  })

}

# Run the application
shinyApp(ui = ui, server = server)

