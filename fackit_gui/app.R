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
library(RColorBrewer)

library(plotly)

library(pheatmap)

library(FACkit)
library(MEM) ## TODO Add MEM to FACkit properly. Rewrite the function to make it more compatible and better.

# load the faster fftRtsne if it exists in user lib, else load the standard Rtsne
if(length(find.package("fftRtsne", quiet = T)) != 0){
  library(fftRtsne)
}else{
  library(Rtsne)
}

## TODO MEM: allow use of character vectors as the cluster column
## TODO MEM: allow a vector of column names to choose the desired markeres instead of shitty interactive mode
## TODO MEM: alter plotting function... or change to use pheatmap.... something.


## TODO Add method to detect fftwRtsne package and load it instead of Rtsne.

## TODO Make the layout pretty with shinydashboards.
## TODO Add documentation guides in the boxes
## TODO implement download of workspace
## TODO implement download of figures
## TODO implement dynamic report creation and download
ui <- dashboardPage(dashboardHeader(title="FACkit Analysis"),
                    dashboardSidebar(sidebarMenu(menuItem("Home", tabName = "home", icon=icon("home", lib = "font-awesome"), selected=TRUE),
                                                 menuItem("Data Import", tabName = "data_import", icon=icon("import", lib = "glyphicon")),
                                                 menuItem("Transformation", tabName = "transformation", icon=icon("equalizer", lib = "glyphicon")),
                                                 menuItem("Dimensional Reduction", tabName = "dim_red", icon=icon("sitemap", lib = "font-awesome")))
                    ),
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
                                fluidRow(
                                  box(h4("Set MEM Parameters"),
                                      uiOutput("marker.select"),
                                      uiOutput("mem.groups"),
                                      numericInput("mem.iqr", label = h5("IQR Threshold"), value = NULL, width = "25%"),
                                      actionButton("run.mem", label="Run", icon = icon("magic", lib="font-awesome")),
                                      width = 12
                                  )
                                ),
                                ## TODO implement the MEM step
                                h4("MEM Results"),
                                fluidRow(box(plotOutput("mem.median"), title = "MEM Group Median Expression"),
                                         box(plotOutput("mem.results"), title="MEM Marker Enrichment Score"))
                        ),
                        tabItem(tabName = "dim_red",
                                h2("tSNE"),
                                fluidRow(
                                  h3("Dimensional Reduction - tSNE"),
                                  uiOutput("tsne.ui")
                                ),
                                fluidRow(
                                  box(plotlyOutput(outputId = "tsne.plot", width = "100%", height = "100%"), width = 8, height = "800px"),
                                  box(title = "tSNE Plot Parameters", width = 4,
                                      selectInput("tsne.col", label = "Colour Variable", choices = list(), selected=1, multiple = FALSE))
                                )
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
    ## TODO add method to add numeric/integer columns - auto detect entered data
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
    rownames(raw.data) <- 1:nrow(raw.data)
    if(n.cols > 0)
    {
      metadata <- as.data.frame(do.call("rbind", metadata.list), stringsAsFactors=FALSE)
      rownames(metadata) <- 1:nrow(metadata)
      expdata$markers.raw <- colnames(raw.data)
      expdata$metadata <- colnames(metadata)

      raw.data <- cbind(raw.data, metadata)
      rownames(raw.data) <- 1:nrow(raw.data)
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

    metadata <- matrix(nrow=nrow(norm.data), ncol=length(expdata[["metadata"]]), dimnames=list(c(1:nrow(norm.data)), c(expdata[["metadata"]])),
                       data=as.vector(expdata[["raw.data"]][,expdata[["metadata"]]]), byrow = TRUE)

    norm.data <- cbind(norm.data, metadata)
    rownames(norm.data) <- 1:nrow(norm.data)
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

  # Define Population Bins
  observeEvent(input$bin.click,{
    if(sign(input$bin.click$x) == 1){
      expdata[["bin.defs"]]["pos",input$norm.dist.marker] <- input$bin.click$x
    }
    if(sign(input$bin.click$x) == -1){
      expdata[["bin.defs"]]["neg",input$norm.dist.marker] <- input$bin.click$x
    }
  })

  # Choose Markers for MEM
  output$marker.select <- renderUI({
    checkboxGroupInput(inputId = "marker.select",label = "Select Markers", inline = TRUE,
                       choices = expdata[["markers.raw"]], selected = expdata[["markers.raw"]])
  })
  observe({
    expdata$markers.mem <- input$marker.select ## input$marker.select is a character vector of markers selected in the checkbox.
  })

  # Choose Conditional Column for MEM
  output$mem.groups <- renderUI({
    selectInput(inputId = "mem.groups", label = "Select Conditions", multiple = FALSE, width = "25%",
                choices = expdata[["metadata"]])
  })

  # Run MEM and plot
  observeEvent(input$run.mem, {
    c("formatting data for MEM") %>% print

    mem.data <- expdata[["norm.data"]][,c(expdata[["markers.mem"]])] ## TODO once MEM function is rewritten, replace with: expdata[["norm.data"]][,c(expdata[["tsne.markers"]], input$mem.groups)]
    c("Collected Expression Data, now adding clust column") %>% print
    input$mem.groups %>% print
    colnames(mem.data) %>% print
    head(expdata[["norm.data"]][,c(input$mem.groups)]) %>% print ### BUG This appears to output a list object ... why?!?!
    table(is.na(expdata[["norm.data"]][,c(input$mem.groups)]))
    ## BUG this step is taking way too long for some reason when there are more than one cond column
    mem.data$cluster <- as.numeric(as.factor(expdata[["norm.data"]][,c(input$mem.groups)])) ## TODO once MEM function is rewritten change this so that the expdata is directly input into the call to MEM without changing cluster to a factor and numeric... etc

    head(mem.data) %>% print
    c("Running MEM") %>% print
    if(is.na(input$mem.iqr))
    {
      expdata$mem.res <- MEM(exp_data = mem.data, transform = FALSE, choose.markers = FALSE, choose.ref = FALSE, IQR_thresh = "auto")
    }else{
      expdata$mem.res <- MEM(exp_data = mem.data, transform = FALSE, choose.markers = FALSE, choose.ref = FALSE, IQR_thresh = input$mem.iqr)
    }

    c("Plotting MEM") %>% print
    output$mem.median <- renderPlot({
      pheatmap(mat=expdata[["mem.res"]][["MAGpop"]][[1]], border_color = NA, cluster_rows = TRUE, cluster_cols = TRUE,
               clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean")
    })
    output$mem.results <- renderPlot({
      pheatmap(mat=expdata[["mem.res"]][["MEM_matrix"]][[1]], border_color = NA, cluster_rows = TRUE, cluster_cols = TRUE,
               clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean")
    })
  })

  ## define tsne UI based on available pacakge (fftRtsne or Rtsne)
  output$tsne.ui <- renderUI({
    if(length(find.package("fftRtsne", quiet = T)) != 0)
    {
      box(
        h4("Set tSNE Parameters:"),
        checkboxGroupInput(inputId = "tsne.markers",label = "Select Markers", inline = TRUE,
                           choices = expdata[["markers.raw"]], selected = expdata[["markers.raw"]]),
        numericInput(inputId = "tsne.dim", label="tSNE Dimensions", value=2, min=1, step=1, width = "25%"),
        numericInput(inputId = "tsne.perp", label="Perplexity", value=30, min=1, step=1, width="25%"),
        numericInput(inputId = "tsne.iter", label="Max Iterations", value=1000, min=1, step=1, width="25%"),
        radioButtons(inputId = "tsne.mode", label = "tSNE Mode", inline = TRUE, c("FFT"=TRUE, "BH"=FALSE), selected = "FFT"),
        radioButtons(inputId = "tsne.tree", label = "tSNE NN Mode", inline = TRUE, c("ANNOY"=TRUE, "Vantage-Point"=FALSE), selected = "Vantage-Point"),
        numericInput(inputId = "tsne.stop.lying.iter", label="Early Exagg. Phase End", value = 250, min = 1, step=1, width="25%"),
        numericInput(inputId = "tsne.early.exag", label="Early Exagg. Coefficient", value = 12.0, min = 0, step=0.5, width="25%"),
        numericInput(inputId = "tsne.start.late.exag", label="Start Late Exag. At Iter:", value = -1, min = -1, step=1, width="25%"),
        numericInput(inputId = "tsne.late.exag", label="Late Exagg. Coefficient", value = 1.0, min = 0, step=0.5, width="25%"),
        numericInput(inputId = "seed", label="System Seed", value = 42, min = 0, step=1, width="25%"),
        actionButton(inputId = "run.tsne", label = "Run tSNE", icon = icon("magic", lib="font-awesome")), width=12
      )

    }else{
      box(
        h4("Set tSNE Parameters:"),
        checkboxGroupInput(inputId = "tsne.markers",label = "Select Markers", inline = TRUE,
                           choices = expdata[["markers.raw"]], selected = expdata[["markers.raw"]]),
        numericInput(inputId = "tsne.dim", label="tSNE Dimensions", value=2, min=1, step=1, width = "25%"),
        numericInput(inputId = "tsne.perp", label="Perplexity", value=30, min=1, step=1, width="25%"),
        numericInput(inputId = "tsne.iter", label="Max Iterations", value=1000, min=1, step=1, width="25%"),
        numericInput(inputId = "tsne.stop.lying.iter", label="Early Exagg. Phase End", value = 250, min = 1, step=1, width="25%"),
        numericInput(inputId = "tsne.early.exag", label="Early Exagg. Coefficient", value = 12.0, min = 0, step=0.5, width="25%"),
        numericInput(inputId = "seed", label="System Seed", value = 42, min = 0, step=1, width="25%"),
        actionButton(inputId = "run.tsne", label = "Run tSNE", icon = icon("magic", lib="font-awesome"))
      )

    }
  })

  observeEvent(input$run.tsne,{
    expdata$tsne.markers <- input$tsne.markers
    updateSelectInput(session, "tsne.col",
                      choices = as.vector(c(expdata[["tsne.markers"]], expdata[["metadata"]])),
                      selected = as.vector(c(expdata[["tsne.markers"]], expdata[["metadata"]])))

    set.seed(input$seed)

    if(length(find.package("fftRtsne", quiet = T)) != 0){
      "fftRtsne" %>% print
      tsne <- fftRtsne(X = expdata[["norm.data"]][,c(expdata[["tsne.markers"]])],
                       dims = input$tsne.dim, perplexity = input$tsne.perp, check_duplicates = FALSE, max_iter = input$tsne.iter,
                       fft_not_bh = input$tsne.mode, ann_not_vptree = input$tsne.tree, stop_lying_iter = input$tsne.stop.lying.iter,
                       exaggeration_factor = input$tsne.early.exag, no_momentum_during_exag = FALSE, start_late_exag_iter = input$tsne.start.late.exag,
                       late_exag_coeff = input$tsne.late.exag, rand_seed = input$seed)

      expdata$tsne <- data.frame(tsne1=tsne[,1], tsne2=tsne[,2])

    }else{
      "Rtsne" %>% print
      tsne <- Rtsne(X = expdata[["norm.data"]][,c(expdata[["tsne.markers"]])], dims = input$tsne.dim, perplexity = input$tsne.perp, check_duplicates = FALSE, max_iter = input$tsne.iter,
                    stop_lying_iter = input$tsne.stop.lying.iter, exaggeration_factor = input$tsne.early.exag)

      expdata$tsne <- data.frame(tsne1=tsne$Y[,1], tsne2=tsne$Y[,2])
    }



    "finished tsne" %>% print
    output$tsne.plot <- renderPlotly({
      ## TODO Make plotly output a figure with 1:1 aspect ratio!
        plot_ly(x = expdata[["tsne"]][,1], y = expdata[["tsne"]][,2], alpha = 0.5, type="scattergl", mode = "markers",
                hoverinfo="none", marker = list(size = 3)) %>%
          layout(xaxis = list(title="tSNE-1"),
                 yaxis = list(title="tSNE-2"))
    })

  })

  observe({
    output$tsne.plot <- renderPlotly({
      ## TODO Make plotly output a figure with 1:1 aspect ratio!

      if(input$tsne.col %in% expdata[["tsne.markers"]]){
        plot_ly(x = expdata[["tsne"]][,1], y = expdata[["tsne"]][,2], colors = viridis(100), alpha = 0.5, color = expdata[["norm.data"]][,input$tsne.col], type="scattergl", mode = "markers",
                hoverinfo="none", marker = list(size = 3)) %>%
          layout(xaxis = list(title="tSNE-1"),
                 yaxis = list(title="tSNE-2"))
      }else{
        if(length(unique(expdata[["norm.data"]][,input$tsne.col])) <= 12){
          plot_ly(x = expdata[["tsne"]][,1], y = expdata[["tsne"]][,2], colors = colorblind_pal()(length(unique(expdata[["norm.data"]][,input$tsne.col]))), alpha = 0.5,
                  color = expdata[["norm.data"]][,input$tsne.col], type="scattergl", mode = "markers",
                  hoverinfo="none", marker = list(size = 3)) %>%
            layout(xaxis = list(title="tSNE-1"),
                   yaxis = list(title="tSNE-2"))
        }else{
          qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
          col_vector <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

          plot_ly(x = expdata[["tsne"]][,1], y = expdata[["tsne"]][,2], colors = col_vector[sample.int(n = length(col_vector), size = length(unique(expdata[["norm.data"]][,input$tsne.col])), replace = F)], alpha = 0.5,
                  color = expdata[["norm.data"]][,input$tsne.col], type="scattergl", mode = "markers",
                  hoverinfo="none", marker = list(size = 3)) %>%
            layout(xaxis = list(title="tSNE-1"),
                   yaxis = list(title="tSNE-2"))
        }

      }
    })
  })


}

# Run the application
shinyApp(ui = ui, server = server)

