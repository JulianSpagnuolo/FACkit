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

library(vipor)
library(plotly)

library(pheatmap)

library(FACkit)
library(MEM)
## TODO Add MEM to FACkit properly. Rewrite the function to make it more compatible and better.
## TODO MEM: allow use of character vectors as the cluster column
## TODO MEM: allow a vector of column names to choose the desired markeres instead of shitty interactive mode
## TODO MEM: alter plotting function... or change to use pheatmap.... something.

# load the faster fftRtsne if it exists in user lib, else load the standard Rtsne
if(length(find.package("fftRtsne", quiet = T)) != 0){
  library(fftRtsne)
}else{
  library(Rtsne)
}
find.package("FACkit")
source("~/Documents/Rproj/fackit/fackit_gui/enrichTest-module.R")

## TODO Make the layout pretty with shinydashboards.
## TODO Add documentation guides in the boxes
## TODO implement download of workspace
## TODO implement download of figures
## TODO implement dynamic report creation and download
## TODO Fix all datatables sig figs and layout options - currently is not working.
## TODO implement server side bookmarking of state, including the associated datafiles.
## TODO implement dynamic report rendering from results.
ui <- dashboardPage(skin = "blue",
                    dashboardHeader(title="FACkit Analysis"),
                    dashboardSidebar(sidebarMenu(menuItem("Home", tabName = "home", icon=icon("home", lib = "font-awesome"), selected=TRUE),
                                                 menuItem("Data Import", tabName = "data_import", icon=icon("import", lib = "glyphicon")),
                                                 menuItem("Transformation", tabName = "transformation", icon=icon("equalizer", lib = "glyphicon")),
                                                 menuItem("Dimensional Reduction", tabName = "dim_red", icon=icon("sitemap", lib = "font-awesome")),
                                                 menuItem("Cluster Enrichment", tabName = "clust_enrich", icon=icon("search", lib = "font-awesome")),
                                                 downloadButton(outputId = "fackit.download", label="Download Data", icon=icon("download", lib="font-awesome")) # TODO centre the button in the sidebar, make it pretty!
                                                 )
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
                                fluidRow(box(width = 12, title = h4("Select Data Folder"),
                                             column(fileInput("file1", "Choose CSV File", multiple = TRUE,
                                                              accept = c("text/csv","text/comma-separated-values,text/plain",".csv")),
                                                    width=4),
                                             column(numericInput(inputId="n.cond.cols", label="Number of Condition Columns", value=3, min=0, max=Inf, step=1, width="25%"),
                                                    width=8)
                                             )
                                ),
                                fluidRow(
                                  box(width = 12, title = h4("Add Conditional Data:"),
                                      DTOutput("table")
                                      ),
                                  box(width = 12, title = h4("Check Marker Names:"),
                                      DTOutput("column.names")
                                      )
                                  ),

                                fluidRow(
                                  box(
                                    actionButton(inputId = "upload", label = "Upload Data", icon = icon("upload", lib="font-awesome"))
                                    )
                                  )
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
                                fluidRow(
                                  box(width = 12, title = h4("Cutoff Values"),
                                      DTOutput("cutoffs")
                                      )
                                  ),
                                fluidRow(
                                  ### TODO Move to box with cutoffs DT
                                  box(width = 12, title = h4("Run Transformation:"),
                                      column(width = 6,
                                             numericInput(inputId = "asincofac", label = h5("Arc Sin Cofactor"), value = 25, min=0, max = Inf,width = "25%")
                                             ),
                                      column(width = 6,
                                             actionButton("transform", label = "Apply Transformation", icon = icon("ok",lib="glyphicon"))
                                             )
                                      )
                                ),
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
                                  box(width = 12,
                                      h4("Set MEM Parameters"),
                                      fluidRow(uiOutput("marker.select")), ## TODO this element needs to move right a bit.
                                      fluidRow(
                                        column(width = 6,
                                               uiOutput("mem.groups"),
                                               actionButton("run.mem", label="Run", icon = icon("magic", lib="font-awesome"))
                                               ),
                                        column(width = 6,
                                               numericInput("mem.iqr", label = "IQR Threshold", value = NULL, width = "25%"))
                                        )
                                      )
                                ),
                                h4("MEM Results"),
                                fluidRow(box(plotOutput("mem.median"), title = "MEM Group Median Expression"),
                                         box(plotOutput("mem.results"), title="MEM Marker Enrichment Score"))
                        ),
                        tabItem(tabName = "dim_red",
                                h2("tSNE"),
                                fluidRow(
                                  uiOutput("tsne.ui")
                                  ),
                                fluidRow(
                                  column(width=7,
                                         box(plotlyOutput(outputId = "tsne.plot"), width = 12, height="700px")
                                  ),
                                  column(width=5,
                                         box(title = "tSNE Plot Parameters", width = 12,
                                             fluidRow(
                                               column(width=6,
                                                      selectInput("tsne.col", label = "Colour Variable", choices = list(), selected=1, multiple = FALSE, width = "60%")),
                                               column(width=6,
                                                      radioButtons("db.tsne.dim", label = "tSNE Dimensions To Use:", choices = c("2D"="tsne", "1D"="tsne1d"), inline = TRUE))
                                               ),
                                             fluidRow(
                                               column(width=6,
                                                      actionButton("db.knn.run", label = "Run DBScan kNN dist", icon=icon("magic")))
                                             )
                                         ),
                                         box(plotOutput("db.knn.plot"), width = 12)
                                  )
                                ),
                                h2("DBscan Clustering"), ## TODO Implement DBScan clustering
                                fluidRow( ## TODO Implement cluster refinement
                                  box(plotOutput(outputId = "db.opt.plot"), width = 8),
                                  box(title = "DBscan Parameter Scanning", width = 4,
                                      column(width = 6,
                                             numericInput("db.opt.eps.start", label="Epsilon Start", min = 0, value = 0.01, width="50%"),
                                             numericInput("db.opt.eps.step", label="Epsilon Step Size", min = 0, value = 0.001, step = 0.001, width="50%"),
                                             numericInput("db.opt.mpts.start", label = "Min Points Start", min = 1, value = 3, step = 1, width="50%"),
                                             actionButton("db.opt.run", label = "Run", icon = icon("magic", lib="font-awesome"))
                                      ),
                                      column(width = 6,
                                             numericInput("db.opt.eps.end", label="Epsilon End", min = 0, value = 0.04, width="50%"),
                                             br(), br(), br(), br(),
                                             numericInput("db.opt.mpts.end", label = "Min Points End", min = 1, value = 7, step = 1, width="50%")
                                      )
                                  )
                                ),
                                fluidRow(
                                  box(title = "Final DBscan Parameters", width=12,
                                      column(width = 6,
                                             numericInput("db.eps", value = 0, min = 0, label = "Epsilon", width="25%"),
                                             actionButton("db.scan.run", label = "Run", icon = icon("magic", lib="font-awesome"))
                                             ),
                                      column(width = 6,
                                             numericInput("db.mpts", value = 0, min = 0, label = "Min Pts", width="25%")
                                             )
                                      )
                                ),
                                h2("Cluster Refinement"),
                                fluidRow(
                                  column(width=6,
                                         box(plotlyOutput(outputId = "db.clust.plot"), width=12, height="700px")),
                                  column(width=6,
                                         box(title = "tSNE Colour Param - NOT IMPLEMENTED!", width = 12, height = "200px",
                                             selectInput("db.tsne.col", label = "Colour Variable", choices = list(), selected=1, multiple = FALSE)),
                                         box(plotlyOutput(outputId = "db.clust.detail.plot"), width=12, height="500px")
                                  )
                                ),
                                fluidRow(
                                  box(title = "Run Reclustering", width = 12,
                                      uiOutput("reclust.markers"),
                                      actionButton(inputId = "reclust.run", label="Run", icon = icon("magic", lib="font-awesome")))
                                ),
                                fluidRow(
                                  column(width = 6,
                                         box(plotlyOutput(outputId = "reclust.plot"), width=12, height = "650px")),
                                  column(width=6,
                                         box(plotlyOutput(outputId = "reclust.detail.plot"), width=12))
                                )
                                ),
                        tabItem(tabName = "clust_enrich",
                                h2("Identification of Clusters Enriched in Annotations of Interest"),
                                tags$div(
                                  id="enrich_test_params", class="row",
                                  fluidRow(
                                    box(title = "Enrichment Testing Parameters", width = 12,
                                        column(width = 6,
                                               uiOutput(outputId = "enrich.clust"),
                                               checkboxInput(inputId = "enrich.equal",  value = FALSE, width = "100%",
                                                             label = "Use Equal Proportions for Category Enrichment?")),
                                        column(width = 6,
                                               uiOutput("enrich.category"),
                                               actionButton(inputId = "enrich.run", label = "Run", icon=icon("magic", lib="font-awesome")))
                                    )
                                  )
                                )

                        )
                      )
                    )
)


server <- function(input, output, session) {
  # TODO Modularise the shiny app.
  data.folder <- reactiveValues()
  expdata <- reactiveValues()

  ## TODO get the enrichment results, hierarchical clustering data back from the enrichment test module and save along with the rest of expdata.
  ## TODO insert important inputs into expdata prior to save (i.e. any marker choices, enrichment test params, tsne params, clustering params, mem....etc)
  ## TODO implement loading of saved expdata RDS file and update all inputs based on it .... maybe bookmarking will be better here.
  output$fackit.download <- downloadHandler(
    filename = function(){paste('fackit_results', Sys.time(), '.Rds', sep='')},
    content = function(file) {
      saveRDS(reactiveValuesToList(expdata),file = file)
    }
  )

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

  output$table <- renderDT(data.folder[["files"]], editable = TRUE, selection = "none", server = TRUE, options=list(dom="ltip", paging=FALSE))
  proxy = dataTableProxy("table")
  observeEvent(input$table_cell_edit, {
    info = input$table_cell_edit
    i = info$row
    j = info$col
    v = info$value
    data.folder[["files"]][i, j] <<- DT::coerceValue(v, data.folder[["files"]][i, j])
    replaceData(proxy, data.folder[["files"]], resetPaging = FALSE)
  })


  # Check Column Names
  observe({
    path <- input$file1$datapath
    if(is.null(path)){return(NULL)}
    header <- vector(length=length(path), mode="list")
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

  output$column.names <- renderDT(data.folder[["col.names"]], editable = TRUE, selection = "none", server = TRUE, options=list(dom="ltip", paging=FALSE))
  proxy.cols = dataTableProxy("column.names")
  observeEvent(input$column.names_cell_edit, {
    info = input$column.names_cell_edit
    i = info$row
    j = info$col
    v = info$value
    data.folder[["col.names"]][i, j] <<- DT::coerceValue(v, data.folder[["col.names"]][i, j])
    replaceData(proxy.cols, data.folder[["col.names"]], resetPaging = FALSE)
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

    str(raw.data) %>% print

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
  ## TODO make this DT show only 3 sig figs.
  output$cutoffs <- renderDT(expdata[["cutoffs"]], editable = TRUE, selection = "none", server = TRUE, options=list(dom="ltip", paging=FALSE, digits=3))
  proxy.cols = dataTableProxy("cutoffs")
  observeEvent(input$cutoffs_cell_edit, {
    info = input$cutoffs_cell_edit
    i = info$row
    j = info$col
    v = as.numeric(info$value)
    expdata[["cutoffs"]][i, j] <<- DT::coerceValue(v, expdata[["cutoffs"]][i, j])
    replaceData(proxy.cols, expdata[["cutoffs"]], resetPaging = FALSE)
  })

  ## Run Transformation
  observeEvent(input$transform, {
    norm.data <- facsnorm(x=expdata[["raw.data"]][,colnames(expdata[["cutoffs"]])], cutoffs = as.numeric(expdata[["cutoffs"]][1,]), asinCofac = input$asincofac, method = "arcsin")

    norm.data <- cbind(norm.data, expdata[["raw.data"]][,expdata[["metadata"]]])
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
            hoverinfo="skip", marker = list(size = 3, color = 'rgba(0, 0, 0, .5)')) %>%
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
    selectInput(inputId = "mem.groups", label = "Select Conditions", multiple = FALSE, width = "25%", choices = expdata[["metadata"]])
  })

  # Run MEM and plot
  observeEvent(input$run.mem, {
    c("formatting data for MEM") %>% print

    ## BUG this now fails if there is only a single Cond column.
    mem.data <- expdata[["norm.data"]][,c(expdata[["markers.mem"]])] ## TODO once MEM function is rewritten, replace with: expdata[["norm.data"]][,c(expdata[["tsne.markers"]], input$mem.groups)]
    mem.data$cluster <- as.numeric(as.factor(expdata[["norm.data"]][,c(input$mem.groups)])) ## TODO once MEM function is rewritten change this so that the expdata is directly input into the call to MEM without changing cluster to a factor and numeric... etc

    c("Running MEM") %>% print
    if(is.na(input$mem.iqr))
    {
      expdata$mem.res <- MEM(exp_data = mem.data, transform = FALSE, choose.markers = FALSE, choose.ref = FALSE, IQR_thresh = "auto")
    }else{
      expdata$mem.res <- MEM(exp_data = mem.data, transform = FALSE, choose.markers = FALSE, choose.ref = FALSE, IQR_thresh = input$mem.iqr)
    }

    str(mem.data) %>% print

    str(expdata[["mem.res"]][["MAGpop"]][[1]]) %>% print

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
        radioButtons(inputId = "tsne.mode", label = "tSNE Mode", inline = TRUE, c("FFT"=TRUE, "BH"=FALSE)),
        radioButtons(inputId = "tsne.tree", label = "tSNE NN Mode", inline = TRUE, c("Vantage-Point"=FALSE, "ANNOY"=TRUE)),
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

    set.seed(input$seed)

    if(length(find.package("fftRtsne", quiet = T)) != 0){
      "fftRtsne" %>% print
      str(expdata[["norm.data"]]) %>% print
      tsne <- fftRtsne(X = expdata[["norm.data"]][,c(expdata[["tsne.markers"]])],
                       dims = input$tsne.dim, perplexity = input$tsne.perp, check_duplicates = FALSE, max_iter = input$tsne.iter,
                       fft_not_bh = input$tsne.mode, ann_not_vptree = input$tsne.tree, stop_lying_iter = input$tsne.stop.lying.iter,
                       exaggeration_factor = input$tsne.early.exag, no_momentum_during_exag = FALSE, start_late_exag_iter = input$tsne.start.late.exag,
                       late_exag_coeff = input$tsne.late.exag, rand_seed = input$seed)

      expdata$tsne <- data.frame(tsne1=tsne[,1], tsne2=tsne[,2])

      tsne <- fftRtsne(X = expdata[["norm.data"]][,c(expdata[["tsne.markers"]])],
                       dims = 1, perplexity = input$tsne.perp, check_duplicates = FALSE, max_iter = input$tsne.iter,
                       fft_not_bh = input$tsne.mode, ann_not_vptree = input$tsne.tree, stop_lying_iter = input$tsne.stop.lying.iter,
                       exaggeration_factor = input$tsne.early.exag, no_momentum_during_exag = FALSE, start_late_exag_iter = input$tsne.start.late.exag,
                       late_exag_coeff = input$tsne.late.exag, rand_seed = input$seed)

      expdata$tsne1d <- data.frame(tsne1=tsne[,1])

    }else{
      "Rtsne" %>% print
      tsne <- Rtsne(X = expdata[["norm.data"]][,c(expdata[["tsne.markers"]])], dims = input$tsne.dim, perplexity = input$tsne.perp, check_duplicates = FALSE, max_iter = input$tsne.iter,
                    stop_lying_iter = input$tsne.stop.lying.iter, exaggeration_factor = input$tsne.early.exag)

      expdata$tsne <- data.frame(tsne1=tsne$Y[,1], tsne2=tsne$Y[,2])

      tsne <- Rtsne(X = expdata[["norm.data"]][,c(expdata[["tsne.markers"]])], dims = 1, perplexity = input$tsne.perp, check_duplicates = FALSE, max_iter = input$tsne.iter,
                    stop_lying_iter = input$tsne.stop.lying.iter, exaggeration_factor = input$tsne.early.exag)

      expdata$tsne1d <- data.frame(tsne1=tsne$Y[,1])
    }

    "finished tsne" %>% print

    updateSelectInput(session, "tsne.col",
                      choices = as.vector(c(expdata[["tsne.markers"]], expdata[["metadata"]])),
                      selected = as.vector(c(expdata[["tsne.markers"]], expdata[["metadata"]]))[1])

    ## TODO implement or remove... what is this???
  #  updateSelectInput(session, "db.tsne.col",
  #                    choices = as.vector(c(expdata[["tsne.markers"]], expdata[["metadata"]])),
  #                    selected = as.vector(c(expdata[["tsne.markers"]], expdata[["metadata"]]))[1])

    output$reclust.markers <- renderUI({
      checkboxGroupInput(inputId = "reclust.markers",label = "Select Markers", inline = TRUE,
                         choices = expdata[["markers.raw"]], selected = expdata[["markers.raw"]])
    })

  })
  ## TODO Override legend plotting params in plotly (alpha and size are too low for categorical plotting.) Not currently possible??
  ## TODO Make plotting options for the 1D tsne
  observeEvent(input$tsne.col, {
    if(input$tsne.col %in% expdata[["tsne.markers"]]){
      output$tsne.plot <- renderPlotly({
        plot_ly(x = expdata[["tsne"]][,1], y = expdata[["tsne"]][,2], colors = viridis(100), alpha = 0.5, color = expdata[["norm.data"]][,input$tsne.col], type="scattergl", mode = "markers",
                hoverinfo="skip", marker = list(size = 3), width = 600, height = 600) %>%
          layout(xaxis = list(title="tSNE-1"), yaxis = list(title="tSNE-2"), legend=list(markers = list(size=6, alpha=1), font=list(size=12)), scene = list(aspectratio = list(x = 1, y = 1)))
      })
    }else{
      output$tsne.plot <- renderPlotly({
        plot_ly(x = expdata[["tsne"]][,1], y = expdata[["tsne"]][,2], alpha = 0.5,
                color = expdata[["norm.data"]][,input$tsne.col], type="scattergl", mode = "markers", hoverinfo="skip", marker = list(size = 3), width = 600, height = 600) %>%
          layout(xaxis = list(title="tSNE-1"), yaxis = list(title="tSNE-2"), legend=list(markers = list(size=6, alpha=1), font=list(size=12)), scene = list(aspectratio = list(x = 1, y = 1)))
      })
    }
  })

  observeEvent(input$db.knn.run, {
    output$db.knn.plot <- renderPlot({
      kNNdistplot(x = expdata[[input$db.tsne.dim]], k = 4)
    })
  })

  observeEvent(input$db.opt.run, {
    expdata$db.opt <- dbscan.opt(data = expdata[[input$db.tsne.dim]], eps.start = input$db.opt.eps.start, eps.end = input$db.opt.eps.end,
                                 step.size = input$db.opt.eps.step, minPts.start = input$db.opt.mpts.start, minPts.end = input$db.opt.mpts.end)

    output$db.opt.plot <- renderPlot({
      ggplot(expdata[["db.opt"]], aes(x=eps, y=n.clust, colour=noise.pts)) +geom_point() +scale_colour_gradientn(colours=magma(100), trans="log") +facet_wrap(~minPts)
    })
  })

  observeEvent(input$db.scan.run, {
    expdata$dbscan <- dbscan(x=expdata[[input$db.tsne.dim]], eps = input$db.eps, minPts = input$db.mpts)

    # this just to make life a little easier later, not having to worry about factor weirdness.
    expdata[["dbscan"]]$cluster <- as.character(expdata[["dbscan"]]$cluster)
    expdata[["norm.data"]]$db.clust <- as.character(expdata[["dbscan"]]$cluster)

    ## plotly output for exploring the dbscan output.
    if(input$db.tsne.dim == "tsne"){
      output$db.clust.plot <- renderPlotly({
        plot_ly(x=expdata[["tsne"]][which(expdata[["dbscan"]]$cluster != 0),1], y=expdata[["tsne"]][which(expdata[["dbscan"]]$cluster != 0),2],
                color=expdata[["dbscan"]]$cluster[which(expdata[["dbscan"]]$cluster != 0)], key=expdata[["dbscan"]]$cluster[which(expdata[["dbscan"]]$cluster != 0)],
                hoverinfo="none", type = "scattergl", mode = "markers", marker = list(size = 3), width = 600, height = 600, source = "db.clust.plot") %>%
          layout(showlegend=FALSE, xaxis = list(title="tSNE-1"), yaxis = list(title="tSNE-2"), legend=list(markers = list(size=6, alpha=1), font=list(size=12)), scene = list(aspectratio = list(x = 1, y = 1)))
      })
    }else{
      output$db.clust.plot <- renderPlotly({
        ## TODO Make plotting options for 1D tsne
        ## TODO Convert to jitter boxplot without the box.
        plot_ly(x=expdata[["tsne1d"]][which(expdata[["dbscan"]]$cluster != 0),], y=1,
                color=expdata[["dbscan"]]$cluster[which(expdata[["dbscan"]]$cluster != 0)], key=expdata[["dbscan"]]$cluster[which(expdata[["dbscan"]]$cluster != 0)],
                hoverinfo="none", type="scattergl", mode = "markers", marker = list(size = 3), width = 600, height = 600, source="db.clust.plot") %>%
          layout(showlegend=FALSE, xaxis = list(title="tSNE-1"), yaxis = list(title="tSNE-2"), legend=list(markers = list(size=6, alpha=1), font=list(size=12)), scene = list(aspectratio = list(x = 1, y = 1)))
      })
    }
  })

  output$db.clust.detail.plot <- renderPlotly({
    coi <- event_data(event = "plotly_click", source = "db.clust.plot")$key[[1]]
    if(is.null(coi) == TRUE){return(NULL)}

    some.data <- expdata[["norm.data"]][which(expdata[["norm.data"]]$db.clust == coi),]
    some.data <- melt(some.data, measure.vars=expdata[["markers.raw"]])

    some.data$x.jitt <- vipor::offsetX(y = some.data$value, x = some.data$variable)
    some.data$tickval <- NA
    temp.markers <- unique(some.data$variable)
    for(i in 1:length(temp.markers)){
      some.data[which(some.data$variable == temp.markers[i]),"x.jitt"] <- some.data[which(some.data$variable == temp.markers[i]),"x.jitt"] + i
      some.data[which(some.data$variable == temp.markers[i]),"tickval"] <- i
    }

    ## TODO add title to plots
    ## HACK this implementation works, but the tick labels can be a bit messy.
    plot_ly(x=some.data[,"x.jitt"], y=some.data[,"value"], type="scattergl", mode="markers", marker=list(size=4, alpha=0.4), hoverinfo="skip") %>%
      layout(title = paste("Marker Expression in Cluster", coi, sep = " - "), xaxis=list(tickmode="array", ticktext=unique(some.data[,"variable"]), tickvals=unique(some.data[,"tickval"]), tickangle=-90), showlegend=FALSE)
  })


  observeEvent(input$reclust.run, {
    ## TODO figure out what to do with the noise clusters in dbscan - need to recluster amongst the split clusts.
    c("Running clust.split") %>% print
    expdata$clust.split <- clust.split(x = expdata[["norm.data"]], markers = input$reclust.markers, clusters = expdata[["dbscan"]]$cluster)
    c("Running binmat") %>% print
    expdata$bin.list <- binmat(data = expdata[["norm.data"]], cluster.col = "db.clust", markers = input$reclust.markers, split.list = expdata[["clust.split"]], thresh = 0) ## TODO alter bin mat to accept a cluster col that is not part of the data frame
    c("Running split.merge") %>% print

    ## TODO currently has an issue where snlocation cannot find location param problem - occurs mostly with smaller data sets and some system seed values - need way to provide snlocation with predetermined sys.seed
    expdata$split.merge <- splitmerge(x = expdata[["norm.data"]], markers = input$reclust.markers, clust.col = "db.clust", bin.list = expdata[["bin.list"]], noise.clust.id = "0")[,c("split.clusts","super.clusts")]

    c("plotting") %>% print

    if(input$db.tsne.dim == "tsne"){
      output$reclust.plot <- renderPlotly({
        plot_ly(x=expdata[["tsne"]][which(expdata[["split.merge"]]$super.clusts != 0),1], y=expdata[["tsne"]][which(expdata[["split.merge"]]$super.clusts != 0),2],
                color=expdata[["split.merge"]]$super.clusts[which(expdata[["split.merge"]]$super.clusts != 0)], key=expdata[["split.merge"]]$super.clusts[which(expdata[["split.merge"]]$super.clusts != 0)],
                hoverinfo="none", type = "scattergl", mode = "markers", marker = list(size = 3), width = 600, height = 600, source = "reclust.plot") %>%
          layout(showlegend=FALSE, xaxis = list(title="tSNE-1"), yaxis = list(title="tSNE-2"), legend=list(markers = list(size=6, alpha=1), font=list(size=12)), scene = list(aspectratio = list(x = 1, y = 1)))
      })
    }else{
      output$reclust.plot <- renderPlotly({
        ## TODO Make plotting options for 1D tsne
        ## TODO Convert to jitter boxplot without the box.
        plot_ly(hoverinfo="none", type = "scattergl", mode = "markers", source = "reclust.plot") %>%
          add_markers(x=jitter(x=expdata[["tsne1d"]][which(expdata[["split.merge"]]$super.clusts != 0),1]), y=1,
                      color=expdata[["split.merge"]]$super.clusts[which(expdata[["split.merge"]]$super.clusts != 0)],
                      key=expdata[["split.merge"]]$super.clusts[which(expdata[["split.merge"]]$super.clusts != 0)],
                      marker = list(size = 3, alpha=0.4), hoverinfo = "none", showlegend = TRUE) %>%
          layout(xaxis = list(title="tSNE-1"))
      })
    }
  })

  output$reclust.detail.plot <- renderPlotly({
    coi <- event_data(event = "plotly_click", source = "reclust.plot")$key[[1]]
    if(is.null(coi) == TRUE){return(NULL)}

    some.data <- expdata[["norm.data"]][which(expdata[["split.merge"]]$super.clusts == coi),]
    some.data <- melt(some.data, measure.vars=expdata[["markers.raw"]], id.vars=expdata[["metadata"]])

    some.data$x.jitt <- vipor::offsetX(y = some.data$value, x = some.data$variable)
    some.data$tickval <- NA
    temp.markers <- unique(some.data$variable)
    for(i in 1:length(temp.markers)){
      some.data[which(some.data$variable == temp.markers[i]),"x.jitt"] <- some.data[which(some.data$variable == temp.markers[i]),"x.jitt"] + i
      some.data[which(some.data$variable == temp.markers[i]),"tickval"] <- i
    }

    ## HACK this implementation works, but the tick labels can be a bit messy.
    plot_ly(x=some.data[,"x.jitt"], y=some.data[,"value"], type="scattergl", mode="markers", marker=list(size=4, alpha=0.4), hoverinfo="skip") %>%
      layout(title = paste("Marker Expression in Cluster", coi, sep = " - "), xaxis=list(tickmode="array", ticktext=unique(some.data[,"variable"]), tickvals=unique(some.data[,"tickval"]), tickangle=-90), showlegend=FALSE)
  })

  ## TODO fix choice names for DB clust and super.clust/split.merge clusts
  output$enrich.clust <- renderUI({
    if(is.null(expdata[["split.merge"]])){
      selectInput(inputId = "enrich.clust", label = "Select Clustering", multiple = FALSE, width = "25%", choices = c("db.clust"))
    }else{
      selectInput(inputId = "enrich.clust", label = "Select Clustering", multiple = FALSE, width = "25%", choices = c("db.clust","super.clusts"))
    }
  })

  output$enrich.category <- renderUI({
    checkboxGroupInput(inputId = "enrich.category", label = "Select Conditions to Perform Enrichment Testing On:", inline = TRUE,
                       choices = expdata[["metadata"]])
  })


  observeEvent(input$enrich.run, {
    input$enrich.category %>% print
    conds <- input$enrich.category

    lapply(conds,FUN = function(x) {
      insertUI(
        selector = "#enrich_test_params",
        where = "afterEnd",
        ui = tagList(
          h4(paste(x,"Enrichment Results", sep=" ")),
          enrichTest.UI(paste0("enrich.test", x))
        )
      )

      ## TODO make if else controlling cbind of split merge
      ## TODO make if else controlling cbind of the correct tsne object in case 1D tsne is used instead of 2D
      ## BUG if clust.col param is changed, the whole thing breaks... interactive plots will continue to show the first elements selected in the first run of the analysis
      ## also, the first analysis ui's will go blank.
      ## HACK What is the callModule returning in df? How can I access the enrichment results tables and return to the user as download?
      ## TODO return enrichment tables to expdata values (one per cond test), return any hierarchical clustering data.
      df <- callModule(enrichTest.module,
                       paste0("enrich.test", x),
                       data=reactive(cbind(expdata[["norm.data"]], expdata[["tsne"]], expdata[["split.merge"]])),
                       cats = x,
                       clust.col = input$enrich.clust,
                       ## tsne.dim = input$tsne.dim,
                       equal.props = input$enrich.equal,
                       markers = expdata[["markers.raw"]]
      )

      expdata[[paste0("enrich.test.", x)]] <- df()

    })

  })


}

# Run the application
shinyApp(ui = ui, server = server)
