library(shiny)
library(shinydashboard)
library(shinyjs)
library(plotly)
library(DT)
library(vipor)
library(stringr)
library(parallel)

library(ggplot2)
library(pals)
library(ggsci)
library(ggthemes)
library(pheatmap)
library(FACkit)
library(scales)

library(amap)

enrichTest.UI <- function(id) {
  ns <- NS(id)

  # TODO make sure that the second row of boxes are on their own row
  fluidPage(
    fluidRow(
      column(width = 6,
             box(width = 12, height="650px", ##
                 title = paste("Clusters Enriched in",str_split(id, pattern = "enrich.test", simplify = T)[,2], sep=" "),
                 DTOutput(outputId = ns("enrich.tab"))
                 )
      ),
      column(width=6,
             box(width = 12, height="650px", ##
                 title = paste("Proportion of",str_split(id, pattern = "enrich.test", simplify = T)[,2],"in Selected Clusters",sep=" "),
                 plotlyOutput(outputId = ns("enrich.prop"))
             )
      )
    ),
    fluidRow(
      column(width = 6,
             box(width = 12, height = "700px",
                 title = paste("tSNE Location of Clusters Enriched in",str_split(id, pattern = "enrich.test", simplify = T)[,2],sep=" "),
                 plotlyOutput(outputId = ns("enrich.tsne"))
             )
             ),
      column(width = 6,
             box(width = 12, height = "700px",
                 plotlyOutput(outputId = ns("enrich.tsne.detail"))
                 )
             )
      ),
    fluidRow(
      box(title = "Heatmap Parameters", width = 12,
          uiOutput(outputId = ns("enrich.hm.markers")),
          actionButton(inputId = ns("enrich.hm.run"), label="Run Heatmap", icon=icon("magic", lib="font-awesome"))),
      box(title = "Heatmap of Cells in Selected Clusters", width = 12,
          plotOutput(outputId = ns("enrich.hm"))
      )
      )

  )



}

enrichTest.module <- function(input, output, session, data, cats, clust.col, tsne.dim, equal.props, markers) {
  # Yields the data frame with an additional column "selected_"
  # that indicates whether that observation is brushed
  ns <- NS(paste0("enrich.test",cats))


  ## TODO Figure out method to return hclust data from the heatmapping module... only enrichment.results table is returned
  ## hm hclust must be put into a reactive({}) somehow but only called when input$enrich.hm.run is changed (maybe EventReactive could work)
  ## else, the hm could be placed into a separate module??? how to get insertUI modules to talk to each other???
  enrichment.results <- reactive({
    ## TODO make noise clust.id dynamic somehow  - only needed when you add another clustering method.
    enrichTest(x= data(), clust.col = clust.col, noise.clust.id = "0", cat.col = cats, equal.props = equal.props, alternative = "greater")
  })

  output$enrich.hm.markers <- renderUI({
    checkboxGroupInput(inputId = ns("enrich.hm.markers"), label="Select Markers to Use in Heatmap", choices = markers, selected=markers, inline=TRUE)
  })

  ## TODO customise the filtering and layout options - dropbox for category column, numeric filtering for the pval, prop and FDR columns.
  ## TODO change the output columns to use the rounded values outputted by the enrichTest function [,c("cluster","category","cluster.size","prop","bin.pval","bin.FDR)]
  output$enrich.tab <- renderDT(enrichment.results()[,c("cluster","category","cluster.size","proportion","binomial.pval","binomial.FDR")], server = TRUE, options=list(digits=3, dom="ltip"))

  output$enrich.prop <- renderPlotly({
    if(is.null(input$enrich.tab_rows_selected)){return(NULL)}

    c("plotting") %>% print
    cois <- as.character(enrichment.results()[input$enrich.tab_rows_selected,"cluster"])

    prop.mat <- as.data.frame(table(
      data()[which(data()[,clust.col] %in% enrichment.results()[input$enrich.tab_rows_selected,"cluster"]),c(clust.col,cats)]
    ), stringsAsFactors=FALSE)

    prop.mat[,clust.col] <- as.character(prop.mat[,clust.col])
    prop.mat <- prop.mat[which(prop.mat[,clust.col] %in% cois),]

    for(i in 1:length(cois)){
      prop.mat[which(prop.mat[,clust.col] == cois[i]),"Freq"] <- signif(c((prop.mat[which(prop.mat[,clust.col] == cois[i]),"Freq"]/sum(prop.mat[which(prop.mat[,clust.col] == cois[i]),"Freq"])) * 100), digits = 4)
    }

    plot_ly(x=prop.mat[,clust.col], y=prop.mat[,"Freq"], color=prop.mat[,cats], type="bar", hoverinfo="skip") %>% layout(barmode="stack")


  })

  output$enrich.tsne <- renderPlotly({

    ## TODO Need to get input indicating if using 2D or 1D tsne - CHANGE THE SCRIPT!!
    ##data()[which(data()[,clust.col] %in% enrichment.results()[input$enrich.tab_rows_selected,"cluster"]),]

    if(is.null(input$enrich.tab_rows_selected)){
      plot_ly(x=data()[,"tsne1"], y=data()[,"tsne2"], color=data()[,cats], colors=c(colorblind_pal()(length(unique(data()[,cats]))),"#CFCFCF"),
              key=data()[,clust.col], source=paste(cats,"enrich.tsne.plot",sep="."), hoverinfo="none", type="scattergl", mode="markers",
              marker=list(size=3, alpha=0.4), width = 600, height = 600) %>%
        layout(showlegend=TRUE, xaxis = list(title="tSNE-1"), yaxis = list(title="tSNE-2"), legend=list(markers = list(size=6, alpha=1), font=list(size=12)),
               scene = list(aspectratio = list(x = 1, y = 1)))
    }else{
      # TODO Find better method to colour unselected clusters grey to highlight the selected clusters.
      data.sub <- data()[,c(clust.col, cats, "tsne1","tsne2")]

      data.sub[,cats] <- as.character(data.sub[,cats])

      data.sub[!(data()[,clust.col] %in% enrichment.results()[input$enrich.tab_rows_selected,"cluster"]),cats] <- "not.selected"

      ## Add if/else statement using input$tsne.dim to control 1D/2D tsne scatter plot drawing.

      ## TODO Make return coloured plot without the greying of unselected clusters.
      plot_ly(x=data.sub[,"tsne1"], y=data.sub[,"tsne2"], color=data.sub[,cats], colors=c(colorblind_pal()(length(unique(data()[,cats]))),"#CFCFCF"),
              key=data.sub[,clust.col], source=paste(cats,"enrich.tsne.plot",sep="."), hoverinfo="none", type="scattergl", mode="markers",
              marker=list(size=3, alpha=0.4), width = 600, height = 600) %>%
        layout(showlegend=TRUE, xaxis = list(title="tSNE-1"), yaxis = list(title="tSNE-2"), legend=list(markers = list(size=6, alpha=1), font=list(size=12)),
               scene = list(aspectratio = list(x = 1, y = 1)))
    }
  })

  output$enrich.tsne.detail <- renderPlotly({
    coi <- event_data(event = "plotly_click", source = paste(cats,"enrich.tsne.plot",sep="."))$key[[1]]
    if(is.null(coi) == TRUE){return(NULL)}

    some.data <- data()[which(data()[,clust.col] == coi),]
    some.data <- melt(some.data, measure.vars=markers)

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

  ## TODO implement heatmapping with check on nrows - if > 45,000 do not run and open dialog box warning. Use action button to start the heatmapping once clusters are selected.
  ## open "are you sure" dialog box if heatmapping > 20,000 rows warning of time required - allow user to choose to stop the heatmap
  ##

  observeEvent(input$enrich.hm.run, {

    c("Running Heatmap") %>% print

    if(is.null(input$enrich.tab_rows_selected)){return(NULL)}

    input$enrich.hm.markers %>% print

    ## TODO this is a bit laggy, need to find out why this runs slow and fix it.
    data.sub <- data()[which(data()[,clust.col] %in% enrichment.results()[input$enrich.tab_rows_selected,"cluster"]),c(input$enrich.hm.markers, cats)]

    str(data.sub) %>% print


    nrows.selected <- nrow(data.sub)
    if(nrows.selected > 45000){
      ## TODO add in a downsampling method for clustering giant groups of cells - use equal numbers of each cluster or use a probability weight based on proportions of each cluster/cond
      showModal(
        modalDialog(title = "ERROR!",
                    "You are trying to heatmap more than 45,000 points - This is impossible!\n Select smaller clusters",
                    size = "s")
        )
      return(NULL)
    }
    if(nrows.selected < 45000 & nrows.selected > 20000){
      ## TODO alter the modal to have 2 buttons, "Yes" and "Cancel".
      showModal(
        modalDialog(title = "WARNING!",
                    "Heatmapping this many cells could take a while, are you sure you want to run this?",
                    size="s")
      )

      ## TODO add hcluster par params to the UI.
      cell.hclust <- hclusterpar(x = data.sub[,input$enrich.hm.markers], method = "euclidean", nbproc = parallel::detectCores() - 1, link = "complete")
      marker.hclust <- hclusterpar(x = t(data.sub[,input$enrich.hm.markers]), method = "euclidean", nbproc = parallel::detectCores() - 1, link = "complete")

      ## TODO add method to add other annotation layers to the heatmap... like other conds, etc.
      ## TODO add method to determine colour palette to use base on the number of unique elements in cats.
      cond <- colorblind_pal()(length(unique(data()[,cats])))
      names(cond) <- unique(data()[,cats])

      anno_layer <- list(cats=cond)

      output$enrich.hm <- renderPlot({
        pheatmap(mat=t(data.sub[,input$enrich.hm.markers]), cluster_rows = marker.hclust, cluster_cols = cell.hclust, treeheight_row = 25, treeheight_col = 30,
                 legend = TRUE, show_colnames = FALSE, cellheight = 20, annotation_colors = anno_layer,
                 annotation_col = as.data.frame(matrix(data=data.sub[,cats], nrow = nrow(data.sub),
                                                       ncol = length(anno_layer), byrow = T,
                                                       dimnames = list(c(rownames(data.sub)), c(names(anno_layer))))))
      })
    }else{
      ## TODO add hclusterpar params to the UI.
      ## TODO add parallel to package imports
      cell.hclust <- hclusterpar(x = data.sub[,input$enrich.hm.markers], method = "euclidean", nbproc = detectCores() - 1, link = "complete")
      marker.hclust <- hclusterpar(x = t(data.sub[,input$enrich.hm.markers]), method = "euclidean", nbproc = detectCores() - 1, link = "complete")

      ## TODO add method to add other annotation layers to the heatmap... like other conds, etc.
      ## TODO add method to determine colour palette to use base on the number of unique elements in cats.
      cond <- colorblind_pal()(length(unique(data()[,cats])))
      names(cond) <- unique(data()[,cats])

      anno_layer <- list(cats=cond)

      output$enrich.hm <- renderPlot({
        pheatmap(mat=t(data.sub[,input$enrich.hm.markers]), cluster_rows = marker.hclust, cluster_cols = cell.hclust, treeheight_row = 25, treeheight_col = 30,
                 legend = TRUE, show_colnames = FALSE, cellheight = 20, annotation_colors = anno_layer,
                 annotation_col = as.data.frame(matrix(data=data.sub[,cats], nrow = nrow(data.sub),
                                                       ncol = length(anno_layer), byrow = T,
                                                       dimnames = list(c(rownames(data.sub)), c(names(anno_layer))))))
      })
    }


  })

  return(enrichment.results)
}

