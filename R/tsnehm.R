#' Compute t-SNE Heatmap
#'
#' Experimental version of t-SNE heatmaps codes
#'
#' @param expression_matrix FACS data. Columns should be markers, rows are cells.
#' @param markers markers of interest
#' @param tsne_embedding 1D t-SNE embedding
#' @param cell_labels Labels for each cell. Typically these are the cluster assignments for each cell, as obtained by 
#'                    dbscan on the 2D t-SNE. This way, the columns of the heatmap will have a color assigned to each of them,
#'                    and they can be mapped to a corresponding location on the 2D t-SNE
#' @param breaks Number of bins
#' @param plot Logical. Whether to plot the heatmap (TRUE). If false, a list containing a matrix of the heatmap data and the column annotations (group_label) will be returned. Default is TRUE
#' @import pheatmap

tsnehm <- function(data, markers, tsne_embedding, cell_labels, breaks=100, plot=TRUE){
  
  goidf <- data.frame(x=tsne_embedding, data[,markers])
  group <- base::split(x=goidf, f=cut(goidf$x, breaks=breaks))
  
  
  bin.meds <- matrix(nrow=length(markers), ncol = breaks, dimnames = list(c(markers), c(1:breaks)))
  for(i in 1:length(group))
  {
    bin.meds[,i] <- apply(group[[i]][,-1], MARGIN = 2, FUN = function(x) {median(x)})
  }
  
  
  #assign a label to each column based on which of the cell_labels is the most common ### JS Not necessary for FACS.
  dbscan_tsne <- data.frame(x=tsne_embedding, y=cell_labels)
  dbscan_group <- split (dbscan_tsne, cut(dbscan_tsne$x,breaks = breaks))
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  group_labels <-unlist(lapply(dbscan_group, function(x) Mode(x[,2])))
  colnames(bin.meds) <- sprintf("bin:%s, label: %d", 1:breaks, group_labels)
  group_labels[is.na(group_labels)] <- NA
  names(group_labels) <- sprintf("bin:%s, label: %d", 1:breaks, group_labels)
  if(isTRUE(plot))
  {
    pheatmap(mat = bin.meds, cluster_rows = T, cluster_cols = F,
             annotation_col = as.data.frame(group_labels), annotation_names_col = F, show_colnames = F)
  }
  else{
    results <- list()
    results$bin.meds <- bin.meds
    results$group_labels <- group_labels
    return(results)
  }
}