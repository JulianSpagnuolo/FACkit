clust.split <- function(x, markers, clusters)
{
  #' @author Julian Spagnuolo
  #' @param x data.frame or matrix of numeric expression data
  #' @param markers character vector of column names in x to use in the cluster splitting
  #' @param clusters character vector of cluster ids. length must equal nrow of x.
  #' 
  #' @importFrom diptest dip.test
  
  
  clust.ids <- unique(clusters)
  
  to.split <- vector(mode="list", length=length(clust.ids))
  names(to.split) <- clust.ids
  
  for(i in 1:length(unique(clust.ids)))
  {
    markers.to.split <- vector()
    for(n in 1:length(markers))
    {
      p <- dip.test(x=x[which(clusters == clust.ids[i]),markers[n]], simulate.p.value = TRUE, B=100)$p.value
      if(p < 0.05) # if expression of marker n is bimodal in the cluster then cluster needs to be split using adapt clust method
      {
        bimod.true <- sign(quantile(x[which(clusters == clust.ids[i]),markers[n]], c(0.25, 0.95)))
        if(isTRUE(all(bimod.true == -1) | all(bimod.true == 1)))
        {
          bimod.true <- FALSE
        }
        else{ bimod.true <- TRUE}
        if(isTRUE(bimod.true))
        {
          markers.to.split <- append(markers.to.split, markers[n])
        }
      }
      if(length(markers.to.split != 0))
      {
        to.split[[clust.ids[i]]] <- markers.to.split
      }
      else{
        to.split[[clust.ids[i]]] <- NA
      }
    }
  }
  to.split <- to.split[!is.na(to.split)]
  return(to.split)
}