nodeSize <- function(somtree, reset=FALSE, maxNodeSize=15)
  #'
  #' @title Scaled Node Sizes for growing SOM
  #' @author Julian Spagnuolo
  #' @description
  #' @param somtree results from calling som.tree
  #' @param reset Logical. If TRUE, node size is reset to maxNodeSize, default = FALSE
  #' @param maxNodeSize Maximum size for any node. Default = 15
  #'
  #'
  #'
  #'
  #'
  #'
  #'
  #'
{
  if(reset == TRUE)
  {
    somtree.res$mst$size <- rep(maxNodeSize, nrow(somtree.res$map$mapped))
  }
  else
  {
    t <- table(somtree.res$map$mapped$bmn)
    t <- (t/sum(t)) * nrow(somtree.res$map$mapped)
    shift <- min(t)
    scale <- max(t - shift)
    rescaled <- maxNodeSize * (t - shift)/scale
    somtree.res$map$nodes$size <- numeric(nrow(somtree.res$map$mapped))
    somtree.res$map$nodes$size[as.numeric(names(t))] <- rescaled
  }
  somtree.res
}



