nodeSize <- function(somtree, reset=FALSE, maxNodeSize=0.1)
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
    somtree.res$map$nodes$size <- as.vector(rep(maxNodeSize, nrow(somtree.res$map$nodes$codes)))
  }
  else
  {
    t <- table(somtree.res$map$mapped$bmn)
    t <- sqrt(t)
    scale <- max(t)
    rescaled <- maxNodeSize * t/scale
    somtree.res$map$nodes$size <- numeric(nrow(somtree.res$map$nodes$codes))
    somtree.res$map$nodes$size[as.numeric(names(t))] <- as.vector(rescaled)
  }
  somtree.res
}


