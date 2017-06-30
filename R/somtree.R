som.tree <- function(data, markers, shape=c("rect","hex"), maxit=500, cores=1)
{
  #' @author Julian Spagnuolo
  #' @title Build Growing Self-organising Map and Construct a Tree
  #' @param data data.frame or matrix containing expression data
  #' @param markers vector of markers/columns in data to build the SOM with
  #' @param shape vector determining the shape of SOM to build, either rectangular, "rect", or hexagonal, "hex".
  #' @param maxit number of iterations the SOM building function will use, default = 500.
  #' @param cores number of cores to using in distance computation. Default = 1, set higher if you have more than 1 cpu to use.
  #'
  #'
  #' @importFrom amap Dist

  results <- vector(mode="list")
  cat("Training SOM\n")
  t1 <- Sys.time()
  som <- train.gsom(data=data[,markers], iterations = maxit, nhood = shape, keepdata = FALSE)
  results$som <- som
  cat("Mapping data to SOM\n")
  results$map <- map.gsom(gsom_object = som, df=data[,markers], retaindata = FALSE)

  cat("Building tree\n")
  t2 <- Sys.time()
  adj.mat <- Dist(results$map$nodes$codes, method="euclidean", nbproc=cores)
  results$full.graph <- graph.adjacency(as.matrix(adj.mat), mode="undirected", weighted = TRUE)
  tx <- Sys.time() - t2
  cat("Tree built in ", round(as.numeric(tx), digits = 2),"seconds\n")

  results$mst <- minimum.spanning.tree(results$full.graph)

  tx <- Sys.time() - t1
  cat("completed in", round(as.numeric(tx), digits=2), "seconds\n")
  return(results)


}
