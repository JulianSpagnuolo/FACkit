som.tree <- function(data, markers, shape=c("rect","hex"), maxit=500, cores=1, graph.type="mst")
{
  #' @author Julian Spagnuolo
  #' @title Build Growing Self-organising Map and Construct a Tree
  #' @param data data.frame or matrix containing expression data
  #' @param markers vector of markers/columns in data to build the SOM with
  #' @param shape vector determining the shape of SOM to build, either rectangular, "rect", or hexagonal, "hex".
  #' @param maxit number of iterations the SOM building function will use, default = 500.
  #' @param cores number of cores to using in distance computation. Default = 1, set higher if you have more than 1 cpu to use.
  #' @param graph.type Type of graph to draw, either minimum spanning tree ("mst") or consensus hierarchical tree ("hrg"). Default is "mst".
  #'
  #'
  #' @importFrom amap Dist

  cat("Training SOM\n")
  som <- train.gsom(data=data[,markers], iterations = maxit, nhood = shape, keepdata = FALSE)
  cat("Mapping data to SOM\n")
  mapped <- map.gsom(gsom_object = som, df=data[,markers], retaindata = FALSE)

  cat("Building tree\n")
  adj.mat <- Dist(mapped$nodes$codes, method="euclidean", nbproc=cores)
  full.graph <- graph.adjacency(adj.mat, mode="undirected")
  if(graph.type == "hrg")
  {
    graph <- consensus_tree(graph=full.graph)
    graph <- layout_with_fr(graph=graph, niter = maxit)
  }
  if(graph.type == "mst")
  {
    graph <- minimum.spanning.tree(full.graph)
    graph <- layout_with_fr(graph=graph, niter=maxit)
  }
  return(graph)


}
