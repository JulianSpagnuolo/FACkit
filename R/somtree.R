som.tree <- function(data, markers, shape=c("rect","hex"), maxit=500, cores=1, alpha=0.9, beta=0.5, spread=0.9, alg=c("kohonen"), dim=c(NULL,NULL))
{
  #' @author Julian Spagnuolo
  #' @title Build Growing Self-organising Map and Construct a Tree
  #' @param data data.frame or matrix containing expression data
  #' @param markers vector of markers/columns in data to build the SOM with
  #' @param shape vector determining the shape of SOM to build, either rectangular, "rect", or hexagonal, "hex".
  #' @param maxit number of iterations the SOM building function will use, default = 500.
  #' @param cores number of cores to using in distance computation. Default = 1, set higher if you have more than 1 cpu to use.
  #' @param alpha discount factor for the learning rate during the growing phase of the training. Values should be between 0 and 1. Default = 0.9
  #' @param beta propagation rate. Determines the rate at which the error of a node, that cannot grow any nodes, is passed on to its neighbours. Suggested values may range between 0 and 1. Default = 0.5
  #' @param spread numeric value between 0 and 1, controls the rate at which new nodes are created in the growing SOM algorithm, lower values decrease the rate, higher ones increase it. Default is 0.9
  #' @param alg character vector of length 1. Determines which SOM algorithm to use, choices include either Kohonen SOM, "kohonen" or Growing SOM, "grow" - GrowingSOM is not implemented - awaiting a bug fix and rerelease.
  #' @param dim integer vector of length 2. Determines dimensions of the SOM in the x and y dimension, respectively. These parameters must be set for the kohonen algorithm. If using the GrowSOM algorithm, defining dim will force the growing SOM to use a gridsize of dimension dim[1]*dim[2]. Default is c(NULL,NULL).
  #'
  #'
  #'
  #' @importFrom amap Dist
  results <- vector(mode="list")
  #if(alg == "grow")
  #{
  #  results$algorithm <- "GrowSOM"
  #  cat("Training Growing SOM\n")
  #
  #  if(is.null(dim) == TRUE)
  #  {
  #    gridsize <- FALSE
  #  }
  #  else
  #  {
  #    gridsize <- dim[1]*dim[2]
  #  }
  #  t1 <- Sys.time()
  #  results$som <- train.gsom(data=data[,markers], iterations = maxit, nhood = shape, keepdata = FALSE, alpha = alpha, beta=beta, spreadFactor = spread)
  #  cat("Mapping data to SOM\n")
  #  results$map <- map.gsom(gsom_object = som, df=data[,markers], retaindata = FALSE)
  #
  #  cat("Building tree\n")
  #  t2 <- Sys.time()
  #  adj.mat <- Dist(results$map$nodes$codes, method="euclidean", nbproc=cores)
  #  results$full.graph <- graph.adjacency(as.matrix(adj.mat), mode="undirected", weighted = TRUE)
  #  tx <- Sys.time() - t2
  #  cat("Tree built in ", round(as.numeric(tx), digits = 2),"seconds\n")
  #
  #  results$mst <- minimum.spanning.tree(results$full.graph)
  #
  #  tx <- Sys.time() - t1
  #  cat("completed in", round(as.numeric(tx), digits=2), "seconds\n")
  #}
  if(alg == "kohonen")
  {
    if(is.null(dim))
    {
      cat("You MUST define dim parameter for kohonen algorithm!!\n")
      stop()
    }
    results$algorithm <- "kohonen"
    t1 <- Sys.time()
    cat("Training Kohonen SOM\n")
    if(cores > 1)
    {
      mode = "pbatch"
    }
    else
    {
      mode = "online"
    }
    if(shape == "hex")
    {
      shape = "hexagonal"
    }
    if(shape == "rect")
    {
      shape = "rectangular"
    }
    results$som <- kohonen::som(X=data[,markers], rlen=maxit, alpha=c(0.05, 0.01), mode=mode, cores=cores, dist.fcts="euclidean", grid=somgrid(xdim=dim[1], ydim=dim[2], topo=shape), keep.data=F)
    cat("Mapping data to SOM\n")
    results$map <- kohonen::map(x=results$som, newdata=data[,markers])

    cat("Building tree\n")
    t2 <- Sys.time()
    adj.mat <- Dist(results$som$codes[[1]], method="euclidean", nbproc=cores)

    results$full.graph <- graph.adjacency(as.matrix(adj.mat), mode="undirected", weighted = TRUE)
    tx <- Sys.time() - t2
    cat("Tree built in ", round(as.numeric(tx), digits = 2),"seconds\n")

    results$mst <- minimum.spanning.tree(results$full.graph)

    tx <- Sys.time() - t1
    cat("completed in", round(as.numeric(tx), digits=2), "seconds\n")
  }

  return(results)


}
