adapt.clust <- function(data, nbins=4, bin.defs=NULL, markers, dist="eucl", growfact=3, dist.thresh=4, maxit=100, minpts = 5)
{
  #' @author Julian Spagnuolo
  #' @title Clustering by Binning
  #' @description Cluster large datasets that cannot ordinarily be clustered by distance calculation because of size limitations. It works by creating bins into which observations are placed by nearest distance, currently euclidean distance is implemented.
  #' @param data data.frame or matrix containing the data to be clustered
  #' @param nbins number of bins created at the start. If set to NULL and bin.mat is given, then starting nodes will be initialised using a perumuation matrix based on the positive and negative expression values given in bin.mat. Default is 4.
  #' @param bin.defs 2*length(markers) numeric matrix giving the expression values for positive and negative phenotypes for each marker. If set to NULL, the number of initial bins will be determined by nbins using the density distribution of the data. Default is NULL.
  #' @param markers vector of marker names in data to be clustered
  #' @param dist character vector indicating which distance metric to be used; one of "eucl" = euclidean, "jacc" = Jaccard.
  #' @param growfact numeric indicates the number of new nodes created if average dist is above threshold
  #' @param dist.thresh numeric minimum average distance threshold required to prevent a node splitting and forming new nodes.
  #' @param maxit integer maximum allowed interations of the algorithm
  #' @param minpts minimum node membership required for resampling the node def from the old node members.
  #'
  #' @import gtools

  ## First create the bins
  cat("Initialising\n")

  if(is.null(bin.defs))
  {
    bin.mat <- matrix(data = NA, nrow = nbins, ncol = length(markers), dimnames = list(c(1:nbins),markers))
    for(i in markers)
    {
      x <- density(data[,i], n = round(sqrt(nrow(data)), digits = 0))
      bin.mat[,i] <- sample(x = x$x, size = nbins, prob = x$y)
    }
  }
  else
  {
    cat("\nSetting up initial nodes")
    bin.mat <- permutations(n=2, r=length(markers), v=c(0,1), repeats.allowed=T)
    for(i in 1:length(markers))
    {
      bin.mat[which(bin.mat[,i] == 0),i] <- bin.defs["0",i]
      bin.mat[which(bin.mat[,i] == 1),i] <- bin.defs["1",i]
    }
    nbins <- nrow(bin.mat)
  }

  cat("\n",nrow(bin.mat),"Initial nodes created")
  cat("\nPopulating nodes")
  node.dists <- bindist(binmat=bin.mat, data=as.matrix(data[,markers]))

  cat("\nCulling initial nodes")
  keep <- vector()
  for(i in 1:length(node.dists$node.pops))
  {
    if(!is.null(node.dists$node.pops[[i]]))
    {
      keep <- append(x = keep, values = i)
    }
  }
  bin.mat <- bin.mat[keep,]
  cat("\n",nrow(bin.mat),"Initial nodes populated")

  nodes <- vector(mode="list", length=nrow(bin.mat))
  for(i in 1:length(keep))
  {
    nodes[[i]]$node <- bin.mat[i,]
    nodes[[i]]$datapoint <- node.dists$node.pops[[keep[i]]]
    nodes[[i]]$dist <- node.dists$node.dists[[keep[i]]]
  }
  nbins <- nrow(bin.mat)
  rm(node.dists)

  cat("\nEntering Grow Phase")
  iter <- 0
  while(iter < maxit)
  {
    cat("\nGrowing")
    # reset the number of bins and the bin matrix

    nbins <- 0
    bin.mat <- matrix(data=NA, nrow=0, ncol=length(markers), dimnames = list(c(),markers))

    # Figure out how many new nodes to grow
    for(i in 1:length(nodes))
    {
      if(length(nodes[[i]]$dist) != 0)
      {
        if(mean(nodes[[i]]$dist) > dist.thresh)
        {
          nbins <- nbins + growfact
          temp <- matrix(data=NA, nrow=growfact, ncol=length(markers), dimnames=list(c(),markers))
          for(j in markers)
          {
            if(length(nodes[[i]]$datapoint) > minpts)
            {
              x <- density(x = data[nodes[[i]]$datapoint,j]) # if node membership is high enough resample node def from members
              temp[,j] <- sample(x$x, size = growfact, prob = x$y)
            }
            else{ temp[,j] <- sample(data[nodes[[i]]$datapoint,j], size = growfact, replace = T) }
          }
          bin.mat <- rbind(bin.mat, temp) # add the new nodes to the bin matrix
        }

        else # if node distance is safe keep it.
        {
          nbins <- nbins + 1
          temp <- matrix(data=NA, nrow=1, ncol=length(markers), dimnames=list(c(),markers))
          for(j in markers)
          {
            if(length(nodes[[i]]$datapoint) > minpts)
            {
              x <- density(x = data[nodes[[i]]$datapoint,j]) # if node membership is high enough resample node def from members
              temp[,j] <- sample(x$x, size = 1, prob = x$y)
            }
            else{ temp[,j] <- nodes[[i]]$node[j] } #or keep it if node membership is too low.
          }
          bin.mat <- rbind(bin.mat, temp) # add the node to the bin matrix
        }
      }
    }
    rm(nodes)

    # Reset names of bin.mat and remove old nodes.
    dimnames(bin.mat) <- list(1:nrow(bin.mat), markers)

    node.dists <- bindist(binmat=bin.mat, data=as.matrix(data[,markers]))

    keep <- vector()
    for(i in 1:length(node.dists$node.pops))
    {
      if(!is.null(node.dists$node.pops[[i]]))
      {
        keep <- append(x = keep, values = i)
      }
    }
    bin.mat <- bin.mat[keep,]

    nodes <- vector(mode="list", length=nrow(bin.mat))
    for(i in 1:length(keep))
    {
      nodes[[i]]$node <- bin.mat[i,]
      nodes[[i]]$datapoint <- node.dists$node.pops[[keep[i]]]
      nodes[[i]]$dist <- node.dists$node.dists[[keep[i]]]
    }
    rm(node.dists)
    cat("Size of SOM is now ",length(nodes),"\n")

    iter <- iter + 1
  }



  return(nodes)

}
