adapt.clust <- function(data, nbins=4, markers, dist="eucl", growfact=3, dist.thresh=4, maxit=100, minpts = 5)
{
  #' @author Julian Spagnuolo
  #' @title Clustering by Binning
  #' @description Cluster large datasets that cannot ordinarily be clustered by distance calculation because of size limitations. It works by creating bins into which observations are placed by nearest distance, currently euclidean distance is implemented.
  #' @param data data.frame or matrix containing the data to be clustered
  #' @param nbins number of bins created at the start
  #' @param markers vector of marker names in data to be clustered
  #' @param dist character vector indicating which distance metric to be used; one of "eucl" = euclidean, "jacc" = Jaccard.
  #' @param growfact numeric indicates the number of new nodes created if average dist is above threshold
  #' @param dist.thresh numeric minimum average distance threshold required to prevent a node splitting and forming new nodes.
  #' @param maxit integer maximum allowed interations of the algorithm
  #' @param minpts minimum node membership required for resampling the node def from the old node members.
  #'
  #'

  ## First create the bins
  cat("Initialising\n")
  bin.mat <- matrix(data = NA, nrow = nbins, ncol = length(markers), dimnames = list(c(1:nbins),markers))
  for(i in markers)
  {
    x <- density(data[,i], n = round(sqrt(nrow(data)), digits = 0))
    bin.mat[,i] <- sample(x = x$x, size = nbins, prob = x$y)
  }

  ## Match the data points to the best bin

  nodes <- vector(mode="list", length=nbins)
  for(i in 1:length(nodes))
  {
    nodes[[i]]$node <- bin.mat[i,]
    nodes[[i]]$datapoint <- vector("integer")
    nodes[[i]]$dist <- vector("numeric")
  }


  for(i in 1:nrow(data))
  {
    dists <- apply(X=bin.mat, MARGIN = 1, FUN= function(x,y) {sqrt(sum((x-y)^2))}, y=data[i,markers])
    nodes[[as.numeric(names(dists[which(dists == min(dists))]))]]$datapoint <- append(x = nodes[[as.numeric(names(dists[which(dists == min(dists))]))]]$datapoint, values = i)
    nodes[[as.numeric(names(dists[which(dists == min(dists))]))]]$dist <- append(x = nodes[[as.numeric(names(dists[which(dists == min(dists))]))]]$dist, values = as.numeric(dists[which(dists == min(dists))]))
  }


  iter <- 0
  while(iter < maxit)
  {
    cat("Growing\n")
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

    # Reset names of bin.mat and remove old nodes.
    dimnames(bin.mat) <- list(1:nrow(bin.mat), markers)
    rm(nodes)
    nodes <- vector(mode="list", length=nbins)
    for(i in 1:length(nodes))
    {
      nodes[[i]]$node <- bin.mat[i,]
      nodes[[i]]$datapoint <- vector("integer")
      nodes[[i]]$dist <- vector("numeric")
    }
    cat("Size of SOM is now ",length(nodes),"\n")
    for(i in 1:nrow(data))
    {
      dists <- apply(X=bin.mat, MARGIN = 1, FUN= function(x,y) {sqrt(sum((x-y)^2))}, y=data[i,markers])
      nodes[[as.numeric(names(dists[which(dists == min(dists))]))]]$datapoint <- append(x = nodes[[as.numeric(names(dists[which(dists == min(dists))]))]]$datapoint, values = i)
      nodes[[as.numeric(names(dists[which(dists == min(dists))]))]]$dist <- append(x = nodes[[as.numeric(names(dists[which(dists == min(dists))]))]]$dist, values = as.numeric(dists[which(dists == min(dists))]))
    }



    iter <- iter + 1
  }



  return(nodes)

}
