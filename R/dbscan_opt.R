dbscan.opt <- function(data, eps.start, eps.end, step.size, minPts.start=NULL, minPts.end=NULL)
{
  #' @param data continuous numeric data to pass to dbscan
  #' @param eps.start numeric, epsilon to begin scan from
  #' @param eps.end numeric, epsilon to scan to
  #' @param step.size numeric, rate at which to increase the epsilon parameter
  #' @param minPts.start integer, minimum k nearest neighbours for a cluster, if NULL defaults to ncol(data)+1.
  #' @param minPts.end integer, minimum k nearest neighbours to scan to. If defined, function will scan all epsilons between eps.start and eps.end for all minPts between minPts.start and minPts.end. Default is NULL.
  #' @return data.frame. eps = vector of epsilons scanned through; minPts = vector of minPts scanned through; n.clust = number of clusters found in each epsilon step for minPts; noise.pts = number of points considered to be noise at each step; log.noise.pts = log10 of the number of noise points (included for plotting).
  
  if(is.null(minPts.start))
  {
    minPts.start <- ncol(data)+1
  }
  
  
  n <- seq(eps.start, eps.end, step.size)
  
  if(is.null(minPts.end))
  {
    minPts <- rep(minPts.start, times=length(n), each=1)
  }
  
  if(!is.null(minPts.end))
  {
    minPts <- rep(minPts.start:minPts.end, each=length(n))
  }
  
  n <- rep(n, length(minPts.start:minPts.end))
  
  results <- data.frame(eps=n, minPts=as.factor(minPts), n.clust=vector(length = length(n)), noise.pts=vector(length=length(n)))
  
  for(i in 1:length(n))
  {
    scan.res <- dbscan(x=data, eps = n[i], minPts=minPts[i], borderPoints = T)
    results[i,]$n.clust <- length(unique(scan.res[which(scan.res$cluster != 0)]$cluster))
    results[i,]$noise.pts <- as.numeric(table(scan.res$cluster == 0)[2])
  }
  results$log.noise.pts <- log(results$noise.pts)
  return(results)
}