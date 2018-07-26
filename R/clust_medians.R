clust.medians <- function(x, markers, clust.col, noise.clust.id=NULL)
{
  #' @title Cluster Medians
  #' @author Julian Spagnuolo
  #' @export

  if(!is.null(noise.clust.id))
  {
    clusts <- unique(x[which(x[,clust.col] != noise.clust.id),clust.col])
  }
  else{
    clusts <- unique(x[,clust.col])
  }

  clust.defs <- matrix(nrow=length(clusts), ncol=length(markers), dimnames = list(c(clusts), c(markers)))

  for(i in 1:length(clusts))
  {
    for(n in 1:length(markers))
    {
      # Some clusters have only one member, thus the skew-normal location (median) function (snlocation) will not work
      if(nrow(x[which(x[,clust.col] == clusts[i]),]) > 1)
      {
        clust.defs[clusts[i],markers[n]] <- snlocation(markdat = x[which(x[,clust.col] == clusts[i]),markers[n]])
      }
      else{
        clust.defs[clusts[i],markers[n]] <- x[which(x[,clust.col] == clusts[i]),markers[n]]
      }
    }
  }

  return(clust.defs)
}
