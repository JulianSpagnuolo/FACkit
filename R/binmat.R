binmat <- function(data, cluster.col, markers, split.list, thresh = 0)
{
  #' @import gtools
  
  bin.list <- vector(mode="list", length=length(split.list))
  names(bin.list) <- names(split.list)
  for(i in 1:length(split.list))
  {
    splitters <- split.list[[i]]
    
    bin.defs <- matrix(nrow=2, ncol=length(splitters), dimnames = list(c("0","1"), c(splitters)))
    
    not.split.defs <- matrix(nrow=2, ncol=length(markers[!(markers %in% splitters)]),
                             dimnames=list(c("0","1"),c(markers[!(markers %in% splitters)])))
    
    for(n in 1:length(markers))
    {
      sub.data <- data[which(data[,cluster.col] == names(split.list)[i]),markers[n]]
      
      if(markers[n] %in% splitters)
      {
        bin.defs["0",markers[n]] <- snlocation(markdat=sub.data[which(sub.data < thresh)])
        bin.defs["1",markers[n]] <- snlocation(markdat=sub.data[which(sub.data > thresh)])
      }
      else{
        const.loc <- snlocation(markdat=sub.data)
        if(sign(const.loc) == 1)
        {
          not.split.defs["1",markers[n]] <- const.loc
          not.split.defs["0",markers[n]] <- NA
        }
        else{
          not.split.defs["1",markers[n]] <- NA
          not.split.defs["0",markers[n]] <- const.loc
        }
        rm(const.loc)
      }
    }
    bin.defs <- cbind(bin.defs,not.split.defs)
    bin.defs <- bin.defs[,markers] ## reordering
    
    rm(not.split.defs)
    rm(sub.data)
    
    cat("\nSetting up bins ", i)
    bin.mat <- permutations(n=2, r=length(splitters), v=c(0,1), repeats.allowed=T)
    colnames(bin.mat) <- splitters
    for(n in splitters)
    {
      bin.mat[which(bin.mat[,n] == 0),n] <- bin.defs["0",n]
      bin.mat[which(bin.mat[,n] == 1),n] <- bin.defs["1",n]
    }
    
    not.split.bins <- matrix(nrow=nrow(bin.mat), ncol=length(markers[!(markers %in% splitters)]),
                             dimnames=list(c(1:nrow(bin.mat)),c(markers[!(markers %in% splitters)])))
    for(n in markers[!(markers %in% splitters)])
    {
      not.split.bins[,n] <- bin.defs[!(is.na(bin.defs[,n])), n]
    }
    
    bin.mat <- cbind(bin.mat, not.split.bins)
    rm(not.split.bins)
    
    bin.mat <- bin.mat[,markers]
    #nbins <- nrow(bin.mat)
    
    
    bin.list[[i]] <- bin.mat
    rm(bin.mat)
  }
  
  return(bin.list)
}

