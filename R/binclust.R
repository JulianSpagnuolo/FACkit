binclust <- function(data, markers, nbins=100, distmet="eucl", percentiles=FALSE)
{
  #' @author Julian Spagnuolo
  #' @title Bin Clustering
  #' @description Clustering by nearest bin
  #' @param data matrix or data.frame containing the data to be clustered
  #' @param markers character vector of the columns names of data that should be used for clustering
  #' @param nbins integer. the number of bins to be created for each dimension in data
  #' @param distmet character vector. the distance metric to be used for matching the closest bin for each marker
  #' @param percentiles logical. Whether to use percentiles to create the bins for each marker. If FALSE, bins will be evenly spaced across the range for each dimension.
  #'
  #'
  #'

  # Create the bins
  bin.mat <- matrix(data=NA, nrow=nbins, ncol=length(markers), dimnames = list(c(1:nbins), markers))
  for(i in markers)
  {
    if(percentiles == FALSE)
    {
      bin.mat[,i] <- seq.int(from=min(data[,i]), to=max(data[,i]), length.out = nbins) # for evenly spaced bins
    }
    if(percentiles == TRUE)
    {
      bin.mat[,i] <- as.vector(quantile(data[,i], seq.int(from=0,to=1, length.out = nbins))) # for bins spaced by the distribution of the data
    }
  }

  for(i in 1:nrow(data)) # iterate through all data points
  {
    # for each datapoint make a container for the nearest bin found
    temp.bin <- vector(length = length(markers))
    names(temp.bin) <- markers
    for(j in markers) # and dimensions
    {
      # find the nearest bin for marker j
      if(distmet == "eucl")
      {
        # Finds the closest j bin for data point i,j
        # Dont need to take sqrt of the sum of squares since this is only a direct point-to-point difference.
        dist.index <- which(abs(bin.mat[,j] - data[i,j]) == min(abs(bin.mat[,j] - data[i,j]))) # this could be improved with a sweep???
        temp.bin[j] <- bin.mat[dist.index,j]
      }
    }

    # initialise the bin list
    binlist <- list()

    if(length(binlist) == 0) # initialise the first entry
    {
      binlist <- list(bin=temp.bin, members=c(i))
    }
    else
    {
      for(n in 1:length(binlist)) # iterate through the binlist to check if the bin is already there, else add it
      {
        if(isTRUE(identical(binlist[[n]]$bin, temp.bin)))
        {
          binlist[[n]]$members <- append(x = binlist[[n]]$members, values = i)
        }
        else
        {
          binlist <- list(binlist, list(bin=temp.bin, members=c(i)))
        }
      }
    }
  }

  return(binlist)
}
