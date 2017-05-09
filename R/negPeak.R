neg.peak <- function(data, markers, thresh)
{
  neg.peaks <- vector(length=length(markers))
  names(neg.peaks) <- markers
  for(i in markers)
  {
    x <- peaksNvalleys(data=data[,i], minDensityThreshold = thresh, gridsize = 14000, fudge=0)

    # take the peak corresponding to the lowest expression.
    neg.peaks[i] <- min(x$dens$x[x$peaks])
  }
  return(neg.peaks)
}
