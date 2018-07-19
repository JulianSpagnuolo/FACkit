peaks.simple <- function(data, markers, npeaks, threshold, DensityThreshold, max.iter=10)
{
  data <- data[,markers]
  peak.mat <- matrix(nrow=npeaks, ncol = length(markers), dimnames = list(c(0,1), markers))

  if(DensityThreshold != length(markers))
  {
    DensityThreshold <- rep(DensityThreshold, length(markers))
  }

  for(i in 1:length(markers))
  {
    cat("\nProcessing", markers[i])
    iter <- 0

    fudge <- 0
    peaks <- peaksNvalleys(data=data[,i], minDensityThreshold=DensityThreshold[i], gridsize=14000, fudge=fudge)

    if(length(peaks$peaks) < 1)
    {
      cat("\n No peaks found for ", i, "!!! \nTry changing gridsize or minDensityThreshold parameters")
    }
    ## Incrementally increase smoothing of the distribution to find expected number of peaks
    if(length(peaks$peaks) != as.integer(npeaks[i]))
    {
      cat("\nFudging bandwidth parameter to increase smoothing")
    }
    while(as.integer(npeaks[i]) != length(peaks$peaks) & iter <= max.iter)
    {
      fudge <- fudge + 0.1
      peaks <- peaksNvalleys(data=data[,i], minDensityThreshold=DensityThreshold[i], gridsize=14000, fudge=fudge)
      iter <- iter + 1
    }
    cat("\nFudge Factor:",fudge)
    if(all(sign(peaks$dens$x[peaks$peaks]) == c(-1,1)) == FALSE)
    {
      warning("\nFailed to adequately find peaks for ",markers[i], immediate. = T)
    }


    peak.mat[,i] <- peaks$dens$x[peaks$peaks]
  }
  return(peak.mat)
}
