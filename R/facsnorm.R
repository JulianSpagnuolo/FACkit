facsnorm <- function(x, cutoffs, asinCofac=25, method=c("arcsin","loglin","logit"))
{
  #' @title Facsnorm
  #' @author Julian Spagnuolo, Tobias Rutishauser
  #' @param x data.frame or matrix of raw data
  #' @param cutoffs vector of cutoff values. values MUST be in the order of the markers in x that they are to be applied to
  #' @param asinCofac spreading parameter for arcsin transformation. Only needs to be set if method == arcsin. Default == 25.
  #' @param method Type of transformation to be applied, either arcsin ("arcsin"), log-linear ("loglin") or logit ("logit") transformation
  if(method == "arcsin")
  {
    cat("Using arcsin Transformation\n")
      if(is.null(names(cutoffs)))
      {
        names(cutoffs) <- names(x[1:ncol(x)])
      }
      asinData <- asinh(sweep(x[1:ncol(x)],2,cutoffs)/asinCofac)
      return(out <- data.frame(asinData))
  }
  if(method == "loglin")
  {
    cat("Using Logical Transformation\n")
    for(i in 1:ncol(x))
    {
      if(sign(min(x[,i])) == -1)
      {
        x[,i] <- x[,i]-min(x[,i])
      }
      if(sign(min(x[,i])) == 1)
      {
        x[,i] <- x[,i]+min(x[,i])
      }
    }
    x <- sweep(x, 2, STATS = 1, FUN = "+")
    return(log10(x))
  }
  if(method == "logit")
  {
    cat("Using logit transformation\n")
    # Linearise the data with 1 as minimum (avoids creating infinite values when log transforming)
    for(i in 1:ncol(x))
    {
      if(sign(min(x[,i])) == -1)
      {
        x[,i] <- x[,i]-min(x[,i])
      }
      if(sign(min(x[,i])) == 1)
      {
        x[,i] <- x[,i]+min(x[,i])
      }
    }
    # add 1 to each value to avoid log-transforming 0 (gives and infinite number)
    x <- sweep(x, 2, STATS = 1, FUN = "+")
    # put data on scale between 0 and 1 and perform logit transformation
    for(i in 1:ncol(x))
    {
      x[,i] <- x[,i]/max(x[,i])
      x[,i] <- log(x[,i]/1-x[,i])
    }
    return(x)
  }
}
