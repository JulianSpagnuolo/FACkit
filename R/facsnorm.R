facsnorm <- function(x, cutoffs, asinCofac=25, method=c("arcsin","logical"))
{
  #' @title Facsnorm
  #' @author Julian Spagnuolo, Tobias Rutishauser
  #' @param x data.frame or matrix of raw data
  #' @param cutoffs vector of cutoff values. If the vector is unnamed the values MUST be in the order of the markers in x that they are to be applied to
  #' @param asinCofac spreading parameter for arcsin transformation. Only needs to be set if method == arcsin. Default == 25.
  #' @param method Type of transformation to be applied, either arcsin or logical transformation
  if(method == "arcsin")
  {
      if(is.null(names(cutoffs)))
      {
        names(cutoffs) <- names(x[1:ncol(x)])
      }
      asinData <- asinh(sweep(x[1:ncol(x)],2,cutoffs)/asinCofac)
      return(out <- data.frame(asinData))
  }
  if(method == "logical")
  {
    for(i in 1:ncol(x))
    {
      if(sign(min(x[i])) == -1)
      {
        x[i] <- x[i]-(min(x[i])+1)
        x[i] <- log10(x[i])
      }
      if(sign(min(x[i])) == 1)
      {
        x[i] <- x[i]+(min(x[i])+1)
        x[i] <- log10(x[i])
      }
    }
    return(x)
  }

}
