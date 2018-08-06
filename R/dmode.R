dmode <- function(x, gridsize, fudge, bw.method=c("bw.select","dpik"), ...)
{
  #' @title dmode
  #' @author Julian Spagnuolo
  #'
  #' @param x vector or column of data frame in which to find the value at which the density of the distribution is highest
  #' @param gridsize gridsize
  #' @param fudge a cofactor of the bandwidth, only used in cases where distributions need to be more smooth than returned by default
  #' @param bw.method default is "bw.select". This chooses the bandwidth selection method. The only other choose is "dpik". See details
  #' @param ... other parameters to be passed to the bandwidth selection function
  #'
  #' @export
  #'
  #'

  if(bw.method == "bw.select")
  {
    den <- KernSmooth::bkde(x, kernel="normal", bandwidth=bw.select(x, gridsize=gridsize, ...)+fudge, gridsize=gridsize)
    ( den$x[den$y==max(den$y)] )
  }
  else if(bw.method == "dpik")
  {
    den <- KernSmooth::bkde(x, kernel="normal", bandwidth=KernSmooth::dpik(x=x, gridsize=gridsize, ...), gridsize=gridsize)
    ( den$x[den$y==max(den$y)] )
  }
}
