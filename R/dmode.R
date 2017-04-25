dmode <- function(x, gridsize, fudge=0, ...) {
  den <- KernSmooth::bkde(x, kernel="normal", bandwidth=bw.select(x, gridsize=gridsize, ...)+fudge, gridsize=gridsize)
  ( den$x[den$y==max(den$y)] )
}
