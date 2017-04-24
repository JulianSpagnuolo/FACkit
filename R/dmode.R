dmode <- function(x, gridsize, ...) {
  den <- KernSmooth::bkde(x, kernel="normal", bandwidth=bw.select(x, gridsize=gridsize, ...), gridsize=gridsize)
  ( den$x[den$y==max(den$y)] )
}
