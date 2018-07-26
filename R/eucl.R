eucl <- function(x, y)
{
  #' @title Euclidean Distance of a Point to a Vector
  #' @author Julian Spagnuolo
  #' @param x a numeric vector representing the origin from which to calculate the distance to y
  #' @param y a data.frame or matrix representing the points to which you want to calculate the distance from x. ncol(y) must equal length(x).
  #'
  #'
  #' @export
  d <- vector(length=nrow(y))
  for(i in 1:nrow(y))
  {
    d[i] <- sqrt(sum((x-y[i,])^2))
  }
  return(d)
}
