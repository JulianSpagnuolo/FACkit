med.cov <- function(expdata, markers, use.median=TRUE){
  #' @title Median Covariance
  #' @author Julian Spagnuolo
  #'
  #' @param expdata data.frame. Contains expression data
  #' @param markers Character vector, identifies vectors in data containing the marker expression data.
  #' @param use.median Logical. If TRUE, the median of columns in expdata is used to calculate the covariance. If FALSE, the mean is used (identical to stats::cov). Default is TRUE.
  #'
  #'
  #' @export

  n <- nrow(expdata) - 1L

  if(!isTRUE(use.median)){
    mu <- colMeans(expdata[,markers])
    cov <- sweep(as.matrix(expdata[,markers]), 2L, mu)
  }else{
    meds <- apply(expdata[,markers], MARGIN = 2, function(x) {median(x)})
    cov <- sweep(as.matrix(expdata[,markers]), 2L, meds)
  }

  cov <- crossprod(cov) / n

  return(cov)
}
