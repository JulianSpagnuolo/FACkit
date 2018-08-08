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

  if(!isTRUE(use.median)){
    cov <- crossprod(sweep(as.matrix(expdata[,markers]), 2L, colMeans(expdata[,markers]))) / (nrow(expdata[,markers]) - 1L)
  }else{
    meds <- apply(expdata[,markers], MARGIN = 2, function(x) {median(x)})
    cov <- crossprod(sweep(as.matrix(expdata[,markers]), 2L, meds)) / (nrow(expdata[,markers]) - 1L)
  }
  return(cov)
}
