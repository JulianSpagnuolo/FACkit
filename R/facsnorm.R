facsnorm <- function(x, cutoffs, asinCofac=25) {for (i in 1:length(x)) {
  names(cutoffs) <- names(x[1:ncol(x)])
  asinData <- asinh(sweep(x[1:ncol(x)],2,cutoffs)/asinCofac)
  return(out <- data.frame(asinData))
}
}
