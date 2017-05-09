
peaksNvalleys <- function (data, minDensityThreshold, gridsize, fudge,...)
{
  # From FlowSOM internal function, modified to use bkde from KernSmooth and the modified bw.select function
  dens <- KernSmooth::bkde(data, bandwidth=bw.select(data, scalest="mad",gridsize=gridsize)+fudge, gridsize=gridsize)
  secondDerivative <- diff(sign(diff(dens$y)))
  peaks <- which(secondDerivative == -2)
  peaks <- peaks[dens$y[peaks] > minDensityThreshold]
  tmp <- split(dens$y, rep(c(peaks, 0), c(peaks[1], diff(peaks),length(dens$y) - max(peaks))))
  valleys <- unlist(lapply(tmp[-c(1, 2)], function(l) {which.min(l)[1]})) + peaks[-length(peaks)]
  list(dens = dens, peaks = peaks, valleys = valleys)
}
