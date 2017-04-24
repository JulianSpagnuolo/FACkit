bw.select <- function (x, scalest = "mad", level = 2L, kernel = "normal",
          canonical = FALSE, gridsize = 401L, range.x = range(x),
          truncate = TRUE)
{
  dmode <- function(x, ...) {
    den <- KernSmooth::bkde(x,kernel="normal", bandwidth=KernSmooth::dpik(x, gridsize=gridsize), gridsize=gridsize)
    ( den$x[den$y==max(den$y)] )
  }

  if (level > 5L)
    stop("Level should be between 0 and 5")
  kernel <- match.arg(kernel, c("normal", "box", "epanech",
                                "biweight", "triweight"))
  del0 <- if (canonical)
    1
  else switch(kernel, normal = 1/((4 * pi)^(1/10)), box = (9/2)^(1/5),
              epanech = 15^(1/5), biweight = 35^(1/5), triweight = (9450/143)^(1/5))
  n <- length(x)
  M <- gridsize
  a <- range.x[1L]
  b <- range.x[2L]
  gpoints <- seq(a, b, length = M)
  gcounts <- KernSmooth:::linbin(x, gpoints, truncate)
  scalest <- match.arg(scalest, c("minim", "stdev", "iqr","mad"))
  scalest <- switch(scalest, stdev = sqrt(var(x)),
                    iqr = (quantile(x, 3/4) - quantile(x, 1/4))/1.349,
                    minim = min((quantile(x, 3/4) - quantile(x, 1/4))/1.349, sqrt(var(x))),
                    mad = mad(x, center=dmode(x)))
  if (scalest == 0)
    stop("scale estimate is zero for input data")
  sx <- (x - mean(x))/scalest
  sa <- (a - mean(x))/scalest
  sb <- (b - mean(x))/scalest
  gpoints <- seq(sa, sb, length = M)
  gcounts <- KernSmooth:::linbin(sx, gpoints, truncate)
  psi4hat <- if (level == 0L)
    3/(8 * sqrt(pi))
  else if (level == 1L) {
    alpha <- (2 * (sqrt(2))^7/(5 * n))^(1/7)
    bkfe(gcounts, 4L, alpha, range.x = c(sa, sb), binned = TRUE)
  }
  else if (level == 2L) {
    alpha <- (2 * (sqrt(2))^9/(7 * n))^(1/9)
    psi6hat <- KernSmooth::bkfe(gcounts, 6L, alpha, range.x = c(sa,
                                                    sb), binned = TRUE)
    alpha <- (-3 * sqrt(2/pi)/(psi6hat * n))^(1/7)
    KernSmooth::bkfe(gcounts, 4L, alpha, range.x = c(sa, sb), binned = TRUE)
  }
  else if (level == 3L) {
    alpha <- (2 * (sqrt(2))^11/(9 * n))^(1/11)
    psi8hat <- KernSmooth::bkfe(gcounts, 8L, alpha, range.x = c(sa,
                                                    sb), binned = TRUE)
    alpha <- (15 * sqrt(2/pi)/(psi8hat * n))^(1/9)
    psi6hat <- KernSmooth::bkfe(gcounts, 6L, alpha, range.x = c(sa,
                                                    sb), binned = TRUE)
    alpha <- (-3 * sqrt(2/pi)/(psi6hat * n))^(1/7)
    KernSmooth::bkfe(gcounts, 4L, alpha, range.x = c(sa, sb), binned = TRUE)
  }
  else if (level == 4L) {
    alpha <- (2 * (sqrt(2))^13/(11 * n))^(1/13)
    psi10hat <- KernSmooth::bkfe(gcounts, 10L, alpha, range.x = c(sa,
                                                      sb), binned = TRUE)
    alpha <- (-105 * sqrt(2/pi)/(psi10hat * n))^(1/11)
    psi8hat <- KernSmooth::bkfe(gcounts, 8L, alpha, range.x = c(sa,
                                                    sb), binned = TRUE)
    alpha <- (15 * sqrt(2/pi)/(psi8hat * n))^(1/9)
    psi6hat <- KernSmooth::bkfe(gcounts, 6L, alpha, range.x = c(sa,
                                                    sb), binned = TRUE)
    alpha <- (-3 * sqrt(2/pi)/(psi6hat * n))^(1/7)
    KernSmooth::bkfe(gcounts, 4L, alpha, range.x = c(sa, sb), binned = TRUE)
  }
  else if (level == 5L) {
    alpha <- (2 * (sqrt(2))^15/(13 * n))^(1/15)
    psi12hat <- KernSmooth::bkfe(gcounts, 12L, alpha, range.x = c(sa,
                                                      sb), binned = TRUE)
    alpha <- (945 * sqrt(2/pi)/(psi12hat * n))^(1/13)
    psi10hat <- KernSmooth::bkfe(gcounts, 10L, alpha, range.x = c(sa,
                                                      sb), binned = TRUE)
    alpha <- (-105 * sqrt(2/pi)/(psi10hat * n))^(1/11)
    psi8hat <- KernSmooth::bkfe(gcounts, 8L, alpha, range.x = c(sa,
                                                    sb), binned = TRUE)
    alpha <- (15 * sqrt(2/pi)/(psi8hat * n))^(1/9)
    psi6hat <- KernSmooth::bkfe(gcounts, 6L, alpha, range.x = c(sa,
                                                    sb), binned = TRUE)
    alpha <- (-3 * sqrt(2/pi)/(psi6hat * n))^(1/7)
    KernSmooth::bkfe(gcounts, 4L, alpha, range.x = c(sa, sb), binned = TRUE)
  }
  scalest * del0 * (1/(psi4hat * n))^(1/5)
}
