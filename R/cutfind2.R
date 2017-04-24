cutfind2 <- function(markers, x, DensityThreshold=NULL, gridsize=14000, max.restarts=1, epsilon=1e-9, maxit=8000, seed=42, auto=TRUE, metric="AIC") {
  require(KernSmooth)
  require(EMMIXskew)
  
  if(length(DensityThreshold) == 1){
    DensityThreshold <- rep(DensityThreshold, length(markers))
    names(DensityThreshold) <- markers
  }
  
  if(is.null(DensityThreshold)){
    cat("Provide a vector of Density Thresholds for each marker (or pick a point below the lowest peak you wish to identify)")
    break()
  }
  
  if(length(DensityThreshold) == length(markers)) {
    names(DensityThreshold) <- markers
  }
  
  if(length(DensityThreshold) != length(markers)) {
    cat("Missing a Density Threshold in your vector")
    break()
  }
  
  # Internal Function 
  ### peaksNvalleys to find peaks and valleys (duh)
  peaksNvalleys <- function (x, minDensityThreshold = 0.05, gridsize, ...) 
  {
    dens <- KernSmooth::bkde(x, bandwidth=bw.select(x, scalest="mad",gridsize=gridsize)+10, gridsize=gridsize)
    secondDerivative <- diff(sign(diff(dens$y)))
    peaks <- which(secondDerivative == -2)
    peaks <- peaks[dens$y[peaks] > minDensityThreshold]
    tmp <- split(dens$y, rep(c(peaks, 0), c(peaks[1], diff(peaks), 
                                            length(dens$y) - max(peaks))))
    valleys <- unlist(lapply(tmp[-c(1, 2)], function(l)
    {which.min(l)[1]}))+peaks[-length(peaks)]
    list(dens = dens, peaks = peaks, valleys = valleys)
  }
  ### dmode to find density peaks
  dmode <- function(x, gridsize, ...) {
    den <- KernSmooth::bkde(x, kernel="normal", bandwidth=bw.select(x, gridsize=gridsize, ...)+10, gridsize=gridsize)
    ( den$x[den$y==max(den$y)] )
  }  
  ### bw.select to get density bandwidth
  bw.select <- function (x, scalest = "mad", level = 2L, kernel = "normal", 
                         canonical = FALSE, gridsize = 401L, range.x = range(x), 
                         truncate = TRUE) 
  {
    require(KernSmooth)
    
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
  cutoffs <- vector(mode="list", length=length(markers))
  names(cutoffs) <- markers
  itmax <- maxit
  for(i in markers){
    cat("\nProcessing", i)
    peaks <- peaksNvalleys(x[,i], minDensityThreshold=DensityThreshold[i], gridsize=gridsize)
    if(length(peaks$peaks) < 1){
      cat("\n No peaks found for ", i, "!!! \nTry changing gridsize or minDensityThreshold parameters")
    }
    if(length(peaks$peaks) == 1 ){
      cat("\n Found", length(peaks$peaks), " peak in ", i, " distribution")
      cutoffs[[i]]$right  <- dmode(x[,i], gridsize=gridsize)+2*mad(x=x[,i],center=dmode(x[,i], gridsize=gridsize))
      cutoffs[[i]]$left  <- dmode(x[,i], gridsize=gridsize)-2*mad(x=x[,i],center=dmode(x[,i], gridsize=gridsize))
      cutoffs[[i]]$peak <-  peaks$dens$x[peaks$peaks]
      cutoffs[[i]]$sigma <- mad(x=x[,i],center=dmode(x[,i], gridsize=gridsize))
      cat("\n", i, " done!")
    }
    if(length(peaks$peaks) > 1){
      cat("\n Found ", length(peaks$peaks), " peaks in ",i, " distribution" )
      cat("\n Mixture modelling ", i, " marker \n")
      k <- length(peaks$peaks)
      
      # create loop to make sure the Modelling function converges
      error <- 1
      restarts <- 1
      set.seed(seed=seed)
      
      ## Find optimal number of components to fit
      ### Make this run in parallel
      if(auto == TRUE)
      {
        cat("Finding optimum number of components using ", metric, "\n")
        if(k > 2)
        {
          comps <- data.frame(row.names=seq(from=k-2, to=k+2, by=1),
                              BIC=vector(length=length(seq(from=k-2, to=k+2, by=1))),
                              AIC=vector(length=length(seq(from=k-2, to=k+2, by=1))),
                              ICL=vector(length=length(seq(from=k-2, to=k+2, by=1))))
        }
        if(k == 2)
        {
          comps <- data.frame(row.names=seq(from=k-1, to=k+2, by=1),
                              BIC=vector(length=length(seq(from=k-1, to=k+2, by=1))),
                              AIC=vector(length=length(seq(from=k-1, to=k+2, by=1))),
                              ICL=vector(length=length(seq(from=k-1, to=k+2, by=1))))
        }
        for(l in 1:length(row.names(comps)))
        {
          mod <- EmSkew(dat=as.matrix(x[,i]), g=as.integer(row.names(comps)[l]), itmax=500, debug=F, distr="mst", epsilon=1e-5)
          comps[l,]$BIC <- mod$bic
          comps[l,]$AIC <- mod$aic
          comps[l,]$ICL <- mod$ICL
        }
        k <- as.integer(row.names(comps[which(comps[,metric] == min(comps[,metric])),]))
        cat("Optimum number of components is ", row.names(comps[which(comps[,metric] == min(comps[,metric])),]),"\n")
      }
      
      mixmod <- EmSkew(dat=as.matrix(x[,i]), g=k, itmax=itmax, epsilon=epsilon, distr="mst", debug=F)
      #while(error == 1 & restarts <= max.restarts)
      #{
      #  set.seed(42)
      #  mixmod <- EmSkew(dat=as.matrix(x[,i]), g=k, itmax=itmax, epsilon=epsilon, distr="mst", debug=F)
      #  error <- mixmod$error
      #  if(mixmod$error == 1)
      #  {
      #    itmax <- itmax + 2000
      #    restarts <- restarts + 1
      #    warning("Model failed to converge, increasing max allowed iterations and restarting\n", immediate.=TRUE)
      #  }
      #}
      ## Fit a miture model to the data
      if(mixmod$error == 0)
      {
        cat("Success! Modeling ", i," as ", length(mixmod$modpts), " skew normal distributions within ",itmax, " iterations\n")
      }
      if(mixmod$error == 1)
      {
        cat("Failed to converge within", itmax,"iterations\n")
      }
      if(mixmod$error == 2)
      {
        cat("Failed to get initial values\n")
      }
      if(mixmod$error == 3)
      {
        cat("Singularity, prepare for the robotic apocalypse\n")
      }
      
      ## Simulate data for individual components
      y <- vector(mode="list", length=length(mixmod$pro))
      for(j in 1:length(mixmod$pro)){
        set.seed(42)
        y[[j]]$model <- rdmst(n=table(mixmod$clust)[j],p=1, mean=mixmod$mu[,j], cov=mixmod$sigma[[j]], del=mixmod$delta[j])
        y[[j]]$left <- mixmod$modpts[j]-2*mad(x=y[[j]]$model, center=mixmod$modpts[j])
        y[[j]]$right <- mixmod$modpts[j]+2*mad(x=y[[j]]$model, center=mixmod$modpts[j])
      }
      cutoffs[[i]]$model <- mixmod
      cutoffs[[i]]$components <- y
      cutoffs[[i]]$comp.opt <- comps
      cat("Finished finding cutoffs for ", i, " marker \n")
    }
    
  }
  cat("\nCutoffs have been FACked!", "\nPlease check accuracy by plotting yo data")
  return(cutoffs)
}