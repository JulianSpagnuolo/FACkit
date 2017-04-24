cutfind2 <- function(markers, x, DensityThreshold=NULL, gridsize=14000, max.restarts=1, epsilon=1e-9, maxit=8000, seed=42, auto=TRUE, metric="AIC") {


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

  cutoffs <- vector(mode="list", length=length(markers))
  names(cutoffs) <- markers
  itmax <- maxit
  for(i in markers){
    cat("\nProcessing", i)
    peaks <- peaksNvalleys(data=x[,i], minDensityThreshold=DensityThreshold[i], gridsize=gridsize, fudge=0.1)
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
          mod <- EMMIXskew::EmSkew(dat=as.matrix(x[,i]), g=as.integer(row.names(comps)[l]), itmax=500, debug=F, distr="mst", epsilon=1e-5)
          comps[l,]$BIC <- mod$bic
          comps[l,]$AIC <- mod$aic
          comps[l,]$ICL <- mod$ICL
        }
        k <- as.integer(row.names(comps[which(comps[,metric] == min(comps[,metric])),]))
        cat("Optimum number of components is ", row.names(comps[which(comps[,metric] == min(comps[,metric])),]),"\n")
      }

      mixmod <- EMMIXskew::EmSkew(dat=as.matrix(x[,i]), g=k, itmax=itmax, epsilon=epsilon, distr="mst", debug=F)
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
        y[[j]]$model <- EMMIXskew::rdmst(n=table(mixmod$clust)[j],p=1, mean=mixmod$mu[,j], cov=mixmod$sigma[[j]], del=mixmod$delta[j])
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
