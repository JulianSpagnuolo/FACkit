cutfind2 <- function(x, markers, npeaks=NULL, DensityThreshold=NULL, gridsize=14000, max.restarts=1, epsilon=1e-9, maxit=8000, seed=42, auto=TRUE, metric="AIC") {

#'  @author Julian Spagnuolo
#'
#'  @param x data.frame or named matrix containing raw data for which you wish to identify cutoff values for subsequent transformation.
#'  @param markers a vector containing the columns in x that you wish to parse through cutfind2
#'  @param npeaks a named vector of integers indicating the number of peaks in each distribution that you wish to identify. This can contain NULL values. See details
#'  @param DensityThreshold Named numeric vector indicating the density above which peaks will be identified. Set this below the smallest peak in the density distribution of your data.
#'  @param gridsize gridsize used in the density estimation steps - this needn't be changed.
#'  @param epsilon accuracy of mixture modelling step (only used if npeaks for marker is > 1)
#'  @param maxit maximum iterations used by mixture modelling step (only used if npeaks for marker is > 1). If this step fails, it will restart (number of restarts set by max.restarts) and add 2000 iterations to maxit.
#'  @param max.restarts maximum restarts allowed in the mixture modelling steps (increases maxit by 2000 each restart). The function will use the parameters returned by the mixture modelling step to re-seed the algorithm.
#'  @param seed This is the system seed used by the mixture modelling step, just to increase the reproducibility of the function
#'  @param auto Logical - if TRUE the function will search for the optimal number of components (distributions) to fit in the mixture model by searching a range of +/- 2 peaks from the number of peaks identified.
#'  @param metric vector - the metric used to choose the optimum mixture model to fit to the data choose one of "BIC", "AIC", or "ICL"
#'

  if(is.null(npeaks))
  {
    npeaks <- rep(NULL, length(markers))
    names(npeaks) <- markers
  }

  # Set universal Density Threshold
  if(length(DensityThreshold) == 1)
  {
    DensityThreshold <- rep(DensityThreshold, length(markers))
    names(DensityThreshold) <- markers
  }

  if(is.null(DensityThreshold))
  {
    cat("Provide a vector of Density Thresholds for each marker (or pick a point below the lowest peak you wish to identify)")
    break()
  }

  if(length(DensityThreshold) != length(markers))
  {
    cat("Missing a Density Threshold in your vector")
    break()
  }

  cutoffs <- vector(mode="list", length=length(markers))
  names(cutoffs) <- markers
  itmax <- maxit
  for(i in markers)
  {
    cat("\nProcessing", i)

    fudge <- 0
    peaks <- peaksNvalleys(data=x[,i], minDensityThreshold=DensityThreshold[i], gridsize=gridsize, fudge=fudge)
    if(length(peaks$peaks) < 1)
    {
      cat("\n No peaks found for ", i, "!!! \nTry changing gridsize or minDensityThreshold parameters")
    }
    ## Incrementally increase smoothing of the distribution to find expected number of peaks
    if(!is.null(npeaks[i]))
    {
      while(length(peaks$peaks) != (npeaks[i]))
      {
        fudge <- fudge + 0.1
        peaks <- peaksNvalleys(data=x[,i], minDensityThreshold=DensityThreshold[i], gridsize=gridsize, fudge=fudge)
      }
    }

    # Find the summary stats for cutoff values

    ## For unimodal distributions find mode-centered mad.
    if(length(peaks$peaks) == 1 )
    {
      cat("\n Found", length(peaks$peaks), " peak in ", i, " distribution")
      cutoffs[[i]]$right  <- dmode(x[,i], gridsize=gridsize)+2*mad(x=x[,i],center=dmode(x[,i], gridsize=gridsize, fudge=fudge))
      cutoffs[[i]]$left  <- dmode(x[,i], gridsize=gridsize)-2*mad(x=x[,i],center=dmode(x[,i], gridsize=gridsize, fudge=fudge))
      cutoffs[[i]]$peak <-  peaks$dens$x[peaks$peaks]
      cutoffs[[i]]$sigma <- mad(x=x[,i],center=dmode(x[,i], gridsize=gridsize))
      cat("\n", i, " done!")
    }
    if(length(peaks$peaks) > 1)
    {
      cat("\n Found ", length(peaks$peaks), " peaks in ",i, " distribution" )
      cat("\n Mixture modelling ", i, " marker \n")
      k <- length(peaks$peaks)

      ## Fit a miture model to the data
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

      # Get restart params if modelling failed to converge in maxit
      if(mixmod$error == 1)
      {
        pro <- mixmod$pro
        mu <- mixmod$mu
        sigma <- mixmod$sigma
        dof <- mixmod$dof
        delta <- mixmod$delta
      }
      # Restart mixture modelling until it converges
      while(error == 1 & restarts <= max.restarts)
      {
        set.seed(42)
        mixmod <- EMMIXskew::EmSkew(dat=as.matrix(x[,i]), g=k, itmax=itmax, epsilon=epsilon, distr="mst", debug=F, init=list(pro, mu, sigma, dof, delta))
        error <- mixmod$error
        if(mixmod$error == 1)
        {
          itmax <- itmax + 2000
          restarts <- restarts + 1
          pro <- mixmod$pro
          mu <- mixmod$mu
          sigma <- mixmod$sigma
          dof <- mixmod$dof
          delta <- mixmod$delta
          warning("Model failed to converge, increasing max allowed iterations and restarting\n", immediate.=TRUE)
        }
      }
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
