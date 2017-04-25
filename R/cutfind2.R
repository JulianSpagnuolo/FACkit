cutfind2 <- function(x, markers, npeaks=NULL, DensityThreshold=NULL, gridsize=14000, max.restarts=1, epsilon=1e-9, maxit=8000, seed=42, auto=TRUE, metric="AIC") {

#'  @title Cutfind
#'  @description finds summary stats for distributions by either mixture modelling or finding the mode and mad of unimodal distributions
#'  @author Julian Spagnuolo
#'
#'  @param x data.frame or named matrix containing raw data for which you wish to identify cutoff values for subsequent transformation.
#'  @param markers a vector containing the columns in x that you wish to parse through cutfind2
#'  @param npeaks a named vector of integers indicating the number of peaks in each distribution that you wish to identify. This can contain NA values. See details
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
    npeaks <- vector(length=length(markers))
    npeaks <- rep(NA, length(markers))
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
    cat("\nLooking for ", as.integer(npeaks[i])," peaks")

    fudge <- 0
    peaks <- peaksNvalleys(data=x[,i], minDensityThreshold=DensityThreshold[i], gridsize=gridsize, fudge=fudge)

    if(length(peaks$peaks) < 1)
    {
      cat("\n No peaks found for ", i, "!!! \nTry changing gridsize or minDensityThreshold parameters")
    }
    ## Incrementally increase smoothing of the distribution to find expected number of peaks
    if(!is.na(npeaks[i]))
    {
      if(as.integer(npeaks[i]) != length(peaks$peaks))
      {
        cat("\nFudging bandwidth parameter to increase smoothing")
      }
      while(as.integer(npeaks[i]) != length(peaks$peaks))
      {
        fudge <- fudge + 0.1
        peaks <- peaksNvalleys(data=x[,i], minDensityThreshold=DensityThreshold[i], gridsize=gridsize, fudge=fudge)
      }
      cat("\nFudge Factor:",fudge)
    }

    # Find the summary stats for cutoff values

    ## For unimodal distributions find mode-centered mad.
    if(length(peaks$peaks) == 1 )
    {
      cat("\n Found", length(peaks$peaks), " peak in ", i, " distribution")
      cutoffs[[i]]$right  <- dmode(x[,i], gridsize=gridsize, fudge=fudge, bw.method="bw.select")+2*mad(x=x[,i],center=dmode(x[,i], gridsize=gridsize, fudge=fudge, bw.method="bw.select"))
      cutoffs[[i]]$left  <- dmode(x[,i], gridsize=gridsize, fudge=fudge, bw.method="bw.select")-2*mad(x=x[,i],center=dmode(x[,i], gridsize=gridsize, fudge=fudge, bw.method="bw.select"))
      cutoffs[[i]]$peak <-  peaks$dens$x[peaks$peaks]
      cutoffs[[i]]$sigma <- mad(x=x[,i],center=dmode(x[,i], gridsize=gridsize, fudge=fudge, bw.method="bw.select"))
      cat("\n", i, " done!")
    }
    if(length(peaks$peaks) > 1)
    {
      cat("\nFound ", length(peaks$peaks), " peaks in ",i, " distribution" )
      cat("\nMixture modelling ", i, " marker \n")
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
        comps <- data.frame(row.names=seq(from=k-1, to=k+1, by=1),
                            BIC=vector(length=length(seq(from=k-1, to=k+1, by=1))),
                            AIC=vector(length=length(seq(from=k-1, to=k+1, by=1))),
                            ICL=vector(length=length(seq(from=k-1, to=k+1, by=1))))
        for(l in 1:length(row.names(comps)))
        {
          mod <- EMMIXskew::EmSkew(dat=as.matrix(x[,i]), g=as.integer(row.names(comps)[l]), itmax=500, debug=F, distr="mst", epsilon=1e-5)
          comps[l,]$BIC <- mod$bic
          comps[l,]$AIC <- mod$aic
          comps[l,]$ICL <- mod$ICL
        }
        ### Choose the simplest model if BIC/AIC for 2 models are close enough.
        if(metric != "ICL")
        {
          for(l in 1:nrow(comps)-1)
          {
            mod.dif <- (comps[l,metric] - comps[l+1,metric])/comps[l,metric]
            if(mod.dif < 5e-5*comps[l,metric])
            {
              k <- as.integer(row.names(comps[l,]))
              cat("Optimum number of components is ", k,"\n")
            }
          }
        }
        if(metric == "ICL")
        {
          k <- as.integer(row.names(comps[which(comps[,metric] == max(comps[,metric])),]))
          cat("Optimum number of components is ", row.names(comps[which(comps[,metric] == max(comps[,metric])),]),"\n")
        }

        #if(metric != "ICL")
        #{
        #  k <- as.integer(row.names(comps[which(comps[,metric] == min(comps[,metric])),]))
        #}
        #else if(metric == "ICL")
        #{
        #  k <- as.integer(row.names(comps[which(comps[,metric] == max(comps[,metric])),]))
        #}
      }
      set.seed(42)
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
      restarts <- 0
      maxit <- itmax + 2000
      while(mixmod$error == 1 & restarts <= max.restarts)
      {
        cat(" Model failed to converge, restarting with params from initial run\n")
        restart <- restart + 1
        set.seed(42)
        mixmod <- EMMIXskew::EmSkew(dat=as.matrix(x[,i]), g=k, itmax=maxit, epsilon=epsilon, distr="mst", debug=F, init=list(pro, mu, sigma, dof, delta))
        # get restart params
        pro <- mixmod$pro
        mu <- mixmod$mu
        sigma <- mixmod$sigma
        dof <- mixmod$dof
        delta <- mixmod$delta

        if(mixmod$error == 1)
        {
          maxit <- maxit + 2000
          cat(" Model failed to converge, increasing max allowed iterations and restarting\n")
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
