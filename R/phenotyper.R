phenotyper <- function(data, markers, threshold, prob=0.95, epsilon=1e-6, itmax=2000, minDensityThreshold, verbose=TRUE, get.start=TRUE, max.restarts=1)
{

  ## First determine initial cluster identities.
  init.clust <- matrix(ncol=ncol(data), nrow=nrow(data), dimnames=list(c(1:nrow(data)), c(markers)))
  if(is.list(threshold) == FALSE)
  {
    for(i in markers)
    {
      init.clust[which(data[,i] < threshold[i]),i] <- 1
      init.clust[which(data[,i] > threshold[i]),i] <- 2
    }
  }
  else(is.list(threshold) == TRUE)
  {
    for(i in markers)
    {
      # Bimodal distributions
      if(length(threshold[[i]] == 1))
      {
        init.clust[which(data[,i] < threshold[[i]]),i] <- 1
        init.clust[which(data[,i] > threshold[[i]]),i] <- 2
      }

      # identify clusters for polymodal distributions
      if(length(threshold[[i]]) > 1)
      {
        for(j in 1:length(threshold[[i]]))
        {
          if(j != length(threshold[[i]]))
          {
            init.clust[which(data[,i] > threshold[[i]][j-1] & data[,i] < threshold[[i]][j+1]),i] <- j
            init.clust[which(data[,i] > threshold[[i]][j] & data[,i] < threshold[[i]][j+1]),i] <- j+1
          }
          if(j == length(threshold[[i]]))
          {
            init.clust[which(data[,i] > threshold[[i]][j] & data[,i]),i] <- j+1
          }
        }
      }
    }
  }

  # Fit models
  models <- vector(mode="list", length=length(markers))
  for(i in markers)
  {

    # Getting initial parameters

    ## Search for peaks, if peaks greater than expected, increase the amount of smoothing using the fudge factor.
    ### This will fail if length peaks is truely greater than expected and the user is forcing a smaller number of peaks
    fudge=0
    peaks <- peaksNvalleys(data=data[,i], minDensityThreshold=minDensityThreshold, gridsize=14000, fudge=fudge)
    g <- 1
    while(length(peaks$peaks) != length(threshold[[i]])+1 & fudge < 1)
    {
      fudge <- fudge+0.1
      peaks <- peaksNvalleys(data=data[,i], minDensityThreshold=minDensityThreshold, gridsize=14000, fudge=fudge)
      g <- length(peaks$peaks)
    }
    if(fudge == 1 & length(peaks$peaks) != length(threshold[[i]])+1)
    {
      g <- 1
      peaks <- peaksNvalleys(data=data[,i], minDensityThreshold=minDensityThreshold, gridsize=14000, fudge=0.5)
    }

    modpts <- matrix(ncol=g, nrow=1)
    dens <- KernSmooth::bkde(x=data[,i], kernel="normal", bandwidth=bw.select(x=data[,i], gridsize=14000), gridsize=14000)
    pro <- dens$y[peaks$peaks]/sum(dens$y[peaks$peaks])
    dof <- rep(0, g)
    sigma <- array(dim=c(1,1,g))
    delta <- matrix(nrow=1, ncol=g)
    for(t in 1:g)
    {
      modpts[,t] <- peaks$dens$x[peaks$peaks[t]]
      sigma[t] <- mad(x=data[which(init.clust[,i] == t),i], center=modpts[,t])
      delta[,t] <- psych::skew(x=data[which(init.clust[,i] == t),i])
    }
    if(g == 1)
    {
      i.clust = NULL
      warning("Could not get parameters for multicomponent modelling of", i, "\n")
    }
    if(g != 1)
    {
      i.clust = init.clust[,i]
    }

    ## Fit the model
    if(verbose == TRUE)
    {
      cat("Fitting ", g ," component model to ", i,"\n")
    }
    error <- 1
    restarts <- 0

     while(error == 1 & restarts <= max.restarts)
    {
      models[[i]] <- EMMIXskew::EmSkew(dat=data[,i],g=g, distr="mst", clust=i.clust, itmax=itmax, epsilon=epsilon, debug=F,
                                       init=list(pro=pro, modpts=modpts, mu=modpts, dof=dof, sigma=sigma, delta=delta), nrandom = 2)
      pro <- models[[i]]$pro
      modpts <- models[[i]]$pro
      dof <- models[[i]]$dof
      sigma <- models[[i]]$sigma
      delta <- models[[i]]$delta
      if(models[[i]]$error == 1 & restarts <= max.restarts)
      {
        itmax <- itmax + 2000
        warning("Failed to converge, increasing max allowed iterations and restarting", immediate.=TRUE)
      }
    }
    if(verbose == TRUE)
    {
      if(models[[i]]$error == 0)
      {
        cat("Success!\n")
      }
      if(models[[i]]$error == 1)
      {
        cat("Failed to converge within", itmax,"iterations\n")
      }
      if(models[[i]]$error == 2)
      {
        cat("Failed to get initial values\n")
      }
      if(models[[i]]$error == 3)
      {
        cat("Singularity, prepare for the robotic apocalypse\n")
      }
    }
    y <- vector(mode="list", length=length(models[[i]]$pro))
    for(j in models$pro)
    {
      y[[j]]$model <- rdmst(n=table(models[[i]]$clust)[j],p=1, mean=models[[i]]$mu[,j], cov=models[[i]]$sigma[[j]], del=models[[i]]$delta[j])
      y[[j]]$density <- ddmst(dat=sort(y[[j]]$model[,1]), n=length(y[[j]]$model[,1]), p=1, mean=models[[i]][,j], del=models[[i]]$delta[j], cov=models[[i]]$sigma[[j]])*models[[i]]$pro[j]
    }
    models[[i]]$sim <- y

    if(verbose == TRUE)
    {
      cat("Calculating phenotype of cells with posterior probability >", prob,"\n")
    }
    if(get.start == TRUE)
    {
      pheno[[i]]$init.pro <- pro
      pheno[[i]]$init.modpts <- modpts
      pheno[[i]]$init.clust <- init.clus[,i]
      pheno[[i]]$init.dof <- dof
      pheno[[i]]$init.sigma <- sigma
      pheno[[i]]$init.delta <- delta
    }
    pheno <- vector(length=nrow(data))
    # Phenotyping Step
    ## Find phenotypes in 1 component models
    if(g == 1)
    {
      pheno[which(models[[i]]$tau[,1] > prob)] <- paste(i, "-", sep="")
      pheno[which(models[[i]]$tau[,1] < prob)] <- paste(i, "+", sep="")
    }
    ## Find phenotypes in multicomponent models
    if(g > 1)
    {
      for(l in 1:ncol(models[[i]]$tau))
      {
        if(l == 1)
        {
          pheno[which(models[[i]]$tau[,l] > prob)] <- paste(i, "-",sep="")
        }
        if(l > 1)
        {
          pheno[which(models[[i]]$tau[,l] > prob & data[,i] > min(models[[i]]$modpts))] <- paste(i, paste(rep("+",l-1), collapse=""),sep="")
        }
      }

    }
    models[[i]]$pheno <- pheno
    models[[i]]$pheno[which(models[[i]]$pheno == FALSE)] <- "outlier"
    if(verbose == TRUE)
    {
      for(z in unique(models[[i]]$pheno))
      {
        cat("Found", table(models[[i]]$pheno == z)[2], z, "cells\n")
      }
    }
    if(verbose == TRUE)
    {
      cat("Finished\n")
    }
  }
  return(models)
}
