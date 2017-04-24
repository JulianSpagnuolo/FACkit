cell.network <- function(data, markers, clusters, phenotypes, state.clusters, major.phenos, halos,
                             distMethod="euclidean", parallel=FALSE, cores=NULL)
{
  require(parallel)
  
    # set the data structure
  data$clusters <- clusters
  data$pheno <- phenotypes
  data$halos <- halos
  
    # Get the mode (peak) fluorescence intensity for each marker in the major phenotype for each landmark cluster
  modes <- matrix(nrow=length(state.clusters), ncol=length(markers), dimnames=list(state.clusters, markers))
    ### Make parallel, switch to using apply...
  for(i in 1:length(state.clusters))
  {
    for(n in markers)
    {
      modes[i,n] <- dmode(x=data[which(data$clusters == state.clusters[i] & data$pheno == major.phenos[i] & data$halos == FALSE),n], gridsize=14000)
    }
  }
    # convert the matrix to dataframe
  modes <- as.data.frame(modes)

    # Create neighbor network
  if(parallel == TRUE)
  {
    # use all but one cores unless cores is explicitly specified (so this can be used in a server queue)
    if(is.null(cores))
    {
      cores <- detectCores() - 1
    }
    cl <- makeCluster(cores)
    clusterExport(cl=cl, c("data"))
    clusterEvalQ(cl=cl, source(file=paste(getwd(),"/neighbor.R",sep="")))
    
    ### neighbor is an internal Mpath function calculating a distance matrix on what is called modes here
    ### neighbor needs to be rewritten to use the Dist function in amap.
    nb <- parApply(t(data[,markers]), 2, function(x) neighbor(x, modes, distMethod=distMethod), cl=cl)
    stopCluster(cl)
  }
  if(parallel == FALSE)
  {
    nb <- apply(t(data[,markers]), 2, function(x) neighbor(x, modes, distMethod=distMethod))
  }
  row.names(nb) <- paste("neighbor", 1:nrow(nb), sep="")
  nb12 <- table(nb["neighbor1",], nb["neighbor2",])
  absent <- rownames(modes)[!rownames(modes) %in% colnames(nb12)]
  if(length(absent) > 0)
  {
    temp <- matrix(0, nrow=nrow(nb12),ncol=length(absent))
    colnames(temp) <- absent
    rownames(temp) <- rownames(nb12)
    nb12 <- cbind(nb12,temp)
  }
  absent <- rownames(modes)[!rownames(modes) %in% rownames(nb12)]
  if(length(absent) > 0)
  {
    temp <- matrix(0, ncol=ncol(nb12), nrow=length(absent))
    rownames(temp) <- absent
    colnames(temp) <- colnames(nb12)
    nb12 <- rbind(nb12, temp)
  }
  nb12 <- nb12[rownames(modes),rownames(modes)]
  return(nb12)
}