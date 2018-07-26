## culling empty sub-clusters, and setting new cluster names.
split.merge <- function(x, markers, clust.col, binlist, noise.clust.id=NULL)
{
  #' @author Julian Spagnuolo
  #' 
  #' @param x data.frame. Contains marker expression and cluster id data.
  #' @param markers Character vector, identifies vectors in data containing the marker expression data.
  #' @param clust.col Character vector, column name identifying the vector of cluster ids.
  #' @param binlist list of matrices. Each member of the list is a matrix defining the sub-clusters within each original clusters that should be split. Output from binmat function.
  #' @param noise.clust.id If a noise cluster is present in the data, this parameter will remove it from sub-clustering. Ie. DBScan identifies noise points as cluster "0". If NULL, this parameter is not used. Default is NULL.
  #'
  #' 
  #' @importFrom FACkit bindist

  
  # step 1: re-cluster data in original clusters amongst the sub-clusters defined by clust.split function (binlist output)
  x$split.clusts <- x[,clust.col]
  
  for(i in 1:length(bin.list))
  {
    node.dist <- bindist(binmat = bin.list[[i]], data=as.matrix(x[which(x[,clust.col] == names(bin.list)[i]),markers]))
    
    keep <- vector()
    for(n in 1:length(node.dist$node.pops))
    {
      if(!is.null(node.dist$node.pops[[n]]))
      {
        keep <- append(x=keep, values=n)
      }
    }
    bin.list[[i]] <- bin.list[[i]][keep,]
    
    # recluster the members if a sub-cluster was removed.
    if(length(keep < nrow(bin.list[[i]])))
    {
      cat("culling ", nrow(bin.list[[i]]) - length(keep), "sub-clusters from cluster", i,"\n")
      node.dist <- bindist(binmat = bin.list[[i]], data=as.matrix(x[which(x[,clust.col] == names(bin.list)[i]),markers]))
    }
    rm(keep)
    
    # sub-cluster member id's are indexed relative to their position in the subset of the data matrix by the bindist function.
    # This step reverses that and renames the cluster as n.1, n.2, etc.
    for(n in 1:length(node.dist$node.pops))
    {
      real.rownames <- rownames(x[which(x[,clust.col] == names(bin.list)[i]),])[node.dist$node.pops[[n]]]
      ##node.dist$node.pops[[n]] <- real.rownames
      x[rownames(x) %in% real.rownames,"split.clusts"] <- paste(names(bin.list)[i],n,sep=".")
    }
  }
  
  
  # Step 2: Define clusters by median marker expression, binarise the data in preparation for merging.
  if(!is.null(noise.clust.id))
  {
    split.clusts <- unique(x[which(x$split.clusts != noise.clust.id),]$split.clusts)
  }
  else{
    split.clusts <- unique(x$split.clusts)
  }
  
  clust.defs <- matrix(nrow=length(split.clusts), ncol=length(markers), dimnames = list(c(split.clusts), c(markers)))
  
  for(i in 1:length(split.clusts))
  {
    for(n in 1:length(markers))
    {
      # Some clusters have only one member, thus the skew-normal location (median) function (snlocation) will not work
      if(nrow(x[which(x$split.clusts == split.clusts[i]),]) > 1)
      {
        clust.defs[split.clusts[i],markers[n]] <- snlocation(markdat = x[which(x$split.clusts == split.clusts[i]),markers[n]])
      }
      else{
        clust.defs[split.clusts[i],markers[n]] <- x[which(x$split.clusts == split.clusts[i]),markers[n]]
      }
    }
  }
  
  # make the binary cluster definitions
  clust.defs <- as.matrix(clust.defs)
  clust.defs[sign(clust.defs) == 1] <- 1
  clust.defs[sign(clust.defs) == -1] <- 0
  
  clust.defs <- matrix(data=as.integer(clust.defs), byrow=F, nrow=nrow(clust.defs), ncol=ncol(clust.defs),dimnames = dimnames(clust.defs))
  
  # Step 3: Super cluster merging
  super.clusts <- unique(clust.defs)
  clust.merge <- vector(mode="list", length=nrow(super.clusts))
  names(clust.merge) <- rownames(super.clusts)
  
  for(i in rownames(clust.defs))
  {
    for(n in rownames(super.clusts))
    {
      if(all(clust.defs[i,] == super.clusts[n,]) == TRUE)
      {
        clust.merge[[n]] <- append(x = clust.merge[[n]], values=i)
      }
    }
  }
  # find the unique binary clusters with more than one member cluster - these are the super clusters.
  super.clusts <- clust.merge[lapply(clust.merge, length) > 1]
  
  # super clusters are named after the first unique cluster i.e. 1.1 and 60 bin clusts are identical and all cells will be labelled as 1.1
  x$super.clusts <- x[,"split.clusts"]
  for(i in 1:length(super.clusts))
  {
    x[which(x$split.clusts %in% super.clusts[[i]]),"super.clusts"] <- names(super.clusts)[i]
  }
  
  return(x)
}