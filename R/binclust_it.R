binclust.it <- function(expdata, markers, clust.col, noise.clust.id = "0", minpts = 4, alpha = 0.001, maxit = 5){
  #' @title Bin Clustering
  #' @author Julian Spagnuolo
  #'
  #' @param x data.frame. Contains marker expression and cluster id data.
  #' @param markers Character vector, identifies vectors in data containing the marker expression data.
  #' @param clust.col Character vector, column name identifying the vector of cluster ids.
  #' @param bin.list list of matrices. Each member of the list is a matrix defining the sub-clusters within each original clusters that should be split. Output from binmat function.
  #' @param noise.clust.id If a noise cluster is present in the data, this parameter will remove it from sub-clustering. Ie. DBScan identifies noise points as cluster "0". If NULL, this parameter is not used. Default is NULL.
  #' @param minpts integer. Minimum number of cells a cluster must contain in order to be kept. Default is 4.
  #' @param alpha numeric. alpha param for qchisq call to identify cluster outliers by mahalanobis distance (see qchisq)
  #' @param maxit integer. Maximum number of iterations to run. Default is 5.
  #'
  #' @export
  #'
  #'
  #'

  cat("Creating Binary Matrix\n")
  bin.mat <- sign(expdata[,markers])

  bin.mat[bin.mat == -1] <- 0
  row.ids <- rownames(bin.mat)

  bin.mat <- apply(bin.mat, MARGIN = 2, function(x){as.integer(x)})
  rownames(bin.mat) <- row.ids

  clust.ids <- data.frame(id=as.character(expdata[,clust.col]), stringsAsFactors = FALSE)
  rownames(clust.ids) <- row.ids

  dist.thresh <- qchisq(p=1-alpha, df = length(markers) - 1)

  cat("Running Initial Reclustering\n")
  # Reclustering loop
  for(n in unique(clust.ids$id)) {
    clusts <- FACkit:::binclust2(binmat = bin.mat[which(clust.ids$id == n),], rowids = row.ids[which(clust.ids$id == n)])

    ## aggregate all clusters not passing min points threshold
    noise <- unlist(clusts[which(lengths(clusts) <= minpts)])

    ## separate clusts passing min points threshold
    clusts <- clusts[which(lengths(clusts) > minpts)]
    names(clusts) <- paste(n, as.character(seq.int(from=1, to=length(clusts), by=1)), sep = ".")

    ## reset the cluster.id for identified noise points to "0"
    clust.ids[which(rownames(clust.ids) %in% noise),1] <- noise.clust.id

    ## rename cluster ids
    for(j in 1:length(clusts)) {
      clust.ids[which(rownames(clust.ids) %in% clusts[[j]]),1] <- names(clusts[j])
    }
  }

  cat("Identifying Outliers\n")
  # Mahalanobis Outlier Detection
  clust.meds <- clust.medians(x=cbind(expdata[,markers], clust.ids), markers=markers, clust.col="id", noise.clust.id="0")
  m.dists <- data.frame(dist=vector(), super.clust=vector())
  for(n in unique(clust.ids$id))
  {
    if(n != noise.clust.id){
      cov.mat <- cov(expdata[which(clust.ids$id == n),markers], use="pairwise")

      x <- mahalanobis(x = expdata[which(clust.ids$id == n),markers], center = clust.meds[n,markers], cov = cov.mat,
                       tol=1e-22, inverted = T)

      x <- data.frame(dist=x)
      x$clust.id <- n
      m.dists <- rbind(m.dists, x)
    }
  }

  cat("Reclassifying Outliers and Small Clusters\n")
  # reclassify all points not passing mahalanobis threshold as noise
  clust.ids[which(rownames(clust.ids) %in% rownames(m.dists[which(m.dists$dist >= dist.thresh),])),1] <- noise.clust.id
  # reclassify clusters which are now too small to pass min points threshold
  clust.freqs <- as.data.frame(table(clust.ids$id), stringsAsFactors = F)
  clust.freqs <- clust.freqs[which(clust.freqs$Freq <= minpts),]
  clust.ids[which(clust.ids$id %in% clust.freqs$Var1),1] <- noise.clust.id

  print(table(m.dists$dist > dist.thresh))
  print(table(clust.ids$id == "0"))

  cat("Starting the Iterator\n")
  # Reclustering Iterator
  for(i in 2:maxit){
    cat("Reclustering\n")
    clusts <- FACkit:::binclust2(binmat = bin.mat[which(clust.ids$id == noise.clust.id),], rowids = row.ids[which(clust.ids$id == noise.clust.id)])

    ## aggregate all clusters not passing min points threshold
    noise <- unlist(clusts[which(lengths(clusts) <= minpts)])

    ## separate clusts passing min points threshold
    clusts <- clusts[which(lengths(clusts) > minpts)]

    new.clust.id <- stringr::str_split(string = as.character(max(as.numeric(clust.ids$id))), pattern = "\\.", simplify = T)[,1]
    names(clusts) <- as.character(seq.int(from=as.numeric(new.clust.id) + 1, to=length(clusts)+as.numeric(new.clust.id), by=1))

    ## reset the cluster.id for identified noise points to "0"
    clust.ids[which(rownames(clust.ids) %in% noise),1] <- noise.clust.id

    ## rename cluster ids
    for(j in 1:length(clusts)) {
      clust.ids[which(rownames(clust.ids) %in% clusts[[j]]),1] <- names(clusts[j])
    }

    cat("Identifying Outliers\n")
    # Mahalanobis Outlier Detection
    ## TODO the outlier detection loop is the slowest bottleneck in this function - convert to cpp or find faster r pkg
    clust.meds <- clust.medians(x=cbind(expdata[,markers], clust.ids), markers=markers, clust.col="id", noise.clust.id="0")
    m.dists <- data.frame(dist=vector(), super.clust=vector())
    for(n in unique(clust.ids$id))
    {
      if(n != noise.clust.id){
        cov.mat <- cov(expdata[which(clust.ids$id == n),markers], use="pairwise")

        x <- mahalanobis(x = expdata[which(clust.ids$id == n),markers], center = clust.meds[n,markers], cov = cov.mat,
                         tol=1e-22, inverted = T)

        x <- data.frame(dist=x)
        x$clust.id <- n
        m.dists <- rbind(m.dists, x)
      }
    }
    cat("Reclassifying Outliers\n")
    # reclassify all points not passing mahalanobis threshold as noise
    clust.ids[which(rownames(clust.ids) %in% rownames(m.dists[which(m.dists$dist >= dist.thresh),])),1] <- noise.clust.id
    # reclassify clusters which are now too small to pass min points threshold
    clust.freqs <- as.data.frame(table(clust.ids$id), stringsAsFactors = F)
    clust.freqs <- clust.freqs[which(clust.freqs$Freq <= minpts),]
    clust.ids[which(clust.ids$id %in% clust.freqs$Var1),1] <- noise.clust.id

    print(table(m.dists$dist >= dist.thresh))
    print(table(clust.ids$id == noise.clust.id))
  }

  return(clust.ids)
}
