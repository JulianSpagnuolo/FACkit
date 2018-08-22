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
  #'
  #' @export
  #'
  #'
  #'

  ## TODO check split clusters for identity - if this makes identical or super close clusters, implement a cluster merging step.

  # First find all clusters not passing minimum size and classify as noise.
  clust.freqs <- as.data.frame(table(expdata[,clust.col]), stringsAsFactors=FALSE)
  clust.freqs <- subset(clust.freqs, Freq < minpts)
  for(i in 1:nrow(clust.freqs)){
    expdata[which(expdata[,clust.col] == clust.freqs[i,"Var1"]), clust.col] <- noise.clust.id
  }
  rm(clust.freqs) # keep the env small.

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
  for(n in unique(clust.ids[which(clust.ids$id != noise.clust.id),"id"])) {
    clusts <- FACkit:::binclust2(binmat = bin.mat[which(clust.ids$id == n),], rowids = row.ids[which(clust.ids$id == n)])

    ## aggregate all clusters not passing min points threshold
    noise <- unlist(clusts[which(lengths(clusts) <= minpts)])
    ## reset the cluster.id for identified noise points to "0"
    clust.ids[which(rownames(clust.ids) %in% noise),1] <- noise.clust.id

    ## separate clusts passing min points threshold
    clusts <- clusts[which(lengths(clusts) > minpts)]
    if(length(clusts) > 1){
      names(clusts) <- paste(n, as.character(seq.int(from=1, to=length(clusts), by=1)), sep = ".")
      ## rename cluster ids
      for(j in 1:length(clusts)) {
        clust.ids[which(rownames(clust.ids) %in% clusts[[j]]),1] <- names(clusts[j])
      }
    }else{
      clust.ids[which(rownames(clust.ids) %in% unlist(clusts)),1] <- n
    }

  }

  cat("Identifying Outliers\n")
  # Mahalanobis Outlier Detection
  clust.meds <- clust.medians(x=cbind(expdata[,markers], clust.ids), markers=markers, clust.col="id", noise.clust.id="0")
  m.dists <- data.frame(dist=vector(), clust.id=vector())
  for(n in unique(clust.ids$id))
  {
    if(n != noise.clust.id){
      cov.mat <- med.cov(expdata = expdata[which(clust.ids$id == n),markers], markers = markers, use.median = TRUE)

      x <- mahalanobis(x=as.matrix(expdata[which(clust.ids$id == n),markers]), center=clust.meds[n,markers],
                       cov = cov.mat, inverted = TRUE, tol=1e-22)

      x <- data.frame(dist = x, clust.id = n)
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

  cat("Reclustering Noise Points To Closest Cluster\n")
  clust.meds <- clust.medians(x=cbind(expdata[,markers], clust.ids), markers=markers, clust.col="id", noise.clust.id="0")
  m.dists <- data.frame(dist=vector(), clust.id=vector(), noise.pt=vector(), stringsAsFactors = FALSE)
  noise <- subset(clust.ids, id == noise.clust.id)

  # Calc mahalanobis dist of noise point to each cluster
  ## TODO this is a bottleneck - can be quite slow... fix by either writing maha function in cpp or trying mvnfast::maha (latter did not work for main part)
  for(n in unique(clust.ids[which(clust.ids$id != noise.clust.id),"id"])){
    cov.mat <- med.cov(expdata = expdata[which(clust.ids$id == n), markers], markers = markers, use.median = TRUE)

    x <- mahalanobis(x = expdata[which(rownames(expdata) %in% rownames(noise)), markers],
                     center = clust.meds[n, markers], cov = cov.mat, inverted = TRUE, tol=1e-22)

    x <- data.frame(dist=x, clust.id=n, noise.pt=rownames(noise), stringsAsFactors = FALSE)
    m.dists <- rbind(m.dists, x)
  }
  # filter results by the distance threshold.
  m.dists <- m.dists[which(m.dists$dist < dist.thresh),]
  # Reclassify the clust ids
  noise.pts <- unique(m.dists$noise.pt)
  for(n in 1:length(noise.pts)){
    current.point <- subset(m.dists, noise.pt == noise.pts[n])
    clust.ids[which(rownames(clust.ids) == noise.pts[n]),"id"] <- current.point[which(current.point$dist == min(current.point$dist)), "clust.id"]
    #clust.ids[which(rownames(clust.ids) == noise.pts[n]),"id"] <- m.dists[which(m.dists$dist == min(m.dists[which(m.dists$noise.pt == noise.pts[n]),"dist"])),"clust.id"]
  }
  print(table(clust.ids$id == "0"))

  cat("Starting the Iterator\n")
  # Reclustering Iterator
  for(i in 2:maxit){
    cat("Identifying Outliers\n")
    # Mahalanobis Outlier Detection
    clust.meds <- clust.medians(x=cbind(expdata[,markers], clust.ids), markers=markers, clust.col="id", noise.clust.id="0")
    m.dists <- data.frame(dist=vector(), super.clust=vector())
    for(n in unique(clust.ids$id))
    {
      if(n != noise.clust.id){
        cov.mat <- med.cov(expdata = expdata[which(clust.ids$id == n), markers], markers = markers, use.median = TRUE)

        x <- mahalanobis(x=as.matrix(expdata[which(clust.ids$id == n),markers]), center=clust.meds[n, markers],
                         cov = cov.mat, inverted = TRUE, tol=1e-22)

        x <- data.frame(dist = x, clust.id = n)
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

    cat("Reclustering Noise Points To Closest Cluster\n")
    clust.meds <- clust.medians(x=cbind(expdata[,markers], clust.ids), markers=markers, clust.col="id", noise.clust.id="0")
    m.dists <- data.frame(dist=vector(), clust.id=vector(), noise.pt=vector(), stringsAsFactors = FALSE)
    noise <- subset(clust.ids, id == noise.clust.id)

    # Calc mahalanobis dist of noise point to each cluster
    ## TODO this is a bottleneck - can be quite slow... fix by either writing maha function in cpp or trying mvnfast::maha (latter did not work for main part)
    for(n in unique(clust.ids[which(clust.ids$id != noise.clust.id),"id"])){
      cov.mat <- med.cov(expdata = expdata[which(clust.ids$id == n), markers], markers = markers, use.median = TRUE)

      x <- mahalanobis(x = expdata[which(rownames(expdata) %in% rownames(noise)), markers],
                       center = clust.meds[n, markers], cov = cov.mat, inverted = TRUE, tol=1e-22)

      x <- data.frame(dist=x, clust.id=n, noise.pt=rownames(noise), stringsAsFactors = FALSE)
      m.dists <- rbind(m.dists, x)
    }
    # filter results by the distance threshold.
    m.dists <- m.dists[which(m.dists$dist < dist.thresh),]
    # Reclassify the clust ids
    noise.pts <- unique(m.dists$noise.pt)
    for(n in 1:length(noise.pts)){
      current.point <- subset(m.dists, noise.pt == noise.pts[n])
      clust.ids[which(rownames(clust.ids) == noise.pts[n]),"id"] <- current.point[which(current.point$dist == min(current.point$dist)), "clust.id"]
      #clust.ids[which(rownames(clust.ids) == noise.pts[n]),"id"] <- m.dists[which(m.dists$dist == min(m.dists[which(m.dists$noise.pt == noise.pts[n]),"dist"])),"clust.id"]
    }
    print(table(clust.ids$id == noise.clust.id))
  }

  return(clust.ids)
}
