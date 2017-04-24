net.opt <- function(network, method="mst", start=NULL)
{
  require(igraph)
  if(method == "mst")
  {
    net <- graph.adjacency(as.matrix(network), weighted = TRUE)
    net <- mst(net, weights = E(net)$weight, algorithm="prim")
    return(net)
  }
  else if(method == "trim")
  {
    net <- network - network
    tot <- unique(c(colnames(network), rownames(network)))
    if(is.null(start))
    {
      max.val <- which(network == max(network), arr.ind=TRUE)
      added <- c(rownames(network)[max.val[,"row"]], colnames(network)[max.val[,"col"]])
    }
    else
    {
      added <- c(start)
    }
    excluded <- tot[which(!tot %in% added)]
    if(length(added) > 0)
    {
      if(length(added) == 2)
      {
        net[added[1], added[2]] <- network[added[1], added[2]]
        net[added[2], added[1]] <- network[added[2], added[1]]
      }
      else
      {
        for(i in 1:length(added))
        {
          for(j in 1:length(added))
          {
            net[added[i], added[j]] <- network[added[i], added[j]]
          }
        }
      }
    }
    while(length(excluded) > 0)
    {
      max <- 0
      max.i <- ""
      max.j <- ""
      for(i in 1:length(added))
      {
        for(j in 1:length(excluded))
        {
          lm.i <- added[i]
          lm.j <- excluded[j]
          sum <- 0
          if(lm.i %in% rownames(network) & lm.j %in% colnames(network))
          {
            sum <- sum + network[lm.i,lm.j]
          }
          if(lm.i %in% colnames(network) & lm.i %in% rownames(network))
          {
            sum <- sum + network[lm.j,lm.i]
          }
          if(sum > max)
          {
            max <- sum
            max.i <- lm.i
            max.j <- lm.j
          }
        }
      }
      added <- c(added, max.j)
      excluded <- tot[which(!tot %in% added)]
      net[max.i, max.j] <- network[max.i, max.j]
      net[max.j, max.i] <- network[max.j, max.i]
    }
    opt.net <- graph.adjacency(as.matrix(net), weighted=TRUE)
    #net.trim <- apply(net > 0, c(1,2), sum)
    #opt.net <- graph.adjacency(as.matrix(net.trim), weighted=TRUE)
    return(opt.net)
  }
}