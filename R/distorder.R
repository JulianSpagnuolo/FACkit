cell.order <- function(data, markers, clust.id, phenos, halos, major.phenos, states)
{
  
  # This function models cell states using the mode of Diffusion Components to determine the "central coordinates" of each state
  # Then the distance of each cell in the cluster core of each state to the cell state model (mode of core cell DC's) is determined
  # to define the maximum distance of core cells in the state
  # The distance between each cell state model is also determined
  # For a non-core cell clustered in state 1 to belong in the branch its distance to state 2 must be less than
  # the sum distance between each state AND the max distance of core cells belonging to each cell state.
  # The the most distant cell clustered to state 1 is set as the first cell and the distance of all cells clustered to state 1 and 2
  # that were determined to belong in the branch is determined and used to order the cells in "pseudotime"
  # the order list is then returned
  # At present time this works for 2 cell states that are located next to each other at the end of the graph tree and
  # do not turn around a corner
  
  
  ## To Test:
  ## Whether this will work for >2 states
  ### - here we can use the two most distant states to find the first cell on which to order all others i.e. A->C (skipping B)
  ### - else ordering could be performed in 2 or more steps. ie. first order cells A->B then B->C and check against A->C
  
  ## Whether this will work for states that turn a corner in the graph tree needs to be determined.
  
  ## It may be nice to keep the location of cell state models in the order list to located them in the pseudotime plot
  ## This may interfere with the subsequent modelling step.
  
  data$clust.id <- clust.id
  data$phenos <- phenos
  data$halos <- halos
  
  # Model the cell states
  cat("Finding centres of state clusters\n")
  modes <- matrix(nrow=length(unique(states)), ncol=length(markers), dimnames=list(states,markers))
  
  for(i in 1:length(states))
  {
    for(n in markers)
    {
      modes[i,n] <- dmode(x=data[which(data$clust.id == states[i] & data$phenos == major.phenos[i] & data$halos == FALSE),n], gridsize=14000)
    }
  }
  stat.dist <- dist(modes[states,])
  
  #for(i in 1:nrow(modes))
  #{
  #  rownames(modes)[i] <- paste("state",rownames(modes)[i],sep="")
  #}
  #stat.name <- rownames(modes)
  
  
  # Get max distance of core cells to state model.
  ## Get the distance of cells within the state cluster cores
  cat("finding the max distance in cluster cores and between states\n")
  #dat <- rbind(modes, data[which(data$clust.id == states[1] & data$phenos == major.phenos[1] & data$halos == FALSE),markers])
  
  #stat1.core <- as.matrix(Dist(dat[which(rownames(dat) != stat.name[2]),markers], method="euclidean", nbproc = cores))[,1]
  stat1.core <- eucl(x=modes[states[1],],y=data[which(data$clust.id == states[1] & data$phenos == major.phenos[1] & data$halos == FALSE),markers])
  stat1.core <- max(stat1.core)
  cat("Max distance in state 1 is: ",stat1.core,"\n")
  
  #dat <- rbind(modes, data[which(data$clust.id == states[2] & data$phenos == major.phenos[2] & data$halos == FALSE),markers])
  #stat2.core <- as.matrix(Dist(dat[which(rownames(dat) != stat.name[1]),markers], method="euclidean", nbproc = cores))[,1]
  stat2.core <- eucl(x=modes[states[2],],y=data[which(data$clust.id == states[2] & data$phenos == major.phenos[2] & data$halos == FALSE),markers])
  stat2.core <- max(stat2.core)
  cat("Max distance in state 2 is: ",stat2.core,"\n")
  
  max.dist <- sum(stat.dist, stat1.core, stat2.core)
  cat("Max distance between states is: ",max.dist,"\n")
  
  # Find cells within the branch
  cat("finding the cells in the branch\n")
  #dat <- rbind(modes, data[which(data$clust.id == states[1]),markers])
  
  #stat1.dist <- as.matrix(Dist(dat[which(rownames(dat) != stat.name[2]),markers], method="euclidean", nbproc = cores))
  stat1.dist <- eucl(x=modes[states[1],], y=data[which(data$clust.id ==  states[1]),markers])
  names(stat1.dist) <- rownames(data[which(data$clust.id == states[1]),markers])
  stat1.dist <- stat1.dist + stat2.core
  
  #dat <- rbind(modes, data[which(data$clust.id == states[2]),markers])
  #stat2.dist <- as.matrix(Dist(dat[which(rownames(dat) != stat.name[1]),markers], method="euclidean", nbproc = cores))
  stat2.dist <- eucl(x=modes[states[2],], y=data[which(data$clust.id ==  states[2]),markers])
  names(stat2.dist) <- rownames(data[which(data$clust.id == states[2]),markers])
  stat2.dist <- stat2.dist + stat1.core
  
  stat1.branch <- names(stat1.dist[which(stat1.dist  <= sum(stat.dist, stat1.core, stat2.core))])
  stat2.branch <- names(stat2.dist[which(stat2.dist  <= sum(stat.dist, stat1.core, stat2.core))])
  
  # Find the first cell
  cat("finding cell 1\n")
  dat <- data[c(stat1.branch, stat2.branch),markers]
  dat <- rbind(modes, dat)
  
  #cell1 <- as.matrix(Dist(dat[which(rownames(dat) != stat.name[2]),], method="euclidean", nbproc = cores))
  cell1 <- eucl(x=modes[states[1],], y=data[c(stat1.branch, stat2.branch),markers])
  names(cell1) <- c(stat1.branch, stat2.branch)
  cell1 <- names(which(cell1 == stat1.core))
  
  # Order cells in branch according to distance to cell 1
  cat("Getting distance to cell 1\n")
  #dat <- rbind(data[cell1,markers],data[which(rownames(data) != cell1 & rownames(data) %in% c(stat1.branch,stat2.branch)),markers])
  
  #ord <- as.matrix(Dist(dat, method="euclidean", nbproc = cores))
  ord <- eucl(x=data[cell1,markers], y=data[which(rownames(data) != cell1 & rownames(data) %in% c(stat1.branch,stat2.branch)),markers])
  names(ord) <- rownames(data[which(rownames(data) != cell1 & rownames(data) %in% c(stat1.branch,stat2.branch)),markers])
  cat("Sorting by distance to cell 1\n")
  ord <- sort(ord, decreasing=FALSE)
  cat("Finished!\n")
  return(ord)
}