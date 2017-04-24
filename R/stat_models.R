state.model <- function(data, markers, clust.id, phenos, halos, major.phenos, states)
{
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
  return(modes)
}