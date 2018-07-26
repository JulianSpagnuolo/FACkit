cutoffs <- function(cut.data, markers, component=NULL)
{
  #' @title cutoffs
  #' @author Julian Spagnuolo
  #' @param cut.data output from cutfind function
  #' @param markers vector of marker names
  #' @param component named vector of the specific component from which the cutoff should be extracted. Names should be identical to markers
  #'
  #' @export

  # Unless specific component is chosen the right-hand cut from the first component is retrieved
  if(is.null(component))
  {
    component <- rep(1, length(markers))
    names(component) <- markers
  }

  cutoffs <- vector(length=length(markers))
  names(cutoffs) <- markers

  for(i in markers)
  {
    if(length(cut.data[[i]]) == 3)
    {
      cutoffs[i] <- cut.data[[i]]$components[[component[i]]]$right
    }
    if(length(cut.data[[i]]) == 4)
    {
      cutoffs[i] <- cut.data[[i]]$right
    }
  }
  return(cutoffs)
}
