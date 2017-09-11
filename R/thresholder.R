thresholder < function(data, markers, thresholds)
  #' @author Julian Spagnuolo
  #' @title Threshold Phenotyping
  #' @param data numeric matrix/data.frame. Contains FACS marker expression data
  #' @param markers character vector. Names of markers in data to phenotype using thresholds
  #' @param thresholds named numeric vector. Expression cut points on which to call positive (1) or negative (0) phenotypes.
{
  results <- matrix(nrow = nrow(data), ncol=length(markers), dimnames = list(rownames(data), markers))

  for(i in 1:nrow(data))
  {
    for(n in markers)
    {
      if(data[i.n] >= threshold[n])
      {
        results[i,n] <- 1
      }
      if(data[i.n] < threshold[n])
      {
        results[i,n] <- 0
      }
    }
  }
  return(results)
}
