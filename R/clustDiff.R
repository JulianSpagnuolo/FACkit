clust.diff <- function(data, markers, cluster.id,test.method="Chisq", p.adj=TRUE, padj.method="BH")
{
  #' @title Cluster Diff
  #' @author Julian Spagnuolo
  #' @param data data.frame of observations, column names must match names in markers
  #' @param markers character vector, length must match number of columns in data
  #' @param cluster.id character vector of length 1 matching the name of the column in data corresponding to the numeric cluster.id vector.
  #' @param test.method Character vector of length 1. Type of test to use in anova step. Default is "Chisq", see \link{\code{help("anova")}} for details.
  #' @param p.adj Logical. Whether to apply multiple testing correction or not. Default is TRUE.
  #' @param padj.method Character vector of length 1. Method of multiple testing correction to apply. Default is "BH" for Benjamini-Hochberg/FDR method, see \code{\link{help("p.adjust")}} for other method.
  #'
  #'
  #'
  #' @export

  # create 3D array for results of dimensions length(clusters)*length(clusters)*length(markers)
  tests <- array(dim = c(length(data[,cluster.id]),length(data[,cluster.id]),length(markers)),
                 dimnames = list(1:length(data[,cluster.id]),1:length(data[,cluster.id]),markers))

  # perform tests between each set of cluster id's.
  for(i in 1:length(data[,cluster.id]))
  {
    for(j in 1:length(data[,cluster.id]))
    {
      # prevents self-cluster testing
      if(j != i)
      {
        # create gaussian glm models for each marker in the set of cluster.id's and test by anova
        # this could be changed to manova or mancova??
        for(m in markers)
        {
          # model for the alt-hypothesis
          x <- glm(formula=data[which(data[,cluster.id] %in% c(i,j)),m]~data[which(data[,cluster.id] %in% c(i,j)), cluster.id], family=gaussian())
          # model for the null-hypothesis
          x0 <- glm(formula=data[which(data[,cluster.id] %in% c(i,j)),m]~1, family=gaussian())
          # test against the null hypothesis.
          # return only the p.value of the test.
          tests[i,j,m] <- anova(x0,x, test=test.method)[2,5]
        }
      }
    }
  }
  # apply multiple testing correction
  if(p.adj == TRUE)
  {
    tests <- array(data=p.adjust(tests,method=padj.method),
                   dim = dim(tests),dimnames = list(1:length(data[,cluster.id]),1:length(data[,cluster.id]),markers))
  }
  return(tests)
}
