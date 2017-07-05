clust.diff <- function(data, markers, cluster.id,test.method="Chisq", p.adj=TRUE, padj.method="BH")
{
  #' @author Julian Spagnuolo
  #' @param data data.frame of observations, column names must match names in markers
  #' @param markers character vector, length must match number of columns in data
  #' @param cluster.id numeric vector of length matching number of rows in data.
  #' @param test.method Character vector of length 1. Type of test to use in anova step. Default is "Chisq", see \link{\code{help("anova")}} for details.
  #' @param p.adj Logical. Whether to apply multiple testing correction or not. Default is TRUE.
  #' @param padj.method Character vector of length 1. Method of multiple testing correction to apply. Default is "BH" for Benjamini-Hochberg/FDR method, see \code{\link{help("p.adjust")}} for other method.
  #'
  #' @seealso \link{\code{help("p.adjust")}}
  #'
  #'
  #'
  #'
  #'

  # Perform safety checks on data
  if(ncol(data) != length(markers))
  {
    cat("Number of markers does not match the number of columns in data!\n")
    break()
  }
  if(colnames(data) != markers)
  {
    cat("Column names of data do not match markers!\n")
    break()
  }
  if(length(cluster.id) != nrow(data))
  {
    cat("Length of cluster.ids does not match the number of rows in data!\n")
    break()
  }

  # create 3D array for results of dimensions length(clusters)*length(clusters)*length(markers)
  tests <- array(dim = c(length(cluster.id),length(cluster.id),length(markers)),
                 dimnames = list(1:length(cluster.id),1:length(cluster.id),markers))

  # add the cluster.id's to the data.frame for testing
  data$cluster.id <- cluster.id

  # perform tests between each set of cluster id's.
  for(i in 1:length(cluster.id))
  {
    for(j in 1:length(cluster.id))
    {
      # prevents self-cluster testing
      if(j != i)
      {
        # create gaussian glm models for each marker in the set of cluster.id's and test by anova
        # this could be changed to manova or mancova??
        for(m in markers)
        {
          # model for the alt-hypothesis
          x <- glm(formula=data[which(data$cluster.id %in% c(i,j)),m]~data[which(data$cluster.id %in% c(i,j)), "cluster.id"], family=gaussian())
          # model for the null-hypothesis
          x0 <- glm(formula=data[which(data$cluster.id %in% c(i,j)),m]~1, family=gaussian())
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
                   dim = dim(tests),dimnames = list(1:length(cluster.id),1:length(cluster.id),markers))
  }
  return(tests)
}
