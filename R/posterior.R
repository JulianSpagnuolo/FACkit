posterior <- function(model.params, data)
{
  #' @title Posterior Probability
  #'
  #' @author Julian Spagnuolo
  #' @param model.params list, output of phenosampler
  #' @param data numeric matrix or data.frame. ncols must equal length(names(model.params))
  #'
  #'
  #'
  #'

  #initialise output storage object
  results <- list()

  for(i in names(model.params))
  {
    ll <- ddmix(dat=data[,i], n = nrow(data), p = 1, g = 2, distr = "mst",
                mu=matrix(data = c(mean(model.params[[i]][,"mu1"]),mean(model.params[[i]][,"mu2"])), nrow = 1, ncol = 2),
                sigma=array(data=c(mean(model.params[[i]][,"sigma1"]),mean(model.params[[i]][,"sigma2"])), dim = c(1,1,2)),
                dof=c(mean(model.params[[i]][,"dof1"]),mean(model.params[[i]][,"dof2"])),
                delta=matrix(data = c(mean(model.params[[i]][,"delta1"]),mean(model.params[[i]][,"delta2"])), nrow = 1, ncol = 2))

    post1 <- exp(ll[,1]-max(ll[,1]))/sum(exp(ll[,1]-max(ll[,1])))
    post2 <- exp(ll[,2]-max(ll[,2]))/sum(exp(ll[,2]-max(ll[,2])))

    results[[i]]$dens1 <- ll[,1]
    results[[i]]$post1 <- post1/max(post1)
    results[[i]]$dens2 <- ll[,2]
    results[[i]]$post2 <- post2/max(post2)
  }
  return(results)
}
