enrichTest <- function(x, clust.col, noise.clust.id=NULL, cat.col, equal.props=FALSE, alternative="greater")
{
  #' @title Enrichment Testing
  #'
  #' @author Julian Spagnuolo
  #'
  #' @param x data.frame/matrix. Data to be tested against, should contain column of cluster IDs and categorical/annotation column.
  #' @param clust.col character vector. Column name in x containing the cluster ID's.
  #' @param noise.clust.id character. ID of noise clusters if present (i.e DBscan ID's noise as "0"). If NULL this will not be used. Default is NULL.
  #' @param cat.col character vector. Column name in x with the categorical (factor) data to be tested for enrichment.
  #' @param equal.props logical/boolean. If TRUE binomial test will use an equal hypothesised probability of success in the binomial test, if FALSE binomial test uses the global proportion (including the noise points if present). Default is FALSE. Note: The hypergeometric test will always use the global proportion.
  #' @param alternative character. Indicates alternative hypothesis of binomial test. Must be one of "two.sided", "greater", or "less". You can specify just the initial letter.Default is "greater"
  #'
  #' @details For the binomial test,
  #'
  #' @export

  # get vector of clusters to test for enrichment.
  unique.clusts <- unique(x[,clust.col])
  # remove the noise points if requested.
  if(!is.null(noise.clust.id))
  {
    unique.clusts <- unique.clusts[!(unique.clusts == noise.clust.id)]
  }

  unique.cats <- unique(x[,cat.col])

  if(isTRUE(equal.props))
  {
    prob <- data.frame(Var1=unique.cats, Freq=1/length(unique.cats))
  }
  else
  {
    prob <- table(x[,cat.col])/nrow(x)
    prob <- as.data.frame(prob) # cat == Var1, prob == Freq
  }

  # initialise results output object
  results <- data.frame(cluster=rep(unique.clusts, times=length(unique.cats)),
                        category=rep(unique.cats, each=length(unique.clusts)), cluster.size=NA,
                        proportion=NA, binomial.pval=NA, binomial.FDR=NA, hg.pval=NA, hg.FDR=NA,
                        stringsAsFactors = TRUE, row.names = NULL)
  results <- results[order(results$cluster),]
  rownames(results) <- 1:nrow(results)

  # begin testing
  for(i in unique.clusts)
  {
    for(n in unique.cats)
    {
      # perform binomial test
      res <- binom.test(x=nrow(x[which(x[,clust.col] == i & x[,cat.col] == n),]),
                        n=nrow(x[which(x[,clust.col] == i),]),
                        p=prob[which(prob$Var1 == n),"Freq"],
                        conf.level = 0.95, alternative = alternative)[c("estimate","p.value")]

      results[which(results[,"cluster"] == i & results[,"category"] == n),c("proportion","binomial.pval")] <- res
      results[which(results[,"cluster"] == i),"cluster.size"] <- nrow(x[which(x[,clust.col] == i),])

      # perform hypergeometric test
      res <- phyper(q=nrow(x[which(x[,clust.col] == i & x[,cat.col] == n),]),
                    m=nrow(x[which(x[,cat.col] == n),]),
                    n=nrow(x)-nrow(x[which(x[,cat.col] == n),]),
                    k=nrow(x[which(x[,clust.col] == i),]),
                    lower.tail=FALSE, log.p=FALSE)
      results[which(results[,"cluster"] == i & results[,"category"] == n),"hg.pval"] <- res
    }
  }
  results$binomial.FDR <- p.adjust(p=results$binomial.pval, method="fdr")
  results$hg.FDR <- p.adjust(p=results$hg.pval, method="fdr")

  results$bin.FDR <- signif(results$binomial.FDR, digits=3)
  results$bin.pval <- signif(results$binomial.pval, digits=3)
  results$prop <- signif(results$proportion, digits=3)

  return(results)
}
