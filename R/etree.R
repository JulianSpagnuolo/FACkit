evtree <- function(df, markers, divisionThreshold = 120, divisionFactor = 4, eta0 = 0.3, sigma0 = 0.8, tau1 = 4, tau2 = 2, bddecay = 0.8, kmeanscount = 2)
#' @author Julian Spagnuolo
#' @title Evolving Tree SOM
#' @description Wrapper for the C++ implementation of the Evolving Tree SOM algorithm by Jussi Pakkanen et al (2004); Neural Processing Letters 20: 199–211.
#'
#'
#' @param df Dataset for training the evolving tree SOM
#' @param divisionThreshold
#' @param divisionFactor
#' @param eta0
#' @param sigma0
#' @param tau1
#' @param tau2
#' @param bddecay
#' @param kmeanscount
#'
#'
{
  out <- .C("etreetrain", data = as.double(as.matrix(df[,markers])),
            as.integer(divisionThreshold), as.integer(divisionFactor), as.double(eta0), as.double(sigma0),
            as.double(tau1), as.double(tau2),  as.double(bddecay), as.integer(kmeanscount),
            PACKAGE = "FACkit")
  return(out)
}




