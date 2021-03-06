% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmode.R
\name{dmode}
\alias{dmode}
\title{dmode}
\usage{
dmode(x, gridsize, fudge, bw.method = c("bw.select", "dpik"), ...)
}
\arguments{
\item{x}{vector or column of data frame in which to find the value at which the density of the distribution is highest}

\item{gridsize}{gridsize}

\item{fudge}{a cofactor of the bandwidth, only used in cases where distributions need to be more smooth than returned by default}

\item{bw.method}{default is "bw.select". This chooses the bandwidth selection method. The only other choose is "dpik". See details}

\item{...}{other parameters to be passed to the bandwidth selection function}
}
\description{
dmode
}
\author{
Julian Spagnuolo
}
