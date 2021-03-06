% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/binclust.R
\name{binclust}
\alias{binclust}
\title{Bin Clustering}
\usage{
binclust(data, markers, nbins = 100, distmet = "eucl",
  percentiles = FALSE)
}
\arguments{
\item{data}{matrix or data.frame containing the data to be clustered}

\item{markers}{character vector of the columns names of data that should be used for clustering}

\item{nbins}{integer. the number of bins to be created for each dimension in data}

\item{distmet}{character vector. the distance metric to be used for matching the closest bin for each marker}

\item{percentiles}{logical. Whether to use percentiles to create the bins for each marker. If FALSE, bins will be evenly spaced across the range for each dimension.}
}
\description{
Clustering by nearest bin
}
\author{
Julian Spagnuolo
}
