\name{spuds_datagen}
\alias{spuds_datagen}
\title{Simulation of Data Sets for Clustering}
\description{
  Generates data from a Gaussian mixture, and then applies a perturbation to each point to render non-convex clusters.
}
\usage{
  spuds_datagen(nits, nclust, dim, scale, curve)
}
\arguments{
  \item{nits}{integer number of data points.}
  \item{nclust}{integer number of clusters.}
  \item{dim}{integer number of dimensions.}
  \item{scale}{(optional) numeric. Larger values result in greater cluster overlap. Default is 1.}
  \item{curve}{(optional) numeric. Larger values result in greater perturbation from Gaussian. Default is number of dimensions.}
}
\value{
  a named list with fields
  $x: the matrix of data
  $c: the vector of class labels
}
\keyword{file}

