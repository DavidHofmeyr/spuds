\name{spuds-package}
\alias{spuds-package}
\docType{package}
\title{
\packageTitle{spuds}
}
\description{
\packageDescription{spuds}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{spuds}
\packageIndices{spuds}
This package provides an implementation of the SPUDS algorithm introduced in Hofmeyr, D. (2017) "Improving spectral clustering using the asymptotic value of the normalised cut". The main function provided in the package is spuds(). The only mandatory argument is a data matrix, all optional arguments are given the default values used in the paper. Details of optional arguments can be seen using the function documentation. A function to simulate data sets like those used in the paper is also provided (spuds_datagen(ndata, nclust, ndimensions)). Publicly available benchmark data sets are also included. Use data(package = 'spuds') to see a list of included data sets.
}
\references{
Hofmeyr, D. (2017) Improving spectral clustering using the asymptotic value of the normalised cut, \emph{ArXiv preprint}, ArXiv 1703.09975.
}
\keyword{ package }

