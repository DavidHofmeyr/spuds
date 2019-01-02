\name{satellite}
\docType{data}
\alias{satellite}
\title{Statlog (landsat satellite) data set.}
\description{
  This dataset contains multi-spectral information taken from 3x3 pixel squares from satellite images. The centre pixel for each case lies in one of seven classes (class 6 is not present in the data set)
1 red soil
2 cotton crop
3 grey soil
4 damp grey soil
5 soil with vegetation stubble
6 mixture class (all types present)
7 very damp grey soil.
}
\usage{satellite}
\format{A list with entries $x (a 6435x36 matrix with each row corresponding to an individual pixel square) and $c (a vector of labels)}
\source{UCI Machine Learning Repository.}
\references{
  Dheeru, D. and E. Karra Taniskidou (2017). UCI machine learning repository. \url{https://archive.ics.uci.edu/ml}
}
\keyword{datasets}
