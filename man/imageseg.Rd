\name{imageseg}
\docType{data}
\alias{imageseg}
\title{Image segmentation data set.}
\description{
  Identification of seven classes of image regions from outdoor scenes
1 = brickface,
2 = sky,
3 = foliage,
4 = cement,
5 = window,
6 = path,
7 = grass.
Instances are statistics from 3x3 RGB pixel squares.
}
\usage{imageseg}
\format{A list with entries $x (a 2310x19 matrix with each row corresponding to a pixel square) and $c (a vector of labels)}
\source{UCI Machine Learning Repository.}
\references{
  Dheeru, D. and E. Karra Taniskidou (2017). UCI machine learning repository. \url{https://archive.ics.uci.edu/ml}
}
\keyword{datasets}
