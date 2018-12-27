\name{spuds}
\alias{spuds}
\title{Spectral Partitioning Using Density Separation}
\description{
  Spectral clustering algorithm which selects the number of clusters based on validation using low density separation.
}
\usage{
spuds(X, c0, scale, sigmult, cplus, cmax, lam, gam, intr)
}
\arguments{
  \item{X}{a numeric matrix (num_data x num_dimensions); the dataset to be clustered.}
  \item{c0}{the initial number of clusters.}
  \item{scale}{(optional) either a numeric value for the scale parameter or a function taking as input the data matrix and which returns a scalar.
Default is the square root of the average of the first k
eigenvalues (evals) of the covariance, multiplied by (4/(2+k)/n)^(1/(k+4)).
k is an estimate of the intrinsic dimensionality.}
  \item{sigmult}{(optional) multiplier for scale value. Default is 1.2}
  \item{cplus}{(optional) increment in cluster number with each iteration. Default is 1. If the number of clusters may be large or the data contain numerous outliers then setting this to a larger value (e.g. 10) should accelerate the runtime. In this instance also increase cmax.}
  \item{cmax}{(optional) maximum number of clusters (including outlier clusters). Default is 50}
  \item{lam}{(optional) Density threshold. default is 1}
  \item{gam}{(optional) minimum non-outlier cluster size. default is n/50}
  \item{intr}{(optional) how to determine the intrinsic dimensionality. For values in (0, 1] chooses dimension which accounts for the corresponding proportion of total variability in the data. For values greater than 1 chooses those dimensions which are greater than the largest eigenvalue divided by intr. Two type string options are also implemented. If intr is set to "kaiser" then dimension is the number of eigenvalues >= 1 (columns are standardised first). If intr is set to "elbow" then a simple elbow rule is used.}
}
\value{
  a vector of cluster labels
}
\examples{
### generate data set using the function provided

X <- spuds_datagen(3000, 10, 10)

### produce clustering solution using SPUDS algorithm

sol <- spuds(X$x)

### assess the quality of the solution

cluster_performance(sol, X$c)
}
\keyword{file}

