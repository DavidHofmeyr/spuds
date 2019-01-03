
#### compute Euclidean norm of vector

norm_vec <- function(v) sqrt(sum(v^2))

##### computes (squared) pairwise distances

squaredists <- function(X){
  Xn <- rowSums(X^2)
  C <- -2*X%*%t(X)
  C <- C+Xn
  C <- t(t(C)+Xn)
  C[which(C<0)] <- 0
  C
}

##### kmeans++ for clustering embedded data

kmeanspp <- function(X, k, nstart){
  n <- nrow(X)
  obj_opt <- Inf
  sol_opt <- list()
  mn <- colMeans(X)
  C0 <- which.min(distmat(mn, X))
  ds0 <- distmat(X[C0,], X)
  for(it in 1:nstart){
    #C <- sample(1:n, 1)
    #ds <- distmat(X[C,], X)
    C <- C0
    ds <- ds0
    for(i in 2:k){
      C <- c(C, sample(1:n, 1, prob = ds^2/sum(ds^2)))
      drep <- distmat(X[C[i],], X)
      wrep <- which(drep<ds)
      ds[wrep] <- drep[wrep]
    }
    sol <- kmeans(X, X[C,])
    if(sol$tot.withinss<obj_opt){
      obj_opt <- sol$tot.withinss
      sol_opt <- sol
    }
  }
  sol_opt
}

#### is.density.separated(...): Determine if a cluster is
#### separated from the remainder
#### parameters:
# ds <- matrix of pairwise distances in dataset
# den <- vector of values proportional to density at each datum
# ixs <- indices of cluster of interest
# X <- matrix of data points (one per row)
# sig <- scale parameter
# lam <- density height threshold parameter

is.density.separated <- function(ds, den, ixs, X, sig, lam, sol, gam){

  ## establish complement of cluster
  ixs2 <- (1:nrow(X))[-ixs]

  ## find boundary points of cluster
  B <- ixs[unique(apply(rbind(c(), ds[ixs2,ixs]), 1, which.min))]

  ## search over boundary points for path of relatively high density
  for(ix in B){
    ix2 <- ixs2[which.min(ds[ix,ixs2])]
    if(sum(sol==sol[ix2])>gam) denthresh <- lam*min(max(den[ixs]), max(den[which(sol==sol[ix2])]))
    else denthresh <- Inf
    if(min(den[ix], den[ix2])>=denthresh){
      #ts <- seq(.1, .9, length = 9)
      #x.path <- ts%*%t(X[ix,]) + (1-ts)%*%t(X[ix2,])
      #d.path <- distmat(x.path, X)
      #den.path <- rowSums(exp(-d.path^2/2/sig^2))
      ts <- c(.5, .4, .6, .3, .7, .2, .8, .1, .9)
      den.path <- numeric(length(ts))
      i <- 1
      while(i <= length(ts)){
        den.path[i] <- sum(exp(-distmat(ts[i]*X[ix,]+(1-ts[i])*X[ix2,], X)^2/2/sig^2))
        if(den.path[i] < denthresh) i <- length(ts)+1
        else i <- i + 1
      }
      if(min(den.path)>=denthresh) return(FALSE) ## a path is found
    }
  }
  ## if no path is found return TRUE
  TRUE
}


#### SPUDS(...): Automatic spectral clustering algorithm based on finding
#### clustering solutions in which clusters are density separated
#### parameters:
## X <- (n x d) data matrix with rows corresponding to data points
## c0 <- initial number of clusters. Default is 20.
## scale <- either a numeric value for the scale parameter or a function
##          taking as input the data matrix and which returns a scalar.
##          Default is the square root of the average of the first k
##          eigenvalues (evals) of the covariance, multiplied by (4/(2+k)/n)^(1/(k+4)).
##          k is an estimate of the intrinsic dimensionality.
## sigmult <- multiplier for default scale value. Default is 1.2
## cplus <- increment in cluster number with each iteration. Default is 1. If the
##          number of clusters may be large or the data contain numerous outliers
##          then setting this to a larger value (e.g. 10) should accelerate the runtime
## cmax <- maximum cluster number. Default is 50 (including outllier clusters)
## lam <- density threshold parameter. Default is 1
## gam <- outlier cluster threshold. Only clusters with a greater number of observations
##        are considered substantial. Default is n/50
## intr <- how to determine the intrinsic dimensionality. For values in (0, 1] chooses
##        dimension which accounts for the corresponding proportion of total variability
##        in the data. For values greater than 1 chooses those dimensions which are greater
##        than the largest eigenvalue divided by intr. Two type string options are also
##        implemented but these tend to select far too few dimensions.

spuds <- function(X, c0 = NULL, scale = NULL, sigmult = 1.2, cplus = NULL, cmax = NULL, lam = NULL, gam = NULL, intr = 100){

  X <- as.matrix(X)

  ## rescale columns of X to have unit variance

  V <- cov(X)
  if(min(diag(V))<1e-10){
    X <- X[,-which(diag(V)<1e-10)]
    V <- cov(X)
  }
  if(max(abs(diag(V)-1))>1e-5){
    X <- t(t(X)/sqrt(diag(V)))
    V <- cov(X)
  }

  ## set up parameter values
  n <- nrow(X)
  d <- ncol(X)

  if(is.null(c0)) c0 <- 20
  if(is.null(cmax)) cmax <- 49
  if(is.null(lam)) lam <- 1
  if(is.null(cplus)) cplus <- 1
  if(is.null(gam)) gam <- ceiling(n/50)

  if(is.null(scale)){
    evals <- eigen(V)$values
    if(intr=='elbow') nf <- get_elbow(evals)
    if(is.numeric(intr)){
      if(intr<1) nf <- min(which(cumsum(evals)>(intr*sum(evals))))
      else nf <- sum(evals>(max(evals)/intr))
    }
    if(intr=='kaiser') nf <- sum(evals>=1)
    sig <- sqrt(mean(evals[1:nf]))*(4/(2+nf)/n)^(1/(nf+4))
    if(!is.null(sigmult)) sig <- sig*sigmult
  }
  else if(is.numeric(scale)) sig <- scale
  else if(is.function(scale)) sig <- scale(X)
  else stop('scale must be a numeric or a function of a data matrix')

  ## compute graph Laplacian

  ds <- squaredists(X)
  W <- exp(-ds/2/sig^2)
  diag(W) <- 0
  D <- rowSums(W)
  if(min(D)<(1e-5*max(D))) D[which(D<(1e-5*max(D)))] <- 1e-5*max(D) ### for better numerical stability
  L <- (1/sqrt(D))%*%t(1/sqrt(D))*W

  ## compute first cluster solution
  nclust <- c0
  e <- try(eigs_sym(L, min(cmax+1, nclust*2))$vectors/sqrt(D))
  if(inherits(e, 'try-error')) e <- eigen(L)$vectors/sqrt(D)

  E <- e[,1:nclust]/apply(matrix(e[,1:nclust], nrow = n), 1, norm_vec)
  sol <- kmeanspp(E, nclust, nstart = 10)$cluster

  ## determine if all non-outlier clusters are separated
  separate <- TRUE
  for(i in 1:nclust){
    if(separate){
      if(sum(sol==i)<=gam) p <- TRUE
      else p <- is.density.separated(ds, D+1, which(sol==i), X, sig, lam, sol, gam)
      if(!p) separate <- FALSE
    }
  }

  ## if all non-outlier clusters are separated then increase the number of clusters until this no longer holds

  if(separate && nclust<cmax){
    repeat{
      nclust <- nclust + cplus
      if(ncol(e)<nclust) e <- eigs_sym(L, min(n, nclust+5*cplus))$vectors/sqrt(D)
      E <- e[,1:nclust]/apply(matrix(e[,1:nclust], nrow = n), 1, norm_vec)
      temp <- kmeanspp(E, nclust, nstart = 10)$cluster

      if(nclust>cmax){
        sol <- temp
        break
      }

      separate <- TRUE
      for(i in 1:nclust){
        if(separate){
          if(sum(temp==i)<=gam) p <- TRUE
          else p <- is.density.separated(ds, D+1, which(temp==i), X, sig, lam, temp, gam)
          if(!p) separate <- FALSE
        }
      }

      if(cplus > 1) sol <- temp

      if(separate) sol <- temp
      else break
    }
  }

  ## if in the initial cluster solution there was a non-outlier cluster not separated from the remainder then
  ## nclust==c0 and we decrease the number of clusters until they are all separated. Additionally, if cplus > 1
  ## and we increased the number of clusters from c0 then we need to refine the solution by decreasing the
  ## number of clusters until all are separated.

  if(nclust==c0 || (nclust>c0 && cplus>1)){
    repeat{
      nclust <- nclust - 1
      E <- e[,1:nclust]/apply(matrix(e[,1:nclust], nrow = n), 1, norm_vec)
      sol <- kmeanspp(E, nclust, nstart = 10)$cluster
      separate <- TRUE
      for(i in 1:nclust){
        if(separate){
          if(sum(sol==i)<=gam) p <- TRUE
          else p <- is.density.separated(ds, D+1, which(sol==i), X, sig, lam, sol, gam)
          if(!p) separate <- FALSE
        }
      }
      if(separate || sum(table(sol)>gam)<=2) break
    }
  }

  ## Finally we merge outlier clusters with their nearest cluster until no outlier clusters remain.
  ## For simplicity we do this sequentially, which could result in a slightly
  ## different solution depending on the order in which the outlier clusters are
  ## considered. As the number of data is few, we are not concerned by this; the
  ## final solutions will not differ substantially.

  outliers <- which(table(sol)<=gam)
  while(length(outliers)>0){
    ix1 <- which(sol==outliers[1])
    inliers <- (1:n)[-ix1]
    ixa <- inliers[which(matrix(ds[ix1, inliers], nrow = length(ix1))==min(ds[ix1, inliers])[1], arr.ind = TRUE)[1,2]]
    sol[ix1] <- sol[ixa]
    sol[which(sol>outliers[1])] = sol[which(sol>outliers[1])] - 1
    outliers <- which(table(sol)<=gam)
  }
  ret = numeric(n)
  for(i in 1:length(unique(sol))) ret[which(sol==unique(sol)[i])] = i
  ret
}


### generates data Gaussian mixture and then applies a distortion to the points to create non-convexities in the clusters
## params:
# nits = number of data
# nclust = number of clusters
# dim = number of dimensions

spuds_datagen = function(nits, nclust, dim, scale = NULL, curve = NULL){
  if(is.null(curve)) curve = dim
  if(is.null(scale)) scale = 1
  data = c()
  labels = c()
  clusts = list()
  for(j in 1:1){
    dats = c()
    labs = c()
    covs = matrix(runif(nclust*dim)^5, nrow = nclust)*sqrt(dim)*runif(1)*4
    covs = covs/norm(covs)*dim
    sphere = matrix(rnorm(100*nclust*dim), ncol = dim)
    MEANS = kmeans(sphere, nclust)$centers*sqrt(dim/10)/scale
    ps = runif(nclust)+1
    ps = ps/sum(ps)
    for(i in 1:nclust){
      ni = ceiling(ps[i]*nits)
      clust = matrix(rnorm(ni*dim), ni, dim)
      clust = t(t(clust)*covs[i,]+MEANS[i,])
      clusts[[i]] = clust
      labs = c(labs, rep(i, ni))
    }
    if(curve>0){
      for(i in 1:nclust){
        dmn = distmat(clusts[[i]], MEANS[-i,])
        mins = apply(dmn, 1, which.min)
        for(l in 1:nrow(clusts[[i]])){
          ws = exp(-dmn[l,]*10)
          cent = colSums(ws*MEANS[-i,])/sum(ws)
          clusts[[i]][l,] = cent+(clusts[[i]][l,]-cent)/(sqrt(sum((clusts[[i]][l,]-cent)^2)))*curve
        }
      }
    }
    for(i in 1:nclust) dats = rbind(dats, clusts[[i]])
    smp = sample(1:length(dats[,1]), nits)
    dats = dats[smp,]
    labs = labs[smp]

    data = rbind(data, dats)
    labels = c(labels, labs)
  }

  X = data + rnorm(nrow(data)*ncol(data))/2

  for(i in 1:dim){
    X[,i] = (X[,i]-mean(X[,i]))/sd(X[,i])
  }

  list(x = X, c = labels)
}


get_elbow <- function(yvals){
  d <- length(yvals)
  angles <- sapply(1:(d-1), function(i){
    a1 <- atan((i-1)/(d-1)*(yvals[1]-yvals[d])/(yvals[1]-yvals[i]))
    a2 <- atan((yvals[i]-yvals[d])/(yvals[1]-yvals[d])*(d-1)/(d-i))
    return(a1+a2)
  })
  which.min(angles)
}

#### computes external cluster validity metrics

cluster_performance = function(assigned, labels, beta = 1){
  n <- length(labels)
  T <- table(assigned, labels)
  RS <- rowSums(T)
  CS <- colSums(T)
  
  ## V-measure
  
  CK <- - sum(apply(T, 1, function(x) return(sum(x[which(x>0)]*log(x[which(x>0)]/sum(x))))))/n
  KC <- - sum(apply(T, 2, function(x) return(sum(x[which(x>0)]*log(x[which(x>0)]/sum(x))))))/n
  K <- - sum(apply(T, 1, function(x) return(sum(x)*log(sum(x)/n))))/n
  C <- - sum(apply(T, 2, function(x) return(sum(x)*log(sum(x)/n))))/n
  if(C!=0){
    h <- 1 - CK/C
  }
  else{
    h <- 0
  }
  if(K!=0){
    c <- 1 - KC/K
  }
  else{
    c <- 0
  }
  if(h==0 && c==0) v.measure <- 0
  else v.measure <- (1+beta)*h*c/(beta*h+c)
  
  ## Purity
  
  purity <- sum(apply(T, 1, function(x) return(max(x))))/n
  
  ## Adjusted Rand Index
  
  O <- sum(sapply(T, function(t) choose(t, 2)))
  E <- (sum(sapply(RS, function(t) choose(t, 2)))*sum(sapply(CS, function(t) choose(t, 2))))/choose(n, 2)
  M <- (sum(sapply(RS, function(t) choose(t, 2))) + sum(sapply(CS, function(t) choose(t, 2))))/2
  adj.rand <- (O-E)/(M-E)
  
  
  ## Normalised Mutual Information
  
  prod <- RS%*%t(CS)
  Tp <- T
  Tp[which(T==0)] <- 1e-10
  IXY <- sum(T*log(Tp*n/prod))
  HX <- sum(RS*log(RS/n))
  HY <- sum(CS*log(CS/n))
  NMI <- IXY/sqrt(HX*HY)
  
  c(adj.rand = adj.rand, purity = purity, v.measure = v.measure, nmi = NMI)
}


