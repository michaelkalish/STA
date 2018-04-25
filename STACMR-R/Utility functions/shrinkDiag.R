shrinkDiag <- function(x, shrink) {
# shrinkDiag <-function(x, shrink)
# x (t*n): t iid observations on n random variables
# sigma (n*n): invertible covariance matrix estimator
#
# Shrinks towards diagonal matrix
# if shrink is specified, then this constant is used for shrinkage
# ****************************************************************
# converted from matlab 12 September 2016
# revised to handle missing data 13 October 2016
# ****************************************************************
 
  if (!is.matrix(x)) {x = as.matrix(x)} 
  if (missing(shrink)) {shrink = -1}
  sigma = var(x); shrinkage=0;
  if (ncol(x) > 1) {
  # de-mean returns
  t = nrow(x)
  meanx = colMeans(x, na.rm=TRUE)
  y = x - t(matrix(rep(meanx,t),ncol(x),t)); y[is.na(y)]=0; y=y^2;
  
  # compute sample covariance matrix
  a=as.matrix(x); a[is.na(x)]=0; a[!is.na(x)]=1; n=t(a)%*%a; a=(n-1)/n
  sample = cov(x,use='pairwise.complete.obs')*a; sample[is.na(sample)]=0
                 
  # compute prior
  if (length(sample)==1) {prior=sample} else {prior=diag(diag(sample))}
               
  if (shrink < 0) { 
    # what we call p 
    y = t(y)%*%y; phiMat= y/t-sample^2
    phi=sum(phiMat)
                                
    # what we call r
    rho=sum(diag(phiMat))
    
    # what we call c
    g = norm((sample-prior), type='F')^2
    
    # compute shrinkage constant
    kappa=(phi-rho)/g
    shrinkage=max(0,min(1,kappa/t))
  } else {shrinkage = shrink}
    
  # compute shrinkage estimator
  sigma=shrinkage*prior+(1-shrinkage)*sample
  }
  
  return (list(sigma=sigma, shrinkage=shrinkage))
}
  