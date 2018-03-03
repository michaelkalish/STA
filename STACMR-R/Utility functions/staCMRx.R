staCMRx <- function (data, model=NULL, E = list(), shrink=-1, tolerance=0, proc=-1, approx=0) {
  # function [x, f, shrinkage] = staCMRx (data, model, E, shrink, tolerance, proc, approx)
  # Multidimensional CMR
  # data is cell array of data or structured output from staSTATS 
  # model is a nvar*k matrix specifying the linear model; default = ones(nvar,1)
  # E is partial order
  # shrink is parameter to control shrinkage of covariance matrix;
  # 0 = no shrinkage; 1 = diagonal matrix; -1 = calculate optimum
  # tolerance ?
  # proc is number of processors (-1 means all possible)
  # approx ?
  # returns:
  # x = best fitting CMR values to y-means
  # fval = fit statistic
  # shrinkage = estimated shrinkage of covariance matrix
  
  # *************************************************************************
  # modified from matlab 14 September 2016
  # minor bug fixed 13 October 2016
  # updated 8 March 2017
  # *************************************************************************
  #
  # set up defaults
  tol <- 10e-5
  if (missing(E)) {E = list()}
  if (missing(shrink)) {shrink = -1}
  if (missing(tolerance)) {tolerance = 0}
  if (missing(proc)) {proc = -1}
  if (missing(approx)) {approx = 0}
  
  # get stats from data (depending on its form)
  if (is(data,"data.frame")) {
    y = gen2list (data) # convert from general format
    y = staSTATS (y, shrink) # get stats
  } else if (is.null(data[[1]]$means)) {y = staSTATS(data, shrink) # in list form, get stats
  } else {y = data} # already in stats form
  
  nvar = length(y)
  if (missing(model) || is.null(model)) {model = matrix(1,nvar,1)} # sta default model
  
  # extract shrinkage parameters (for information only)
  shrinkage = matrix(0, length(y[[1]]$shrinkage), nvar)
  for (ivar in 1:nvar) {shrinkage[,ivar] = y[[ivar]]$shrinkage}
  
  # decompose model if there are null rows to short circuit result
  s = rowSums(abs(model)); i = which(s != 0); nnonzero = length(i)
  xx = vector("list", nnonzero); fval = 0
  if (qr(model[i,])$rank==nnonzero){
    for (k in 1:nnonzero) {xx[[k]] = matrix(y[[i[k]]]$means,length(y[[i[k]]]$means),1)}
  } else {
    out = jCMRx (y[i], model[i,], E) # call jCMRx for non-zero rows of model
    fval = out$fval; 
    for (k in 1:nnonzero) {xx[[k]] = matrix(out$x[,k],nrow(out$x),1)}
  }
  
  # reassemble output
  x = vector("list", nvar); 
  for (k in 1:nnonzero) {x[[i[k]]] = xx[[k]]}
  
  # replace zero model rows with weighted means
  j = which(s == 0) 
  nzero = length(j)
  if (nzero > 0) {
    for (k in 1:nzero) {
      u = j[k]
      w = y[[u]]$weights
      if (is.matrix(w)) {w = diag(w)}
      w = matrix(w/sum(w),length(w),1)
      z = y[[u]]$means%*%w
      a = rep(z,length(y[[u]]$means))
      x[[u]] = matrix(a,length(y[[u]]$means),1)
      d = y[[u]]$means - t(x[[u]])
      fval = fval + d%*%y[[u]]$weights%*%t(d)
    }
  }
  
  if (fval < tol) {fval = 0} # round down
  
  output = list(x, fval, shrinkage)
  names(output) = c("x", "fval", "shrinkage")
  
  return(output)
}