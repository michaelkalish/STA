staMR <- function(data=list(), partial = list(), shrink=-1) {
  # function [xPrime, fit, shrinkage] = staMR (data, partial, shrink)
  # fits monotonic regression model to data according to partial order
  # data is list of lists of data or structured output from staSTATS
  # partial is partial order in list format
  # shrink is parameter to control shrinkage of covariance matrix;
  # 0 = no shrinkage; 1 = diagonal matrix; -1 = calculate optimum
  # returns:
  #   x = best fitting MR values to y-means
  #   f = total fit statistic
  #   shrinkage = shrinkage from applying staSTATS
  # *************************************************************************
  # modified from matlab 13 September 2016
  # *************************************************************************
  #

  # set up defaults
  tol <- 10e-5
  if (missing(partial)) {partial = list()}
  if (missing(shrink)) {shrink = -1}
  
  # get stats from data (depending on its form)
  if (is(data,"data.frame")) {
    y = gen2list (data) # convert from general format
    y = staSTATS (y, shrink) # get stats
  } else if (is.null(data[[1]]$means)) {y = staSTATS(data, shrink) # in list form, get stats
  } else {y = data} # already in stats form
  
  # convert partial order to list if in adjacency matrix form
  if (is(partial,"matrix")) {partial = adj2list(partial)}
  
  # extract shrinkage parameters (for information only)
  nvar = length(y)
  shrinkage = matrix(0, length(y[[1]]$shrinkage), nvar)
  for (ivar in 1:nvar) {shrinkage[,ivar] = y[[ivar]]$shrinkage}

  # do MR for each dependent variable
  xPrime = vector("list", nvar)
  fit = matrix(0, nvar, 1)
  for (ivar in 1:nvar) {
    out = jMR (y[[ivar]]$means, y[[ivar]]$weights, partial)
    xPrime[[ivar]] = out$x
    fit[ivar] = out$fval
  }
  fval = sum(fit)
  if (fval < tol) {fval = 0} # round down
  
  for (i in 1:nvar) {xPrime[[i]]=matrix(xPrime[[i]],length(xPrime[[i]]),1)}
  output = list(xPrime, fval, shrinkage)
  names(output) = c("x", "fval", "shrinkage")
  
  return(output)
}