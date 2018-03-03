staCMR <- function (data, partial = list(), shrink=-1, approx=F) {
  # function [x, f, shrinkage] = staCMR (data, model, partial, shrink)
  # wrapper function for staCMRx
  # Multidimensional CMR
  # data is cell array of data or structured output from staSTATS 
  # partial is partial order
  # shrink is parameter to control shrinkage of covariance matrix;
  # 0 = no shrinkage; 1 = diagonal matrix; -1 = calculate optimum
  # returns:
  # x = best fitting CMR values to y-means
  # fval = fit statistic
  # shrinkage = estimated shrinkage of covariance matrix
  
  # *************************************************************************
  # 21 August 2017
  # *************************************************************************
  #
  # set up defaults
  tol <- 10e-5
  if (missing(partial)) {partial = list()}
  if (missing(shrink)) {shrink = -1}
  
  output = staCMRx (data, model=NULL, E=partial, shrink=shrink, tolerance=0, proc=-1, approx=approx)
  
  return (output)
}