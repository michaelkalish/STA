staCMRFIT <- function (data=NULL, partial = list(), nsample=1, shrink=-1, approx=F) {
# input:
  # nsample = no. of Monte Carlo samples (about 10000 is good)
  # data = data structure (cell array or general)
  # model is a nvar * k matrix specifying the linear model, default = ones(nvar,1))
  # partial = optional partial order model e.g. E={[1 2] [3 4 5]} indicates that
  # condition 1 <= condition 2 and condition 3 <= condition 4 <= condition 5
  # default = none (empty)
  # shrink is parameter to control shrinkage of covariance matrix (if input is not stats form);
  # 0 = no shrinkage; 1 = diagonal matrix; -1 = calculate optimum, default = -1
# output:
  # p = empirical p-value
  # datafit = observed fit of monotonic (1D) model
  # fits = nsample vector of fits of Monte Carlo samples (it is against this
  # distribution that datafit is compared to calculate p)
  # *************************************************************************
  # converted from matlab 18 September 2016
  # *************************************************************************

  if (is(data,"data.frame")) {
    y = gen2list (data) # convert from general format to list format
  } else {y = data} 
  
  tol <- 10e-6
  if (missing(partial)) {partial = list()}
  if (missing(shrink)) {shrink = -1}
  proc = -1
  cheapP = F
  mrTol = 0
  seed = -1
  
  nvar =length(y[[1]])
  model = NULL
  if (missing(model) | is.null(model)) {model = matrix(1,nvar,1)} # sta default model
  
  if (!is.list(partial)) {partial = adj2list(partial)} # convert from adjacency matrix to list
  
  output = jCMRfitsx(nsample, y, model, partial, shrink, proc, cheapP, approx, mrTol, seed) # call java program
  
  output$fits[which(output$fits <= tol)] = 0;
  
  return (output)
}