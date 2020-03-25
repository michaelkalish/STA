staMRFIT <- function (data=NULL, partial = list(), nsample=1, shrink=-1) {
# input:
  # nsample = no. of Monte Carlo samples (about 10000 is good)
  # data = data structure (cell array or general)
  # partial = optional partial order model e.g. E={[1 2] [3 4 5]} indicates that
  # condition 1 <= condition 2 and condition 3 <= condition 4 <= condition 5
  # default = none (empty)
  # shrink is parameter to control shrinkage of covariance matrix (if input is not stats form);
  # 0 = no shrinkage; 1 = diagonal matrix; -1 = calculate optimum, default = -1
# output:
  # p = empirical p-value
  # datafit = observed fit of partial order model
  # fits = nsample vector of fits of Monte Carlo samples (it is against this
  # distribution that datafit is compared to calculate p)
  # pars = nvar list of bootstrap means nsample x ncond
  # *************************************************************************
  # converted from matlab 7 February 2018
  # bootstrap parameter estimates added 24 March 2020
  # *************************************************************************

  if (is(data,"data.frame")) {
    y = gen2list (data) # convert from general format to list format
  } else {y = data} 
  
  tol <- 10e-6
  if (missing(partial)) {partial = list()}
  if (missing(shrink)) {shrink = -1}
  
  nvar = length(y[[1]])
  ngroup = length(y); ncond = length(y[[1]][[1]][1,]) * ngroup
  
  if (!is.list(partial)) {partial = adj2list(partial)} # convert from adjacency matrix to list
  
  output = jMRfits(nsample, y, partial, shrink); # call java program

  output$fits[which(output$fits <= tol)] = 0;
  output$datafit[which(output$datafit <= tol)] = 0;
  
  # unpack bootstrap means
  z = array(0,dim=c(nvar,ncond,nsample))
  for (isample in 1:nsample) {
    z[,,isample] = .jevalArray(output$pars[[isample]],simplify=T)
  }
  a = vector("list", nvar)
  z = aperm(z,c(3,2,1))
  for (ivar in 1:nvar){
    a[[ivar]] = z[,,ivar]
  }
  output$pars = a;
  
  return (output)
}