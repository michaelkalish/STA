staMRFITBN <- function (data=NULL, partial = list(), nsample=1, proc=-1, approx=0) {
# input:
  # nsample = no. of Monte Carlo samples (about 10000 is good)
  # data = data structure (cell array or general)
  # partial = optional partial order model e.g. E={[1 2] [3 4 5]} indicates that
  # condition 1 <= condition 2 and condition 3 <= condition 4 <= condition 5

# output:
  # p = empirical p-value
  # datafit = observed fit of monotonic (1D) model
  # fits = nsample vector of fits of Monte Carlo samples (it is against this
  # distribution that datafit is compared to calculate p)
  # pars = ngroup x nvar list of bootstrap means nsample x ncond
  # *************************************************************************
  # converted from matlab 21 February 2018
  # bootstrap parameter estimates added 24 March 2020
  # *************************************************************************

  if (is(data,"data.frame")) {
    y = gen2listBN (data) # convert from general format to list format
  } else {y = data} 
  
  a = "means" %in% variable.names(as.data.frame(y))
  if (!a) {y = binSTATS(y)}
  
  tol <- 10e-6
  if (missing(partial)) {partial = list()}
  
  nvar = y$nvar
  ncond = y$ncond
  ngroup = y$ngroup
  
  if (!is.list(partial)) {partial = adj2list(partial)} # convert from adjacency matrix to list
  
  output = jMRBNfits(nsample, y, partial, proc) # call java program
  
  output$fits[which(output$fits <= tol)] = 0;
  output$datafit[which(output$datafit <= tol)] = 0;
  
  # unpack bootstrap means
  z = array(0,c(nvar,ncond,nsample,ngroup))
  for (isample in 1:nsample){
    p = .jevalArray(output$pars[[isample]])
    for (igroup in 1:ngroup){
      z[,,isample,igroup] = .jevalArray(p[[igroup]], simplify=T)
    }
  }
  a = vector("list", ngroup)
  z = aperm(z, c(3,2,4,1))
  for (igroup in 1:ngroup){
    temp = vector("list", nvar)
    for (ivar in 1:nvar){
      temp[[ivar]] = z[,,igroup,ivar]
    }
    a[[igroup]] = temp;
  }
  output$pars = a
  
  return (output)
}