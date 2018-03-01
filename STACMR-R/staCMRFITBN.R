staCMRFITBN <- function (data=NULL, partial = list(), nsample=1, model=NULL, proc=-1, approx=0) {
  # input:
  # nsample = no. of Monte Carlo samples (about 10000 is good)
  # data = data structure (cell array or general)
  # model is a nvar * k matrix specifying the linear model, default = ones(nvar,1))
  # partial = optional partial order model e.g. E={[1 2] [3 4 5]} indicates that
  # condition 1 <= condition 2 and condition 3 <= condition 4 <= condition 5
  # default = none (empty)
  
  # output:
  # p = empirical p-value
  # datafit = observed fit of monotonic (1D) model
  # fits = nsample vector of fits of Monte Carlo samples (it is against this
  # distribution that datafit is compared to calculate p)
  # *************************************************************************
  # converted from matlab 7 February 2018
  # *************************************************************************
  
  if (is(data,"data.frame")) {
    y = gen2listBN (data) # convert from general format to list format
  } else {y = data} 
  
  a = "means" %in% variable.names(as.data.frame(y))
  if (!a) {y = binSTATS(y)}
  
  tol <- 10e-6
  if (missing(partial)) {partial = list()}
  
  nvar =length(y)
  model = NULL
  # if (missing(model) | is.null(model)) {model = matrix(1,nvar,1)} # sta default model
  
  if (!is.list(partial)) {partial = adj2list(partial)} # convert from adjacency matrix to list
  
  output = jCMRxBNfits(nsample, y, partial, model, proc, approx) # call java program
  
  output$fits[which(output$fits <= tol)] = 0;
  
  return (output)
}