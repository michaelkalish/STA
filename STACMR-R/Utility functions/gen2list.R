gen2list = function (data=NULL, varnames=list()) {
# gen2cell(data)
  # R version of gen2cell.m
# converts data in "general format" to list format suitable for input to staSTATS
# general format is defined as:
  # column 1 = subject number (nsub)
  # column 2 = between-subjects condition (ngroup)
  # column 3 = dependent variable (nvar)
  # columns 4 to end = values for each within-subjects condition (ncond)
  # output is ngroup x nvar list in which each element is an nsub x ncond matrix of values
  #
  # *************************************************************************
  # written 12 September 2016
  # revised 9 March 2017 to remove missing within variables in a group
  # revised 22 August 2017 to add variable names
  # *************************************************************************
  #
  
  group = data[,2]; ugroup = sort(unique(group)); ngroup = length(ugroup)
  var = data[,3]; uvar = sort(unique(var)); nvar = length(uvar)
  within = as.matrix(data[,4:ncol(data)])
  
  y = vector("list",ngroup)
  for (igroup in 1:ngroup) {
    temp = vector("list", nvar)
    for (ivar in 1:nvar){
      k = which(group==ugroup[igroup] & var==uvar[ivar])
      a = as.matrix(within[k,])
      # delete any variables that all all missing
      n = colSums(is.na(a)); k=which(n==nrow(a)); if (length(k) > 0) {a = a[,-k]}
      # store in 2D list
      y[[igroup]][[ivar]]=a
      if (!missing(varnames)) {names(y[[igroup]][[ivar]])=varnames}
    }
  }
  return (y)
}
  