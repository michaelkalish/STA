binSTATS = function(data=list()) {
# returns statistics on binomial data
# composed of NGROUP list of NVAR list of matrices (unless NGROUP=1) where
# each element is a NCOND x 2 matrix (hits, misses)
# output is a NGROUP list of NVAR list consisting of:
# count = hits and misses
# means = observed means
# n = number of observations
# weights = weights for monotonic regression
#
# *************************************************************************
# converted from Matlab on 20 September 2016
# *************************************************************************
#
  y = data;
  
  if (!is.null(nrow(y))) {y = BNframe2list(y); ngroup=y$ngroup; nvar=y$nvar} # convert data frame to list of lists
  
  if (!with(y, exists('ngroup'))) {ngroup = length(y); nvar = length(y[[1]])
  } else {ngroup = y$ngroup; nvar = y$nvar} # calculate ngroup and nvar
  
  output = vector("list", ngroup) # define output list
  
  y1 = y[[1]][[1]] # y1 = first entry in list
  ncond = nrow(y1)
  
  if (with(y1, exists('means'))) {output = y # already in stats form
  } else {

    for (igroup in 1:ngroup) {
      temp = vector("list", nvar)
      for (ivar in 1:nvar) {   
        b.count = y[[igroup]][[ivar]];
        b.n = rowSums(b.count)
        b.means = b.count[,1]/b.n
        b.weights = diag(b.n)
        
        out=list(b.count, b.means, b.n, b.weights)
        names(out)=c('count', 'means','n','weights')
        temp[[ivar]] = out
      }
      output[[igroup]] = temp
    }
  }
  output$ngroup = ngroup; output$nvar = nvar; output$ncond = ncond
  return(output)
}

  