staMRBN = function (data=list(), partial=list()) {
# fits partial order monotonic regression model to binomial data
# data = nsubj list of nvar lists of ncond x 2 matrices of counts (hits, misses) 
# where ncond = no. of conditions
# partial is the partial order specified as either:
# (a) list, 
# (b) adjacency matrix,
# returns:
# x = nsubj x nvar cell array of ncond model means (proportions)
# fval = nsub-vector of least squares fit of model
# gsquare = nsub vector of g-squared fit of model 
#
# *************************************************************************
# converted from Matlab on 20 September 2016
# *************************************************************************
#
  if (missing(partial) | is.null(partial)) {partial = list()}
  if (is.matrix(partial)) {partial = adj2list(partial)} # convert adjacency matrix to list
  
  y = data
  if (is(data,"data.frame")) {y = gen2listBN(data)} # convert data frame to list of lists
  a = "means" %in% variable.names(as.data.frame(y))
  if (!a) {y = binSTATS(y)}

  #y = binSTATS(data) # get summary stats
  
  x = vector('list', y$ngroup); f = matrix(0, y$ngroup, y$nvar); g = f
  
  for (igroup in 1:y$ngroup) {
    temp = matrix(0, y$ncond, y$nvar)
    for (ivar in 1:y$nvar) {
      means = y[[igroup]][[ivar]]$means
      weights = y[[igroup]][[ivar]]$weights
      count = y[[igroup]][[ivar]]$count
      
      out = jMR (means, weights, partial) # do monotonic regression
      
      out$x[which(out$x < 0)] = 0 # lower limit of 0
      out$x[which(out$x > 1)] = 1 # upper limit of 1
      
      temp[,ivar] = out$x
      f[igroup,ivar] = out$fval
      g[igroup,ivar] = mleBN (out$x, count)
    }
    x[[igroup]] = temp
  }
  f = as.matrix(rowSums(f))
  g = as.matrix(rowSums(g))
  return(list(x=x, fval=f, g.squared=g))
}

mleBN = function (x, count) {
  # calculates g.squared statistic
  nn = sum(count)
  x[which(x <= 0)] = .5/nn
  x[which(x >= 1)] = (nn - .5)/nn
  n = rowSums (count)
  a = t(rbind(x, 1-x))*cbind(n, n)
  g2 = 2*count*log(count/a); g2[which(count <= 0)]=0;
  f = sum(g2); f[which(f < 0)] = 0
  return(f)
}
