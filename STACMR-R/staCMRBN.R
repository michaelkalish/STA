staCMRBN = function (data=list(), partial=list()) {
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
  y = data
  
  if (missing(partial) | is.null(partial)) {partial = list()}
  if (is.matrix(partial)) {partial = adj2list(partial)} # convert adjacency matrix to list
  
  if (is(y,"data.frame")) {y = BNframe2list(y)} # convert data frame to list of lists
  
  out = jCMRBN (y, partial) # fit model

  #g = matrix(0, y$ngroup, y$nvar) # calculate g-squared
  #for (igroup in 1:y$ngroup) {
  #  x = out$x[[igroup]]
  #  for (ivar in 1:y$nvar) {
  #    g[igroup,ivar] = mleBN (x[,ivar], y[[igroup]][[ivar]])
  #  }
  #}
  #g2 = as.matrix(rowSums(g), y$ngroup, y$nvar)
  return(list(x=out$x, fval=out$f, g.squared=out$g.squared))
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