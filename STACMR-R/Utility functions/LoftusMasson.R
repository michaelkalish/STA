LoftusMasson <- function(data) {
  y = data
  nr = nrow(y); nc = ncol(y)
  m.cols = colMeans(y); y.cond = t(matrix(rep(m.cols,nr),nc,nr))
  m.rows = rowMeans(y); y.subj = matrix(rep(m.rows,nc),nr,nc)
  m.mean = mean(colMeans(y)); y.mean = matrix(rep(m.mean,nr*nc),nr,nc)
  ya = y - y.cond - y.subj + y.mean
  ss = sum(ya^2); df = (nr-1)*(nc-1)
  ms.resid = ss/df
  r = diag(rep(ms.resid,nc))
  return(r)
}