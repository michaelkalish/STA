adj2list <- function(adj) {
  ## Converts a partial order model from adjacency matrix to list form 
  ##
  ## Args:
  ##   adj: Adjacency matrix for monotonic regression
  ##
  ## Returns:
  ##   List containing the partial order model
  # ********************************************************
  # Written 12 September 2016
  # *********************************************************

  u = which(adj==1,arr.ind=T); n = nrow(u)
  E = NULL
  if (n > 0) {
    E = vector("list", n)
    for (i in 1:n) {
      E[[i]] = c(u[i,1],u[i,2])
    } 
  }
  return (E)
}
