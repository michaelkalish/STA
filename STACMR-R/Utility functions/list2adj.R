list2adj <- function(nodes, E=list()) {
  ## Converts a partial order model in list form to an adjacency matrix suitable for monotonic regression
  ##
  ## Args:
  ##   nodes: set of elements (usually 1:n)
  ##   E: List containing the partial order model
  ##
  ## Returns:
  ##   Adjacency matrix for monotonic regression
  # ********************************************************
  # Written 12 September 2016
  # based on matlab code cell2adj.m written by Luke Finlay
  # and original R code written by Wai Keen Vong
  # *********************************************************
  
  if (!is.list(E)) {E <- list(E)}
  
  n <- length(nodes)
  
  ## Initialize adjacency matrix with zeros
  adj <- matrix(0, nrow = n, ncol = n)
  
  ## Fill in adjacency matrix
  
  if (n > 1) {
    if (length(E) > 0) {
      for (i in 1:length(E)) {
        if (length(E[[i]]) >= 2) {
          u = E[[i]]
          for (j in 1:(length(u)-1)) {
            adj[u[j],u[j+1]]=1
          }
        }
      }
    }
  }
  return(adj)
}
