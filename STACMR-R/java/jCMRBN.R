jCMRBN <- function(data=list(), E=list()) {
  d <- data
  
  if (with(d, exists('ngroup'))) {
    nSubj = d$ngroup
    nVar = d$nvar
  } else {
    nSubj = length(d)
    nVar = length(d[[1]])
  }
  
  problemMaker <- new(J("au.edu.adelaide.fxmr.model.bin.BinCMRProblemMaker"), nSubj, nVar)

  if (with(d, exists('ngroup'))) {
    for (v in 1:nVar){
      for (s in 1:nSubj){
        #Minus 1 to zero index
        df = as.data.frame(d[[s]][[v]], nrow=2)
        problemMaker$setElement(as.integer(s-1), as.integer(v-1), as.integer(unlist(df)))
      }
    }
  } else {
    for (v in 1:nVar){
      for (s in 1:nSubj){
        df = as.data.frame(d[[s]][[v]], nrow=2)
        problemMaker$setElement(as.integer(s-1), as.integer(v-1), as.integer(unlist(df)))
      }
    }
  }
  
  if (!missing("E") && is.list(E) && length(E) > 0) {
    #3d list, different constrains for each variable
    for(e in E){
      problemMaker$addRangeSet(as.integer(e))
    }
  }
  
  solver <- new(J("au.edu.adelaide.fxmr.model.bin.BinCMRSolver"))
  ps <- problemMaker$getProblems()
  solutions <- solver$solve(problemMaker$getProblems())
  
  x <- vector("list", nSubj)
  f <- matrix(0,nSubj,1)
  g2 <- matrix(0,nSubj,1)
  iter <- vector("list", nSubj)
  i <- 1
  for(s in .jevalArray(solutions)){
    x[[i]] = t(.jevalArray(s$getXStar(),simplify = T))
    f[i] = s$getFStar()
    g2[i] = s$getG2Star()
    iter[[i]] = s$getIter()
    i <- i + 1
  }
  return(list(x=x,fval=f,g.squared=g2))
}
