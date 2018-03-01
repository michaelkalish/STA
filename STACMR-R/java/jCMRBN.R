jCMRBN <- function(data=list(), E=list(), model=NULL, approximate=FALSE) {
  d <- data
  
  if (with(d, exists('ngroup'))) {
    nSubj = d$ngroup
    nVar = d$nvar
  } else {
    nSubj = length(d)
    nVar = length(d[[1]])
  }
  
  if (missing(model) | is.null(model)) {model <- matrix(1, nVar, 1)}
  
  problemMaker <- new(J("au.edu.adelaide.fxmr.model.bin.BinCMRProblemMaker"), nSubj, nVar)
  problemMaker$setModel(.jarray(model));
  
  for (v in 1:nVar){
    for (s in 1:nSubj){
      dcur = d[[s]][[v]]
      if (with(dcur, exists('count'))) {
        df = as.data.frame(dcur$count, nrow=2)
      }else{
        df = as.data.frame(dcur, nrow=2)
      }
      problemMaker$setElement(as.integer(s-1), as.integer(v-1), as.integer(unlist(df)))
    }
  }

  if (!missing("E") && is.list(E) && length(E) > 0) {
    #3d list, different constrains for each variable
    for(e in E){
      problemMaker$addRangeSet(.jarray(as.integer(e)))
    }
  }
  
  solver <- new(J("au.edu.adelaide.fxmr.model.bin.BinCMRxSolver"))
  solver$setOnlyFeas(as.logical(approximate));
  ps <- problemMaker$getProblems()
  solutions <- solver$solve(ps)
  
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
