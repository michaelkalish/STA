jMRBNfits <- function(nsample, data, E=list(),proc=-1) {
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
        df = as.data.frame(d[[s]][[v]]$count, nrow=2)
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
  
  problem <- problemMaker$getBaseProblem()
  fObj <- new(J("au.edu.adelaide.fxmr.model.bin.BinCMRFits"),as.integer(nsample), problem, as.integer(proc),TRUE)
  
  p <- fObj$getP()
  datafit <- fObj$getBaseFitDiff()
  fits <- t(.jevalArray(fObj$getFits(),simplify=T))
  
  return(list(p=p,datafit=datafit,fits=fits))
}
