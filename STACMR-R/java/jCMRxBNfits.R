jCMRxBNfits <- function(nsample, data, E=list(),model,proc=-1, approximate=FALSE) {
  d <- data
  
  if (with(d, exists('ngroup'))) {
    nSubj = d$ngroup
    nVar = d$nvar
  } else {
    nSubj = length(d)
    nVar = length(d[[1]])
  }
  
  if (missing("model")) { 
    model <- matrix(1, nVar, 1)
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
  print
  fObj <- new(J("au.edu.adelaide.fxmr.model.bin.BinCMRxFits"),as.integer(nsample), problem, as.integer(proc),approximate,FALSE)

  p <- fObj$getP()
  datafit <- fObj$getBaseFitDiff()
  fits <- t(.jevalArray(fObj$getFits(),simplify=T))
  
  return(list(p=p,datafit=datafit,fits=fits))
}