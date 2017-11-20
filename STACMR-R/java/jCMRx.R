jCMRx <- function(d, model, E) {
  nVar <- length(d)
  nCond <- length(d[[1]]$means)
  
  if (missing("model")) { 
    model <- matrix(1, nVar, 1)
  }

  problemMaker <- new(J("au.edu.adelaide.fxmr.model.CMRxProblemMaker"))
  
  .jcall(problemMaker,returnSig = "V","setModel",.jarray(model, dispatch = T))
  for (i in 1:nVar){
    .jcall(problemMaker,returnSig = "V","addMeanArray",.jarray(d[[i]]$means, dispatch = T))
    .jcall(problemMaker,returnSig = "V","addWeightArray",.jarray(d[[i]]$weights, dispatch = T))
  }

  if (!missing("E") && length(E) > 0) {
    if (is.list(E) && is.list(E[[1]])){
      #3d list, different constrains for each variable
      nE <- length(E)
      for(iVar in 1:nVar){
        Ecur <- E[[iVar]]
        index <- problemMaker$initAdj()
        if (length(Ecur) > 0){
          for (j in 1:length(Ecur)){
            problemMaker$addAdj(nCond, index, as.integer(Ecur[[j]]))
          }
        }
      }
    }else if (is.list(E)){
      #Same list for all variables
      index = problemMaker$initAdj();
      for (j in 1:length(E)){
        problemMaker$addAdj(nCond, index, as.integer(E[[j]]));
      }
      problemMaker$dupeAdj(nVar);
    }
  }
  
  cmrxSolver <- new(J("au.edu.adelaide.fxmr.model.ParCMRxSolver"))
  problem <- problemMaker$getProblem()
  solution <- cmrxSolver$solve(problem)
  fVal <- .jcall(solution,"getFStar",returnSig = "D")
  xMatrix <- t(.jcall(solution,"getXStar",returnSig = "[[D",evalArray=TRUE,simplify = TRUE))
  
  return(list(fval=fVal, x=xMatrix))
}