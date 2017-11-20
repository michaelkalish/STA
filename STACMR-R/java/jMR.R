jMR <- function(y, W, E) {
  nvar<-length(y)
  yM <- data.matrix(y)
  problemMaker <- new(J("au.edu.adelaide.fxmr.model.mr.MRProblemMaker"), .jarray(y))
  
  if (missing(W)){
    W <- diag(nvar)
  }
  .jcall(problemMaker,returnSig = "V","setWeightArray",.jarray(W, dispatch = T))
  
  if (!missing("E")) {
    for (e in E){
      .jcall(problemMaker,returnSig = "V","addConstraints",.jarray(as.integer(e)))
    }
  }
  
  mrSolver <- new(J("au.edu.adelaide.fxmr.model.mr.MRSolverAJOptimiser"))

  problem <- problemMaker$getProblem()
  
  solution <- mrSolver$solve(problem)
 
  fVal <- .jcall(solution,"getfVal",returnSig = "D")
  xVector <- t(.jcall(solution,"getxVector",returnSig = "[D",evalArray=TRUE,simplify = TRUE))
  
  return(list(fval=fVal, x=xVector))
}