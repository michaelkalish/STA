jMRfits <- function(nsample, y, E=list(), shrink=-1, proc=-1) {
  if (is.data.frame(y)){
    problemMaker <- new(J("au.edu.adelaide.fxmr.model.CMRxFitsGMProblemMaker"))
    problemMaker$setShrink(shrink)
    nVar <- length(unique(y[,2]))
    nCond <- length(unique(y[,3])) * (ncol(y) - 3)
  }else{    
    if (!is.null(y[[1]][["means"]])){
      #Parametric case, probably passed in from staSTATS
      problemMaker <- new(J("au.edu.adelaide.fxmr.model.CMRxFitsProblemMaker"))
      nVar <- length(y)
      nCond <- length(y[[1]]$means)
    }else{
      problemMaker <- new(J("au.edu.adelaide.fxmr.model.CMRxFitsGMProblemMaker"))
      problemMaker$setShrink(shrink)
      nVar <- length(y)
      nCond <- length(y[[1]])
    }
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
            problemMaker$addAdj(as.integer(nCond), index, as.integer(Ecur[[j]]))
          }
        }
      }
    }else if (is.list(E)){
      #Same list for all variables
      index = problemMaker$initAdj();
      for (j in 1:length(E)){
        problemMaker$addAdj(as.integer(nCond), index, as.integer(E[[j]]));
      }
      problemMaker$dupeAdj(nVar);
    }
  }
  
  if (is.data.frame(y)){
    problemMaker$setGM(.jarray(as.matrix(y), dispatch=T));
    sol <- problemMaker$solve(as.integer(nsample),as.integer(proc),FALSE,TRUE);
  }else{
    if (!is.null(y[[1]][["means"]])){
      #Parametric case, probably passed in from staSTATS
      ns = matrix(1,nVar)
      for (iVar in 1:nVar){
        problemMaker$addMeanArray(y[[iVar]]$means)
        w <- y[[iVar]]$weights
        if (is.null(w)){
          w <- diag(length(y[[iVar]]$means));
        }
        problemMaker$addWeightArray(.jarray(w, dispatch = T))
        problemMaker$addCov(.jarray(y[[iVar]]$regcov, dispatch = T))
        if (length(y[[iVar]]$n) == 1){
          ns[[iVar]] <- y[[iVar]]$n
        }else{
          ns[[iVar]] <- y[[iVar]]$n[[1]][[1]]
        }
      }
      problemMaker$setN(as.integer(ns))
      problem <- problemMaker$getProblem()
      
      sol <- new(J("au.edu.adelaide.fxmr.model.CMRxFits"),as.integer(nsample),problem,shrink,as.integer(proc))
    }else{
      for (group in 1:nCond){
        for (iVar in 1:nVar){
          problemMaker$addCell(as.integer(iVar),as.integer(group),.jarray(as.matrix(y[[iVar]][[group]]), dispatch = T))
        }
      }
      sol <- problemMaker$solve(as.integer(nsample),as.integer(proc),FALSE,TRUE);
    }
  }
  
  return(list(p = sol$getP(), fits = sol$getFits(),datafit = sol$getDataFit()))
}