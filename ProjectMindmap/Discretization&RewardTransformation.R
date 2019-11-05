##----------------------------------------------------
## Discretization of continuous phase type distribution
##----------------------------------------------------
## Name: discretization
## Purpose: Computing the discrete phase-type parameters
##          for certain summary statistics for standard
##          coalescent model, that arise by sprinkling
##          mutations on a tree given by a certain rate
##          (Similar to Theorem 3.5 in [HSB])
## Input:
## object = a continuous phase-type distribution object
## mutrate = the mutation rate
## Output: 
## A discrete phase-type distribution object
##----------------------------------------------------
discretization <- function(object, mutrate){
  
  discphasetype(initDist = object$initDist,
                P.mat = solve(diag(nrow = nrow(object$T.mat))-2/mutrate*object$T.mat))
  
}

##----------------------------------------------------
## The reward-transformed distribution
##----------------------------------------------------
## Name: RewTransDistribution
## Purpose: An implementation of Theorem 3.1.33 in 
##         [BN](2017). That is, the function computes
##         the reward-transformed distribution
## Input:
## object = a continuous phasetype distribution object
## rewards = the non-negative reward vector
## Output: 
## A new contphasetype object with the reward transformed
## parameters.
##----------------------------------------------------
RewTransDistribution <- function(object, rewards){
  
  if(sum(rewards < 0) >0 | sum(rewards)==0 ){
    
    stop("The reward vector has to be non-negative, and contain at least one positve number!")
    
  }else{
    
    #We extract the initial distribution and subintensity matrix from the object
    initDist <- object$initDist
    T.mat <- object$T.mat
    
    ## We define the sets S^+ and S^0
    S.plus <- which(rewards>0)
    S.zero <- which(rewards==0)
    ## as well as the matrix Q
    Q.mat <- -T.mat/diag(T.mat)
    diag(Q.mat) <- 0
    
    ## We rearrange the states according to the 
    ## sets S^+ and S^0:
    Q.matpp <- Q.mat[S.plus,S.plus]
    Q.matp0 <- Q.mat[S.plus,S.zero]
    Q.mat0p <- Q.mat[S.zero,S.plus]
    Q.mat00 <- Q.mat[S.zero,S.zero]
    
    ## Now we can define the transition matrix
    ## for the new markov chain
    if(length(S.zero)==1){
      
      P.mat <- Q.matpp+(1-Q.mat00)^(-1)*Q.matp0%*%t(Q.mat0p)
      
    }else if(length(S.zero)==0){
      P.mat <- Q.matpp
      
    }else{
      P.mat <- Q.matpp+Q.matp0%*%solve(diag(1,nrow = nrow(Q.mat00))-
                                         Q.mat00)%*%Q.mat0p
    }
    
    ## We define the new exit vector
    p.vec <- 1-rowSums(P.mat)
    
    ## We also split the original initial distribution
    ## into pi = (pi^+,pi^0)
    pi.vecp <- initDist[S.plus]
    pi.vec0 <- initDist[S.zero]
    
    ## Then, the initial distribution of the new
    ## Markov chain is given by
    if(length(S.zero)==1){
      
      newInitDist <- pi.vecp + (1-Q.mat00)*pi.vec0%*%t(Q.mat0p)
      
    }else if(length(S.zero)==0){
      
      newInitDist <- pi.vecp 
      
    }else{
      newInitDist <- pi.vecp + pi.vec0%*%solve(diag(1,nrow = nrow(Q.mat00))
                                             -Q.mat00)%*%Q.mat0p
    }
    
    ## Now we can define the subintensity matrix
    ## newT.mat as
    if(length(S.plus)==1){
      
      newT.mat <- -T.mat[S.plus,S.plus]/rewards[S.plus]*P.mat
      exitrate <- -T.mat[S.plus,S.plus]/rewards[S.plus]*p.vec
    }else{
      
      newT.mat <- -diag(T.mat[S.plus,S.plus])/rewards[S.plus]*P.mat
      exitrate <- -diag(T.mat[S.plus,S.plus])/rewards[S.plus]*p.vec
    }
    
    
    diag(newT.mat) <- -rowSums(newT.mat)-exitrate
    
  }
  return(contphasetype(initDist = newInitDist, T.mat = newT.mat))
}

##----------------------------------------------------
## Moments of continuous phase-type distribution.
##----------------------------------------------------
## Name: LaplacePhasetype
## Purpose: An implementation of Corollary 3.1.18 in [BN].
##          Computes a certain moment of a continuous
##          phase-type distribution, based on the 
##          Laplace-transform.
## Input:
## object = a continuous phase-type distribution object
## i = the order of the desired moment
## Output: 
## The i'th moment of the continuous phase-type object
##----------------------------------------------------
LaplacePhaseType <- function(object, i){
    
    initDist = object$initDist
    T.mat = object$T.mat
    
    return((-1)^i*factorial(i)*sum(initDist%*%(solve(Tmat)%^%i)))
}