##----------------------------------------------------
## EM-estimation for discrete phase type distribution
##----------------------------------------------------
## Name: emdiscphasetypeestimation
## Purpose: An implementation of section 13.2 in [BN]
##          regarding EM-estimation of the parameters
##          of a discrete phase-type distribution.
## Input:
## observations = a vector of i.i.d. realization of the same
##                DPH_p-distribution
## initDist0 = the initial distribution that the first
##             run of the algorithm uses
## P.mat0 = the subtransition matrix that the first
##          run of the algorithm uses
## n = the number of times the algorithm will run
## Output: 
## A discphasetype object with the estimated parameters.
##----------------------------------------------------

emdiscphasetypeestimation <- function(observations, initDist0, P.mat0, n){
  
  p <- length(initDist0)
  nobs <- length(observations)

  Kmatlist <- function(y, initDist, P.mat){
    
    pvec <- 1 - rowSums(P.mat)
    
    Kmatlist <- list(P.mat%^%(y-2)%*%pvec%*%t(initDist))
    
    for(k in 1:(y-2)){
      Kmatlist[[k+1]] <- Kmatlist[[k]] + P.mat%^%(y-2-k)%*%pvec%*%t(initDist)%*%(P.mat%^%(k))
    }
    
    return(Kmatlist)
  }
  
  estepandmstep <- function(initDist, P.mat){
    
    pvec <- 1 - rowSums(P.mat)
    
    listofKmats <- Kmatlist(max(observations), initDist, P.mat)
    
    normalizer <- 1
    bvec <- rep(0,p)
    Nmat <- matrix(0, nrow = p, ncol = p)
    Nvec <- rep(0,p)
    
    for(i in 1:nobs){
      normalizer <- as.numeric(initDist %*% (P.mat %^% (observations[i]-1)) %*% pvec)
      
      bvec <- bvec + initDist * (P.mat %^%(observations[i]-1) %*% pvec)/normalizer
      
      if(observations[i] > 1){
        
        Nmat <- Nmat + P.mat * t(listofKmats[[observations[i]-1]])/normalizer
        
      }
        
      Nvec <- Nvec + (initDist %*% (P.mat %^%(observations[i]-1))) * pvec/normalizer
      
    }
    
    newinitDist <- as.vector(bvec/sum(bvec))
    
    sumconds <- sum(Nmat) + sum(Nvec)
    
    newP.mat <- Nmat/sumconds
    
    diag(newP.mat) <- 1 - (rowSums(newP.mat) + Nvec/sumconds)
    
    return(discphasetype(newinitDist,newP.mat))
  }
  
  initDist <- initDist0
  P.mat <- P.mat0
    
  for(j in 1:n){
    newparms <- estepandmstep(initDist, P.mat)
    initDist <- newparms[[1]]
    P.mat <- newparms[[2]]
  }
  return(newparms)
}
 


discphasetypetrajectorybasedestimate <- function(listoftrajectories){
  
  # A helping function that calculates the number of
  # transition between the states and store the initial
  # and exiting states.
  trajectorysummarystat <- function(trajectory){
    
    p <- max(trajectory) - 1
    
    trajectorymat <- cbind(trajectory[-length(trajectory)],trajectory[-1])
    Nmat <- matrix(0, ncol = p, nrow = p)
    for(i in 1:p){
      for(j in 1:p){
        Nmat[i,j] <- length(which(trajectorymat[,1] == i & trajectorymat[,2] == j))
      }
    }
    
    return(cbind(1:p == trajectory[1],Nmat,1:p == trajectory[length(trajectory)-1]))
  }
  
  N <- length(listoftrajectories)
  tmp <- lapply(listoftrajectories, trajectorysummarystat)
  p <- nrow(tmp[[1]])
  totalmat <- matrix(0, ncol = p + 2, nrow = p)
  for(i in 1:N){
    totalmat <- totalmat + tmp[[i]]
    }
  initDistEst <- totalmat[,1]/N
  transmatest <- totalmat[,-1]/sum(totalmat[,-1])
  P.matEst <- transmatest[,-(p+1)]
  diag(P.matEst) <- 1 - rowSums(transmatest)
  return(discphasetype(initDist = initDistEst, P.mat = P.matEst))
  
}