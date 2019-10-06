
## Density
dphasetype <- function(u, initDist, Tmat){
  
  evec <- replicate(nrow(Tmat),1)
  return(c(-initDist%*%expm(u*Tmat)%*%Tmat%*%evec))
}

## distribution
pphasetype <- function(u, initDist, Tmat){
  
  evec <- replicate(nrow(Tmat),1)
  return(c(1-initDist%*%expm(u*Tmat)%*%evec))
}

## random numbers from the distribution
rphasetype <- function(n, initDist, Tmat){
  
  tau <- replicate(n,0)
  states <- 1:(length(initDist)+1)
  
  evec <- replicate(nrow(Tmat),1)
  
  TransMat <- cbind(Tmat, -Tmat%*%evec)
  TransMat <- rbind(TransMat, replicate(ncol(TransMat),0))
  for(i in 1:n){
    
    initState <- sample(states[-(length(initDist)+1)], size = 1, prob = initDist )
    
    x <- initState
    while(x != tail(states, 1)){
      
      tau[i] <- tau[i] + rexp(n=1, rate = -diag(TransMat)[x])
      x <- sample(states[-x], size = 1, prob = -TransMat[x,-x]/TransMat[x,x] )
      
    }
    
  }
  
  return(tau)
}


## Laplace transform
LaplacePhaseType <- function(initDist, Tmat, i){
  
  ## For n=2, we have that the initial distribution is piMRCA = 1 and
  ## the Transition probability matrix is Tmat = -1 (for both T_MRCA and T_Total),
  ## hence the Laplace transform is given by
  if(length(initDist) == 1){
    
    res <- factorial(i) ## = (-1)^n*n!*1*(-1)^(-n)*1 = (-1)^n*n!*pi*T^(-n)*e
  }else{
    
    evec <- replicate(nrow(Tmat),1)
    res <- (-1)^i*factorial(i)*initDist%*%(solve(Tmat)%^%i)%*%evec
  }
  
  return(res)
}


##----------------------------------------------------
## The reward-transformed distribution
##----------------------------------------------------
## Name: RewTransDistribution
## Purpose: An implementation of Theorem 3.1.33 in 
##         [BN](2017). That is, the function computes
##         the reward-transformed distribution
## Input:
## pi.vec = the initial distribution of the original
##          phase type distribution
## Tmat = the subintensity matrix from the original
##        phase type distribution
## r.vec = the non-negative reward vector
## Output: 
## alpha.d = the defect size
## alpha.vec = the initial distribution of the reward-
##             transformed distribution
## T.rewmat = the subintensity matrix of the reward-
##            transformed distribution
##----------------------------------------------------

RewTransDistribution <- function(pi.vec,Tmat, r.vec){
  
  if(sum(r.vec < 0)>0){
    
    stop("The reward vector r.vec has to be non-negative!")
    
  }else{
    
    ## We define the sets S^+ and S^0
    S.plus <- which(r.vec>0)
    S.zero <- which(r.vec==0)
    ## as well as the matrix Q
    Q.mat <- -Tmat/diag(Tmat)
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
    pi.vecp <- pi.vec[S.plus]
    pi.vec0 <- pi.vec[S.zero]
    
    ## Then, the initial distribution of the new
    ## Markov chain is given by
    if(length(S.zero)==1){
      
      alpha.vec <- pi.vecp + (1-Q.mat00)*pi.vec0%*%t(Q.mat0p)
      
    }else if(length(S.zero)==0){
      
      alpha.vec <- pi.vecp 
      
    }else{
      alpha.vec <- pi.vecp + pi.vec0%*%solve(diag(1,nrow = nrow(Q.mat00))
                                             -Q.mat00)%*%Q.mat0p
    }
    ## The defect size alpha_d+1 is given by
    alpha.d <- 1- sum(alpha.vec)
    d <- length(S.plus)
    
    ## Now we can define the subintensity matrix
    ## T.rewmat as
    if(length(S.plus)==1){
      
      T.rewmat <- -Tmat[S.plus,S.plus]/r.vec[S.plus]*P.mat
      t.vec <- -Tmat[S.plus,S.plus]/r.vec[S.plus]*p.vec
    }else{
      
      T.rewmat <- -diag(Tmat[S.plus,S.plus])/r.vec[S.plus]*P.mat
      t.vec <- -diag(Tmat[S.plus,S.plus])/r.vec[S.plus]*p.vec
    }
    
    
    diag(T.rewmat) <- -rowSums(T.rewmat)-t.vec
    
  }
  return(list(alpha.vec = alpha.vec, alpha.d = alpha.d, T.rewmat = T.rewmat))
}

##----------------------------------------------------
## The distribution of the number S of segregating 
## sites
##----------------------------------------------------
## Name: DistSegregatingSites
## Purpose: Compute the probability P(S=k) for a 
##          given number k
## Input:
## pi.vec = the initial distribution of the
##          phase type distribution
## Tmat = the subintensity matrix from the
##        phase type distribution
## theta = the mutation rate (theta/2)
## k = a non-negative number
## Output: 
## The probability P(S=k)
##----------------------------------------------------

DistSegregatingSites <- function(pi.vec, Tmat, theta, k){
  
  P <- solve(diag(x=1, nrow = nrow(Tmat))-2/theta*Tmat)
  p <- replicate(n=nrow(Tmat),1)-P%*%replicate(n=nrow(Tmat),1)
  
  return(pi.vec%*%(P%^%k)%*%p)
}