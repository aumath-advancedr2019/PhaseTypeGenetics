##----------------------------------------------------
## The quantile function for of a continuous phase-type distribution
##----------------------------------------------------
## Name: qcontphasetype
## Purpose: Computing the quantile function of a continuous
##          phase-type distribution by numeric inversion
##          of the distribution function
## Input:
## p = the probability at which the quantile function is evaluated
## initDist = the initial distribution
## T.mat = the subintensity matrix
##
## Output:
## The quantile function at p
##----------------------------------------------------
qcontphasetype <- function(p, initDist, T.mat){
  
  m <- uniroot(function(y) pcontphasetype(x = y, initDist = initDist, T.mat = T.mat)-p, c(0, 400))
  return(m$root[1])
  
}

##----------------------------------------------------
## The quantile function for of a discrete phase-type distribution
##----------------------------------------------------
## Name: qdiscphasetype
## Purpose: Computing the quantile function of a discrete
##          phase-type distribution by numeric inversion
##          of the distribution function
## Input:
## p = the probability at which the quantile function is evaluated
## initDist = the initial distribution
## T.mat = the subtransition matrix
##
## Output:
## The quantile function at p
##----------------------------------------------------
qdiscphasetype <- function(p, initDist, T.mat){
  
  m <- uniroot(function(y) pdiscphasetype(x = y, initDist = initDist, T.mat = T.mat)-p, c(0, 400))
  return(round(m$root[1]))
  
}

qphasetype <- function(...){
  UseMethod("qphasetype")
}

##----------------------------------------------------
## The quantile function for of a continuous phase-type distribution
##----------------------------------------------------
## Name: qcontphasetype
## Purpose: Computing the quantile function of a continuous
##          phase-type distribution by numeric inversion
##          of the distribution function
## Input:
## cptd = The continuous phase-type distribution object
## p = the probability at which the quantile function is evaluated
##
## Output:
## The quantile function at p
##----------------------------------------------------
qphasetype.contphasetype <- function(cptd, p){
  
  qcontphasetype(p = p, initDist = cptd$initDist, T.mat = cptd$T.mat)
  
}

##----------------------------------------------------
## The quantile function for of a discrete phase-type distribution
##----------------------------------------------------
## Name: qcontphasetype
## Purpose: Computing the quantile function of a discrete
##          phase-type distribution by numeric inversion
##          of the distribution function
## Input:
## cptd = The discrete phase-type distribution object
## p = the probability at which the quantile function is evaluated
##
## Output:
## The quantile function at p
##----------------------------------------------------
qphasetype.discphasetype <- function(cptd, p){
  
  qdiscphasetype(p = p, initDist = cptd$initDist, T.mat = cptd$T.mat)
  
}



##----------------------------------------------------
## Simulation from a discrete phase-type distribution
##----------------------------------------------------
## Name: rdiscphasetype
## Purpose: Generating outcomes from a discrete
##          phase-type distribution
## 
## Input:
## n = the number of outcomes
## initDist = the initial distribution
## T.mat = the subtransition matrix
##
## Output:
## A vector of length n with outcomes from the discrete
## phase-type distribution
##----------------------------------------------------
rdiscphasetype <- function(n, initDist, T.mat){
  
  # Calculate the number of transient states
  p <- length(initDist)
  
  # Calculate the defect
  defect <- 1 - sum(initDist)
  
  # Initialize the vector with 1's because if the Markov Chain is absorbed immediately then tau would be 1
  tau <- rep(1,n)
  
  # If the defect is positive immediate absorption needs to be a possibility
  # for that we generate n uniform(0,1) variables
  if(defect < 1){
    u <- runif(n)
  }
  
  # We make the rows of the transition matrix corresponding to the transient states from
  # the subtransition matrix
  TransMat <- cbind(T.mat, 1-rowSums(T.mat))
  
  # Now we simulate the Markov Chain until absorption in 'p+1'
  for(i in 1:n){
    if(u[i] <= defect){
      next()
    }
    else{
     initState <- sample(p, size = 1, prob = initDist)
     x <- initState
     while(x != p + 1){
       
       tau[i] <- tau[i] + 1
       x <- sample(p + 1, size = 1, prob = TransMat[x,] )
       
     }
    }
  }
  return(tau)
}

##----------------------------------------------------
## Simulation from a continuous phase-type distribution
##----------------------------------------------------
## Name: rcontphasetype
## Purpose: Generating outcomes from a continuous
##          phase-type distribution
## 
## Input:
## n = the number of outcomes
## initDist = the initial distribution
## T.mat = the subintensity matrix
##
## Output:
## A vector of length n with outcomes from the continuous
## phase-type distribution
##----------------------------------------------------
rcontphasetype <- function(n, initDist, T.mat){
  
  # Calculate the number of transient states
  p <- length(initDist)
  
  # Calculate the defect
  defect <- 1 - sum(initDist)
  
  # Initialize the vector with 0's because if the Markov Chain is absorbed immediately
  # then tau would be 0
  tau <- rep(0,n)
  
  # If the defect is positive immediate absorption needs to be a possibility
  # for that we generate n uniform(0,1) variables
  if(defect < 1){
    u <- runif(n)
  }
  
  # We make the rows of the intensity matrix corresponding to the transient states from
  # the subintensity matrix
  IntenseMat <- cbind(T.mat, 1-rowSums(T.mat))
  
  # Now we simulate the Markov Jump Process until absorption in 'p+1'
  for(i in 1:n){
    if(u[i] <= defect){
      next()
    }
    else{
      initState <- sample(p, size = 1, prob = initDist)
      x <- initState
      while(x != p + 1){
        
        tau[i] <- tau[i] + rexp(1, rate = -IntenseMat[x,x])
        x <- sample((1:(p+1))[-x], size = 1, prob = IntenseMat[x,-x] )
        
      }
    }
  }
  return(tau)
}

