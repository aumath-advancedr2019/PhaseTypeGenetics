##----------------------------------------------------
## The mean of a continuous phase-type distribution
##----------------------------------------------------
## Name: mean.cphasetype
## Purpose: Computing the mean of a continuous phase-
##          type distribution (Corollary 1.2.64 in [BN])
## Input:
## T.mat = the subintensity matrix from the
##        continuous phase type distribution
## initDist = the initial distribution
## Output: 
## The mean E[tau] 
##----------------------------------------------------

mean.cphasetype <- function(initDist, T.mat){
  
  e.vec <- replicate(n=ncol(T.mat),1)
  return(initDist%*%solve(-T.mat)%*%e.vec)
  
}

##----------------------------------------------------
## The mean of a discrete phase-type distribution
##----------------------------------------------------
## Name: mean.dphasetype
## Purpose: Computing the mean of a discrete phase-
##          type distribution (Corollary 1.2.64 in [BN])
## Input:
## T.mat = the subintensity matrix from the
##        discrete phase type distribution
## pi.vec = the initial distribution
## Output: 
## The mean E[tau] 
##----------------------------------------------------

mean.dphasetype <- function(T.mat, pi.vec){
  
  e.vec <- replicate(n=ncol(T.mat),1)
  return(pi.vec%*%solve(diag(x=1, nrow = nrow(T.mat))-T.mat)%*%e.vec)
  
}

## Defining a generic var() function

var <- function(...){
  
  UseMethod("var")
}

## The default 
var.default <- function(x,...){
  
  var(x,...)
}


##----------------------------------------------------
## The variance of a discrete phase-type distribution
##----------------------------------------------------
## Name: var.dphasetype
## Purpose: Computing the variance of a discrete phase-
##          type distribution 
##          (using Theorem 1.2.69 in [BN])
## Input:
## T.mat = the subintensity matrix from the
##        discrete phase type distribution
## pi.vec = the initial distribution
## Output: 
## The variance Var[tau] 
##----------------------------------------------------
var.dphasetype <- function(T.mat,pi.vec){
  
  e.vec <- replicate(n=ncol(T.mat),1)
  secondMoment <- 2*pi.vec%*%T.mat%*%solve((diag(x=1, nrow = nrow(T.mat))-T.mat)%^%2)%*%e.vec
  firstmoment <- pi.vec%*%solve(diag(x=1, nrow = nrow(T.mat))-T.mat)%*%e.vec
  return(secondMoment +firstmoment - firstmoment^2)
}

##----------------------------------------------------
## The variance of a continuous phase-type distribution
##----------------------------------------------------
## Name: var.cphasetype
## Purpose: Computing the variance of a continuous 
##          phase-type distribution 
##          (using Theorem 1.2.69 in [BN])
## Input:
## T.mat = the subintensity matrix from the
##        discrete phase type distribution
## initDist = the initial distribution
## Output: 
## The variance Var[tau] 
##----------------------------------------------------
var.cphasetype <- function(T.mat, initDist){
  
  return(LaplacePhaseType(initDist = initDist, Tmat = T.mat, i=2)-
    LaplacePhaseType(initDist = initDist, Tmat = T.mat, i=1)^2)
}

