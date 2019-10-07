##----------------------------------------------------
## The mean of a continuous phase-type distribution
##----------------------------------------------------
## Name: mean_contphasetype
## Purpose: Computing the mean of a continuous phase-
##          type distribution (Corollary 1.2.64 in [BN])
## Input:
## initDist = the initial distribution
## T.mat = the subintensity matrix from the
##        continuous phase type distribution
## Output: 
## The mean E[tau] 
##----------------------------------------------------

mean_contphasetype <- function(initDist, T.mat){
  
  return(sum(initDist%*%solve(-T.mat)))
  
}

##----------------------------------------------------
## The mean of a continuous phase-type distribution
##----------------------------------------------------
## Name: mean.contphasetype
## Purpose: Computing the mean of a continuous phase-
##          type distribution (Corollary 1.2.64 in [BN])
## Input:
## cptd = the continuous phase-type distribution object 
## Output: 
## The mean E[tau] 
##----------------------------------------------------

mean.contphasetype <- function(cptd){
  
  mean_contphasetype(initDist = cptd$initDist, T.mat= cptd$T.mat)
  
}

##----------------------------------------------------
## The mean of a discrete phase-type distribution
##----------------------------------------------------
## Name: mean_discphasetype
## Purpose: Computing the mean of a discrete phase-
##          type distribution (Corollary 1.2.64 in [BN])
## Input:
## T.mat = the subtransition matrix from the
##        discrete phase type distribution
## initDist = the initial distribution
## Output: 
## The mean E[tau] 
##----------------------------------------------------

mean_discphasetype <- function(dptd){
  
  return(sum(initDist%*%solve(diag(x=1, nrow = nrow(T.mat))-T.mat)) + 1 - sum(initDist))
  
}

##----------------------------------------------------
## The mean of a discrete phase-type distribution
##----------------------------------------------------
## Name: mean.discphasetype
## Purpose: Computing the mean of a discrete phase-
##          type distribution (Corollary 1.2.64 in [BN])
## Input:
## dptd = the discrete phase-type distribution object 
## Output: 
## The mean E[tau] 
##----------------------------------------------------

mean.discphasetype <- function(dptd){
  mean_discphasetype(initDist = dptd$initDist, T.mat=dptd$T.mat)
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
## Name: var.discphasetype
## Purpose: Computing the variance of a discrete phase-
##          type distribution 
##          (using Theorem 1.2.69 in [BN])
## Input:
## dptd = the discrete phase-type distribution object which
## is a list with two entries
## dptd$T.mat = the subtransition matrix from the
##        discrete phase type distribution
## dptd$initDist = the initial distribution
## Output: 
## The variance Var[tau] 
##----------------------------------------------------
var.discphasetype <- function(dptd){
  defect <- 1 - sum(dptd$initDist)
  secondMoment <- 2*sum(dptd$initDist%*%dptd$T.mat%*%solve((diag(x=1, nrow = nrow(dptd$T.mat))-dptd$T.mat)%^%2))
  firstmoment <- sum(dptd$initDist%*%solve(diag(x=1, nrow = nrow(dptd$T.mat))-dptd$T.mat)) + defect
  return(secondMoment + firstmoment - firstmoment^2)
}

##----------------------------------------------------
## The variance of a continuous phase-type distribution
##----------------------------------------------------
## Name: var.contphasetype
## Purpose: Computing the variance of a continuous 
##          phase-type distribution 
##          (using Theorem 1.2.69 in [BN])
## Input:
## cptd = the continuous phase-type distribution object which
## is a list with two entries
## cptd$T.mat = the subintensity matrix from the
##        continuous phase type distribution
## cptd$initDist = the initial distribution
## Output: 
## The variance Var[tau] 
##----------------------------------------------------
var.contphasetype <- function(cptd){
  
  return(LaplacePhaseType(initDist = cptd$initDist, Tmat = cptd$T.mat, i=2)-
    LaplacePhaseType(initDist = cptd$initDist, Tmat = cptd$T.mat, i=1)^2)
}

