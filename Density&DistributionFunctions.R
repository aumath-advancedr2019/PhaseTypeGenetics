##----------------------------------------------------
## The density function for of a continuous phase-type distribution
##----------------------------------------------------
## Name: dcontphasetype
## Purpose: Computing the density function of a continuous
##          phase-type distribution
##          (using Theorem 3.1.7 in [BN])
## Input:
## x = the number at which the density is evaluated
## initDist = the initial distribution
## T.mat = the subintensity matrix
##
## Output:
## The density function at x, f_tau(x)
##----------------------------------------------------
dcontphasetype <- function(x, initDist, T.mat){
  
    return(-sum(initDist %*% expm(x * T.mat) %*% T.mat))
  
}

##----------------------------------------------------
## The density function for of a discrete phase-type distribution
##----------------------------------------------------
## Name: ddiscphasetype
## Purpose: Computing the density function of a discrete
##          phase-type distribution
##          (using Theorem 1.2.58 in [BN])
## Input:
## x = the number at which the density is evaluated
## initDist = the initial distribution
## T.mat = the subtransition matrix
##
## Output:
## The density function at x, f_tau(x)
##----------------------------------------------------
ddiscphasetype <- function(x, initDist, T.mat){
  sum(initDist%*%(T.mat %^% x)%*%(diag(1, nrow = nrow(T.mat))-T.mat))
}


# The generic function
dphasetype <- function(...){
  
  UseMethod("dphasetype")
  
}

##----------------------------------------------------
## The density function for of a continuous phase-type distribution
##----------------------------------------------------
## Name: dphasetype.contphasetype
## Purpose: Computing the density function of a continuous
##          phase-type distribution
##          (using Theorem 3.1.7 in [BN])
## Input:
## cptd = The continuous phase-type distribution object
## x = the number at which the density is evaluated
##
## Output:
## The density function at x, f_tau(x)
##----------------------------------------------------
dphasetype.contphasetype <- function(cptd, x){
  
  dcontphasetype(x = x, initDist = cptd$initDist, T.mat = cptd$T.mat)
  
}

##----------------------------------------------------
## The density function for of a discrete phase-type distribution
##----------------------------------------------------
## Name: dphasetype.discphasetype
## Purpose: Computing the density function of a continuous
##          phase-type distribution
##          (using Theorem 1.2.58 in [BN])
## Input:
## dptd = the discrete phase-type distribution object
## x = the number at which the density is evaluated
##
## Output:
## The density function at x, f_tau(x)
##----------------------------------------------------
dphasetype.discphasetype <- function(dptd,x){
  
  ddiscphasetype(x = x, initDist = dptd$initDist, T.mat = dptd$T.mat)
  
}

##----------------------------------------------------
## The distribution function for of a continuous phase-type distribution
##----------------------------------------------------
## Name: pcontphasetype
## Purpose: Computing the distribution function of a continuous
##          phase-type distribution
##          (using Theorem 3.1.8 in [BN])
## Input:
## x = the number at which the distribution is evaluated
## initDist = the initial distribution
## T.mat = the subintensity matrix
##
## Output:
## The distribution function at x, F_tau(x)
##----------------------------------------------------
pcontphasetype <- function(x, initDist, T.mat){
  
  return(1 - sum(initDist %*% expm(x * T.mat)))
  
}

##----------------------------------------------------
## The distribution function for of a discrete phase-type distribution
##----------------------------------------------------
## Name: pdiscphasetype
## Purpose: Computing the distribution function of a discrete
##          phase-type distribution
##          (using Theorem 1.2.59 in [BN])
## Input:
## x = the number at which the distribution is evaluated
## initDist = the initial distribution
## T.mat = the subtransition matrix
##
## Output:
## The distribution function at x, F_tau(x)
##----------------------------------------------------
pdiscphasetype <- function(x, initDist, T.mat){
  return(1 - sum(initDist%*%(T.mat %^% x)))
}


# The generic function
pphasetype <- function(...){
  
  UseMethod("pphasetype")
  
}

##----------------------------------------------------
## The distribution function for of a continuous phase-type distribution
##----------------------------------------------------
## Name: dphasetype.contphasetype
## Purpose: Computing the distribution function of a continuous
##          phase-type distribution
##          (using Theorem 1.2.59 in [BN])
## Input:
## cptd = The continuous phase-type distribution object
## x = the number at which the distribution is evaluated
##
## Output:
## The distribution function at x, F_tau(x)
##----------------------------------------------------
pphasetype.contphasetype <- function(cptd, x){
  
  pcontphasetype(x = x, initDist = cptd$initDist, T.mat = cptd$T.mat)
  
}

##----------------------------------------------------
## The distribution function for of a discrete phase-type distribution
##----------------------------------------------------
## Name: dphasetype.discphasetype
## Purpose: Computing the distribution function of a continuous
##          phase-type distribution
##          (using Theorem 3.1.8 in [BN])
## Input:
## dptd = the discrete phase-type distribution object
## x = the number at which the distribution is evaluated
##
## Output:
## The distribution function at x, F_tau(x)
##----------------------------------------------------
pphasetype.discphasetype <- function(dptd,x){
  pdiscphasetype(x = x, initDist = dptd$initDist, T.mat = dptd$T.mat)
}



