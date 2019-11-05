phasetypexamples <- list(
"discphasetype1" = discphasetype(initDist = c(.5,.2,.3,0),
                                P.mat = matrix(c(.1,.3,.2,.6,.3,.2,.1,.2,.1,.1,.2,.1,.2,.1,.3,.1), nrow = 4)),
"discphasetype2" = discphasetype(initDist = c(.2,.7), 
                                  P.mat = matrix(c(.5,0,.3,.2), nrow = 2)),
"contphasetype1" = contphasetype(initDist = c(.5,.2,.3,0),
                                 T.mat = matrix(c(-6,1,1,0,2,-8,0,0,2,3,-5,0,0,3,2,-4), nrow = 4)),
"contphasetype2" = contphasetype(initDist = c(.1,.1,.3),
                                T.mat = matrix(c(-2,0,.5,1,-3,1,.25,2,-4), nrow = 3)),
"contphasetype3" = contphasetype(initDist = c(.2,.7), 
                                T.mat = matrix(c(-7,3,2,-4), nrow = 2))
)

#The continuous phase-type distribution corresponding to convolution of n exponential distributions
#with parameters contained in the parm-vector
expconv <- function(parm){
  n <- length(parm)
  if(n == 1){
    return(contphasetype(1,-parm))
  }else{
    T.mat <- -diag(parm)
    if(n == 2){
      T.mat[1,2] <- parm[1]
    }else{
    diag(T.mat[(1:(n-1)),2:n]) <- parm[1:(n-1)]
    }
    return(contphasetype(initDist = c(1,rep(0,n-1)), T.mat = T.mat))
  }
}

#TMRCA for sample size n
TMRCA <- function(n){
  expconv(choose(n:2,2))
}

#TTotal for sample size n
TTotal <- function(n){
  expconv((n:2-1)/2)
}

#STotal
STotal <- function(n, mutrate){
  discretization(object = TTotal, mutrate = mutrate)
}

# Example the contphasetype corresponding to the total length of double-ton branches for sample size = 6
doubletonbranches6 <- contphasetype(initDist = c(1,0,0,0,0), 
                             T.mat= matrix(c(-10,0,0,0,0,6,-3,0,0,0,2,2,-3,0,0,0,.5,0,-1,0,2/3,1/6,1,1,-1),
                                           nrow = 5))
# Via discretization we can obtain the discrete phase-type parameter for the number of
# doubletons + 1
doubletons6 <- discretization(doubletonbranches6, mutrate = 2)
