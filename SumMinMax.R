##----------------------------------------------------
## The sum of discrete phase-type distributions
##----------------------------------------------------
## Name: sum.discphasetype
## Purpose: The phase-type parameters for a sum of 
##          independent discrete phase-type distributed 
##          variables.(using Theorem 1.2.65 in [BN])
##
## Input:
## object1,object2 = the discphasetype objects
##
## Output:
## A discphasetype object with the distribution
## of the sum of the input objects
##----------------------------------------------------
sum.discphasetype <- function(object1,object2){
  
  initDist1 = object1$initDist
  P.mat1 = object1$P.mat
  initDist2 = object2$initDist
  P.mat2 = object2$P.mat
  
  newInitDist <- c(initDist1,rep(0,length(initDist2)))
  newP.mat <- rbind(
    cbind(P.mat1,(1-rowSums(P.mat1))%*%t(initDist2)),
    cbind(matrix(0, nrow = length(initDist2), ncol = length(initDist2)), P.mat2))
  
  discphasetype(initDist = newInitDist, P.mat = newP.mat)
  
}

##----------------------------------------------------
## The minimum of discrete phase-type distributions
##----------------------------------------------------
## Name: min.discphasetype
## Purpose: The phase-type parameters for a minimum of 
##          independent discrete phase-type distributed 
##          variables.(using Theorem 1.2.67 in [BN])
##
## Input:
## object1,object2 = the discphasetype objects
##
## Output:
## A discphasetype object with the distribution
## of the minimum of the input objects
##----------------------------------------------------
min.discphasetype <- function(object1,object2){
  
  initDist1 = object1$initDist
  P.mat1 = object1$P.mat
  initDist2 = object2$initDist
  P.mat2 = object2$P.mat

  newInitDist <- kronecker(initDist1,initDist2)
  newP.mat <- kronecker(P.mat1, P.mat2)
  
  discphasetype(initDist = newInitDist, P.mat = newP.mat)
  
}

##----------------------------------------------------
## The maximum of discrete phase-type distributions
##----------------------------------------------------
## Name: max.discphasetype
## Purpose: The phase-type parameters for a maximum of 
##          independent discrete phase-type distributed 
##          variables.(using Theorem 1.2.67 in [BN])
##
## Input:
## object1,object2 = the discphasetype objects
##
## Output:
## A discphasetype object with the distribution
## of the maximum of the input objects
##----------------------------------------------------
max.discphasetype <- function(object1,object2){
  
  initDist1 = object1$initDist
  P.mat1 = object1$P.mat
  p <- length(initDist1)
  initDist2 = object2$initDist
  P.mat2 = object2$P.mat
  q <- length(initDist2)
  
  newInitDist <- c(kronecker(initDist1,initDist2),rep(0,p+q))
  newP.mat <- rbind(cbind(kronecker(P.mat1,P.mat2),
                          kronecker(P.mat1,1-rowSums(P.mat2)),
                          kronecker(1-rowSums(P.mat1),P.mat2)),
                    cbind(matrix(0, ncol = p*q, nrow = p), P.mat1, matrix(0, ncol = q, nrow = p)),
                    cbind(matrix(0, ncol = p*q, nrow = q), matrix(0, ncol = p, nrow = q), P.mat2))
  
  discphasetype(initDist = newInitDist, P.mat = newP.mat)
  
}

##----------------------------------------------------
## The sum of continuous phase-type distributions
##----------------------------------------------------
## Name: sum.contphasetype
## Purpose: The phase-type parameters for a sum of 
##          independent continuous phase-type distributed 
##          variables.(using Theorem 3.1.26 in [BN])
##
## Input:
## object1,object2 = the contphasetype objects
##
## Output:
## A contphasetype object with the distribution
## of the sum of the input objects
##----------------------------------------------------
sum.contphasetype <- function(object1,object2){
  
  initDist1 = object1$initDist
  T.mat1 = object1$T.mat
  initDist2 = object2$initDist
  T.mat2 = object2$T.mat
  
  newInitDist <- c(initDist1,rep(0,length(initDist2)))
  newT.mat <- rbind(
    cbind(T.mat1, rowSums(-T.mat1) %*% t(initDist2)),
    cbind(matrix(0, nrow = length(initDist2), ncol = length(initDist2)), T.mat2))
  
  contphasetype(initDist = newInitDist, T.mat = newT.mat)
  
}

##----------------------------------------------------
## The minimum of continuous type distributions
##----------------------------------------------------
## Name: min.contphasetype
## Purpose: The phase-type parameters for a minimum of 
##          independent continuous phase-type distributed 
##          variables.(using Corollary 3.1.32 in [BN])
##
## Input:
## object1,object2 = the contphasetype objects
##
## Output:
## A contphasetype object with the distribution
## of the minimum of the input objects
##----------------------------------------------------
min.contphasetype <- function(object1,object2){
  
  initDist1 = object1$initDist
  T.mat1 = object1$T.mat
  initDist2 = object2$initDist
  T.mat2 = object2$T.mat
  
  newInitDist <- kronecker(initDist1,initDist2)
  newT.mat <- kronecker(T.mat1, diag(nrow = nrow(T.mat2))) + 
              kronecker(diag(nrow = nrow(T.mat1)), T.mat2)
  
  contphasetype(initDist = newInitDist, T.mat = newT.mat)
  
}

##----------------------------------------------------
## The maximum of continuous type distributions
##----------------------------------------------------
## Name: max.contphasetype
## Purpose: The phase-type parameters for a maximum of 
##          independent continuous phase-type distributed 
##          variables.(using Corollary 3.1.32 in [BN])
##
## Input:
## object1,object2 = the contphasetype objects
##
## Output:
## A contphasetype object with the distribution
## of the maximum of the input objects
##----------------------------------------------------
max.contphasetype <- function(object1,object2){
  
  initDist1 = object1$initDist
  T.mat1 = object1$T.mat
  p <- length(initDist1)
  initDist2 = object2$initDist
  T.mat2 = object2$T.mat
  q <- length(initDist2)
  
  newInitDist <- c(kronecker(initDist1,initDist2),rep(0,p+q))
  firstblockrow <- cbind(kronecker(T.mat1, diag(nrow = nrow(T.mat2))) +
                           kronecker(diag(nrow = nrow(T.mat1)), T.mat2),
                         kronecker(diag(nrow = nrow(T.mat1)),rowSums(-T.mat2)),
                         kronecker(rowSums(-T.mat1),diag(nrow = nrow(T.mat2))))
  
  newT.mat <- rbind(firstblockrow,
                    cbind(matrix(0, ncol = p*q, nrow = p), T.mat1, matrix(0, ncol = q, nrow = p)),
                    cbind(matrix(0, ncol = p*q, nrow = q), matrix(0, ncol = p, nrow = q), T.mat2))
  
  contphasetype(initDist = newInitDist, T.mat = newT.mat)
  
}
