#' Sum of phase-type distributions
#'
#' The sum of two independent discrete or continuous 
#' phase-type distributed variables with initial distributions
#' \code{initDist1} and \code{initDist2} as well as subtransition/subintensity 
#' matrices equal to \code{P.mat1}/\code{T.mat1} and  \code{P.mat2}/\code{T.mat2}.
#'
#' In the discrete case, the sum of two phase-type distributed 
#' variables \eqn{tau1 ~ DPH_p(\alpha,S)} and \eqn{tau2 ~ DPH_q(\beta,T)} is 
#' again discrete phase-type distributed in the following way
#' \deqn{tau1 + tau2 ~ DPH_{p+q}((\alpha,0),cbind((S, s \beta),(0,T)) )}.
#' In the continous case, the sum of two phase-type distributed 
#' variables \eqn{X ~ PH_p(\alpha,S)} and \eqn{Y ~ PH_q(\beta,T)} is 
#' again continuous and phase-type distributed in the following way
#' \deqn{X + Y ~ PH_{p+q}((\alpha,0),cbind((S, s \beta),(0,T)) ).}
#' 
#' @param object1,object2 two objects of class \code{discphasetype} 
#' or \code{contphasetype} for which the sum should be computed.
#'
#' @return The function \code{sum} returns an object of type \code{discphasetype} 
#' or \code{contphasetype} (depending on the input) holding the ditribution 
#' of the sum of the input objects.
#'
#' @source Mogens Bladt and Bo Friis Nielsen (2017):
#' \emph{ Matrix-Exponential Distributions in Applied Probability}.
#' Probability Theory and Stochastic Modelling (Springer), Volume 81.
#'
#' @seealso \code{\link{sum}}.
#'
#' @examples
#'
#'
#' @describeIn sum computing the sum of two discrete phase-type distributions
sum.discphasetype <- function(object1,object2){
  
  initDist1 = object1$initDist
  P.mat1 = object1$P.mat
  initDist2 = object2$initDist
  P.mat2 = object2$P.mat
  
  newInitDist <- c(initDist1,rep(0,length(initDist2)))
  newP.mat <- rbind(
    cbind(P.mat1,(1-rowSums(P.mat1))%*%t(initDist2)),
    cbind(matrix(0, nrow = length(initDist2), ncol = length(initDist2)), P.mat2))
  
  return(discphasetype(initDist = newInitDist, P.mat = newP.mat))
  
}

#' @describeIn sum computing the sum of two continuous phase-type distributions
sum.contphasetype <- function(object1,object2){
  
  initDist1 = object1$initDist
  T.mat1 = object1$T.mat
  initDist2 = object2$initDist
  T.mat2 = object2$T.mat
  
  newInitDist <- c(initDist1,rep(0,length(initDist2)))
  newT.mat <- rbind(
    cbind(T.mat1, rowSums(-T.mat1) %*% t(initDist2)),
    cbind(matrix(0, nrow = length(initDist2), ncol = length(initDist2)), T.mat2))
  
  return(contphasetype(initDist = newInitDist, T.mat = newT.mat))
  
}
