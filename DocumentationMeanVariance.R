#' Mean and variance of phase-type distributions
#'
#' Mean and variance for both the discrete and continuous 
#' phase-type distribution with initial distribution equal to
#' \code{initDist} and subtransition/subintensity matrix equal to
#' \code{P.mat}/\code{T.mat}.
#'
#' In the discrete case, the phase-type distribution has mean
#' \deqn{E[\tau] = initDist (I-P.mat)^{-1} e + 1 - initDist e},
#' where \code{initDist} is the initial distribution, \code{P.mat} is the subtransition
#' probability matrix and \code{e} is the vector having one in each entry. 
#' Furthermore, the variance can be calculated as
#' \deqn{Var[\tau] = E[\tau(\tau-1)] + E[\tau] - E[\tau]^2}, 
#' where 
#' \deqn{E[\tau(\tau-1)] = 2 initDist P.mat (I-P.mat)^{-2} e + 1 - initDist e}.
#' In the continuous case, the phase-type distribution has mean
#' \deqn{E[\tau] = initDist (-T.mat)^{-1} e},
#' where \code{initDist} is the initial distribution and \code{T.mat} is the subintensity
#' rate matrix. Furthermore, the variance can be calculated in the usual way
#' \deqn{Var[tau] = E[tau^2] - E[tau]^2},
#' where 
#' \deqn{E[\tau^2] = 2 initDist (-T.mat)^{-2} e}.
#'
#' @param object an object for which the mean or variance should be computed.
#' To be able to use these function,the object has to be of
#' class \code{discphasetype} or \code{contphasetype}.
#'
#' @return \code{mean} gives the mean and \code{var} gives the 
#' variance of the phase-type distribution. The length of the output is 1.
#'
#' @source Mogens Bladt and Bo Friis Nielsen (2017):
#' \emph{ Matrix-Exponential Distributions in Applied Probability}.
#' Probability Theory and Stochastic Modelling (Springer), Volume 81.
#'
#' @seealso \code{\link{mean}}, \code{\link{var}}.
#'
#' @examples
#'
#'
mean.discphasetype <- function(object){
  
  initDist <- object$initDist
  P.mat <- object$P.mat
  
  return(sum(initDist%*%solve(diag(x=1, nrow = nrow(P.mat))-P.mat)) + 1 - sum(initDist))
  
}

#' @rdname mean.discphasetype
mean.contphasetype <- function(object){
  
  initDist = object$initDist
  T.mat= object$T.mat
  
  return(sum(initDist%*%solve(-T.mat)))
}

#' @rdname mean.discphasetype
var <- function(...){
  
  UseMethod("var")
}

#' @rdname mean.discphasetype
var.default <- function(x,...){
  
  var(x,...)
}

#' @rdname mean.discphasetype
var.discphasetype <- function(object){
  
  initDist = object$initDist
  P.mat = object$P.mat
  defect <- 1 - sum(initDist)
  
  secondMoment <- 2*sum(initDist%*%P.mat%*%solve((diag(x=1, nrow = nrow(P.mat))-P.mat)%^%2))
  firstmoment <- sum(initDist%*%solve(diag(x=1, nrow = nrow(P.mat))-P.mat)) + defect
  return(secondMoment + firstmoment - firstmoment^2)
}

#' @rdname mean.discphasetype
var.contphasetype <- function(object){
  
  initDist = object$initDist
  T.mat = object$T.mat
  
  secondMoment <- 2*sum(initDist%*%solve(T.mat%^%2))
  firstmoment <- sum(initDist%*%solve(-T.mat))
 
  return(secondMoment-firstmoment^2) 
}
