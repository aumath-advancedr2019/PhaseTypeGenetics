#' Statistical moments of phase-type distributions
#'
#' Computing (factorial) moments of a given order of phase-type distributed random variables.
#' 
#' In the discrete case ( tau ~ DPH(pi,P) ), the factorial moments are given by
#' E[tau(tau-1)*...*(tau-i+1)] = i! * pi %*% P^(i-1) %*% (I-P)^(-i)%*% e.
#' For tau ~ PH(pi, T), the i'th-order moment is defined as
#' E[tau^i] = i!* pi %*% (-T)^(-i)%*%e.
#' In both cases, \code{e} is a vector with one in each entry. 
#' 
#' 
#' @param object either a continuous phase-type distributed object of class \code{contphasetype} or 
#' a discrete phase-type distributed object of class \code{discphasetype}.
#' @param i a positive number stating the order of the desired moment
#' @param all a logical value indicating whether the function should compute
#' all moments up to the given order. The default is equal to FALSE.
#' 
#' @return For \code{all = FALSE}, the function either returns the i'th-order moment (if the object is continuous 
#' phase-type distributed) or the i'th factorial moment (if the object is discrete
#' phase-type distributed). In both cases, the length of the output is one. 
#' For \code{all = TRUE}, the function computes all (factorial) moments up to the given order,
#' hence the output is a vector of length \code{i}.
#'
#' @source Mogens Bladt and Bo Friis Nielsen (2017):
#' \emph{ Matrix-Exponential Distributions in Applied Probability}.
#' Probability Theory and Stochastic Modelling (Springer), Volume 81.
#' 
#'
#' @examples
#'
#'
#' @export

#' @describeIn moments generic
moments <- function(..., all = FALSE){
  
  UseMethod("moments")
}

#' @describeIn momets the default moment function
moments.default <- function(...){
  
  print("The package ‘moments’ can be used to compute sample moments of specified order.")
}

#' @describeIn momets computing moments of a discrete phase-type distribution
moments.discphasetype <- function(object, i, all){
  
  if(i < 1 ) stop("Not a valid order. The number i has to be positive!")
  
  initDist = object$initDist
  P.mat = object$P.mat
  
  res <- NULL
  if(all == FALSE){
    
    res <- factorial(i)*sum(initDist%*%P.mat^(i-1)%*%
           (solve(diag(1, nrow = nrow(P.mat))-P.mat)%^%(i)))
  }else{
    for (k in 1:i) {
      
      res[k] <- factorial(k)*sum(initDist%*%P.mat^(k-1)%*%
                (solve(diag(1, nrow = nrow(P.mat))-P.mat)%^%(k)))
    }
    names(res) <- 1:i
  }
  return(res)
}

#' @describeIn momets computing moments of a continuous phase-type distribution
moments.contphasetype <- function(object, i, all){
  
  if(i < 1 ) stop("Not a valid order. The number i has to be positive!")
  
  initDist = object$initDist
  T.mat = object$T.mat
  
  res <- NULL
  if(all == FALSE){
    
    res <- factorial(i)*sum(initDist%*%(solve(-T.mat)%^%(i)))
  }else{
    for (k in 1:i) {
      
      res[k] <- factorial(k)*sum(initDist%*%(solve(-T.mat)%^%(k)))
    }
    names(res) <- 1:i
  }
  return(res)
}