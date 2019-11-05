#' Minima and maxima of phase-type distributions
#'
#' Computes the minima and maxima of two independent discrete or continuous 
#' phase-type distributed variables with initial distributions
#' \code{initDist1} and \code{initDist2} as well as subtransition/subintensity 
#' matrices equal to \code{P.mat1}/\code{T.mat1} and  \code{P.mat2}/\code{T.mat2}.
#'
#' In the discrete case, the minimum and maximum of two phase-type distributed 
#' variables \eqn{tau1 ~ DPH_p(\alpha,S)} and \eqn{tau2 ~ DPH_q(\beta,T)} are defined 
#' as follows
#' \deqn{min(tau1, tau2) ~ DPH_{pq}( kronecker(\alpha,\beta), kronecker(S,T) )}
#' and 
#' \deqn{max(tau1, tau2) ~ DPH_{pq + p + q}( c(kronecker(\alpha,\beta),0.vec), K )}, 
#' where 
#' \code{0.vec} is a vector of length \eqn{p+q} with zero in each entry and 
#' \deqn{K= rbind( cbind(kronecker(S,T), kronecker(S,t),kronecker(s,T) ),
#'           cbind(matrix(0, ncol = p*q, nrow = p), S , matrix(0, ncol = q, nrow = p)), 
#'           cbind(matrix(0, ncol = p*q, nrow = q), matrix(0, ncol = p, nrow = q), T)).}
#' 
#' In the continous case, the minima and maxima of two phase-type distributed 
#' variables \eqn{X ~ PH_p(\alpha,S)} and \eqn{Y ~ PH_q(\beta,T)} is 
#' given in the following way
#' \deqn{min(X, Y) ~ PH( kronecker(\alpha,\beta), kronecker(S,T) )}
#' and 
#' \deqn{max(X, Y) ~ PH( c(kronecker(alpha,beta),0.vec), K )}
#' where 
#' \code{0.vec} is a vector of length \eqn{p+q} with zero in each entry and 
#' \deqn{K= rbind( cbind(kronecker(S,T), kronecker(I,t),kronecker(s,I) ),
#'           cbind(matrix(0, ncol = p*q, nrow = p), S , matrix(0, ncol = q, nrow = p)), 
#'           cbind(matrix(0, ncol = p*q, nrow = q), matrix(0, ncol = p, nrow = q), T)).}
#' 
#' 
#' @param object1,object2 two objects of class \code{discphasetype} 
#' or \code{contphasetype} for which the maximum or minimum should be computed.
#'
#' @return The function \code{minima} returns an object of type \code{discphasetype} 
#' or \code{contphasetype} (depending on the input) holding the ditribution 
#' of the minimum of the input objects, while \code{maxima} returns an object of type \code{discphasetype} 
#' or \code{contphasetype} that holds the ditribution 
#' of the maximum of the input objects.
#'
#' @source Mogens Bladt and Bo Friis Nielsen (2017):
#' \emph{ Matrix-Exponential Distributions in Applied Probability}.
#' Probability Theory and Stochastic Modelling (Springer), Volume 81.
#'
#' @seealso \code{\link{min}} or \code{\link{max}}.
#'
#' @examples
#'
#' ## Two representations of the total branch length
#' ## are given by
#' T_Total1 <- matrix(c(-1.5, 1.5, 0,
#'                      0, -1, 1,
#'                      0, 0, -0.5), nrow = 3, byrow = TRUE)
#'                      
#' T_Total2 <- matrix(c(-1.5, 1.5, 0, 0,
#'                      0, -1, 2/3, 1/3,
#'                      0, 0, -0.5, 0,
#'                      0, 0, 0, -0.5), nrow = 4, byrow = TRUE)
#'                      
#' ## Defining two objects of type "contphasetype".
#' T_Total1 <- contphasetype(initDist = c(1,0,0), T_Total1)
#' T_Total2 <- contphasetype(initDist = c(1,0,0,0), T_Total2)
#' 
#' ## Computing the minimum and maximum
#' minima(T_Total1, T_Total2)
#' maxima(T_Total1, T_Total2)
minima <- function(...) UseMethod("minima")

#' @rdname minima
minima.default <- function(...){
  
  min(...)
}

#' @rdname minima
minima.discphasetype <- function(object1,object2){
  
  initDist1 = object1$initDist
  P.mat1 = object1$P.mat
  initDist2 = object2$initDist
  P.mat2 = object2$P.mat
  
  newInitDist <- kronecker(initDist1,initDist2)
  newP.mat <- kronecker(P.mat1, P.mat2)
  
  return(discphasetype(initDist = newInitDist, P.mat = newP.mat))
  
}
#' @rdname minima
minima.contphasetype <- function(object1,object2){
  
  initDist1 = object1$initDist
  T.mat1 = object1$T.mat
  initDist2 = object2$initDist
  T.mat2 = object2$T.mat
  
  newInitDist <- kronecker(initDist1,initDist2)
  newT.mat <- kronecker(T.mat1, diag(nrow = nrow(T.mat2))) + 
    kronecker(diag(nrow = nrow(T.mat1)), T.mat2)
  
  return(contphasetype(initDist = newInitDist, T.mat = newT.mat))
  
}

#' @rdname minima
maxima <- function(...) UseMethod("maxima")

#' @rdname minima
maxima.default <- function(...){
  
  max(...)
}

#' @rdname minima
maxima.discphasetype <- function(object1,object2){
  
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
  
  return(discphasetype(initDist = newInitDist, P.mat = newP.mat))
}

#' @rdname minima
maxima.contphasetype <- function(object1,object2){
  
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
  
  return(contphasetype(initDist = newInitDist, T.mat = newT.mat))
}

