#' Minima and maxima of phase-type distributions
#'
#' Computes the minima and maxima of two independent discrete or continuous
#' phase-type distributed variables with initial distributions
#' \code{initDist1} and \code{initDist2} as well as sub-transition/sub-intensity
#' matrices equal to \code{P_Mat1}/\code{T_Mat1} and  \code{P_Mat2}/\code{T_Mat2}.
#'
#' In the discrete case, the minimum and maximum of two phase-type distributed
#' variables \eqn{tau1 ~ DPH_p(\alpha,S)} and \eqn{tau2 ~ DPH_q(\beta,T)} are defined
#' as follows
#' \deqn{min(tau1, tau2) ~ DPH_{pq}( kronecker(\alpha,\beta), kronecker(S,T) ),}
#' and
#' \deqn{max(tau1, tau2) ~ DPH_{pq + p + q}( c(kronecker(\alpha,\beta),0_vec), K ),}
#' where
#' \code{0_vec} is a vector of length \eqn{p+q} with zero in each entry and
#' \deqn{K= rbind( cbind(kronecker(S,T), kronecker(S,t),kronecker(s,T) ),
#'           cbind(matrix(0, ncol = p*q, nrow = p), S , matrix(0, ncol = q, nrow = p)),
#'           cbind(matrix(0, ncol = p*q, nrow = q), matrix(0, ncol = p, nrow = q), T)).}
#'
#' In the continuous case, the minima and maxima of two phase-type distributed
#' variables \eqn{X ~ PH_p(\alpha,S)} and \eqn{Y ~ PH_q(\beta,T)} is
#' given in the following way
#' \deqn{min(X, Y) ~ PH( kronecker(\alpha,\beta), kronecker(S,T) ),}
#' and
#' \deqn{max(X, Y) ~ PH( c(kronecker(alpha,beta),0_vec), K ),}
#' where
#' \code{0_vec} is a vector of length \eqn{p+q} with zero in each entry and
#' \deqn{K= rbind( cbind(kronecker(S,T), kronecker(I,t),kronecker(s,I) ),
#'           cbind(matrix(0, ncol = p*q, nrow = p), S , matrix(0, ncol = q, nrow = p)),
#'           cbind(matrix(0, ncol = p*q, nrow = q), matrix(0, ncol = p, nrow = q), T)).}
#'
#'
#' @param object1,object2 two objects of class \code{discphasetype}
#' or \code{contphasetype} for which the maximum or minimum should be computed.
#'
#' @return The function \code{minima} returns an object of type \code{discphasetype}
#' or \code{contphasetype} (depending on the input) holding the phase-type representation
#' of the minimum of the input objects, while \code{maxima} returns an object of type \code{discphasetype}
#' or \code{contphasetype} that holds the phase-type representation
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
#'
#' @export
minima <- function(object1,object2){

  UseMethod("minima")
}

#' @export
minima.default <- function(object1,object2){

  min(object1,object2)
}

#' @export
minima.discphasetype <- function(object1,object2){

  if(class(object1) != "discphasetype"| class(object2) != "discphasetype") stop("Invalid objects! object1 and object2 must be of class 'discphasetype'.")

  initDist1 = object1$initDist
  P_Mat1 = object1$P_Mat
  initDist2 = object2$initDist
  P_Mat2 = object2$P_Mat

  newInitDist <- kronecker(initDist1,initDist2)
  newP_Mat <- kronecker(P_Mat1, P_Mat2)

  return(discphasetype(initDist = newInitDist, P_Mat = newP_Mat))

}

#' @export
minima.contphasetype <- function(object1,object2){

  if(class(object1) != "contphasetype"| class(object2) != "contphasetype") stop("Invalid objects! object1 and object2 must be of class 'contphasetype'.")

  initDist1 = object1$initDist
  T_Mat1 = object1$T_Mat
  initDist2 = object2$initDist
  T_Mat2 = object2$T_Mat

  newInitDist <- kronecker(initDist1,initDist2)
  newT_Mat <- kronecker(T_Mat1, diag(nrow = nrow(T_Mat2))) +
    kronecker(diag(nrow = nrow(T_Mat1)), T_Mat2)

  return(contphasetype(initDist = newInitDist, T_Mat = newT_Mat))

}

#' @rdname minima
#' @export
maxima <- function(object1,object2){

  UseMethod("maxima")
}

#' @export
maxima.default <- function(object1,object2){

  max(object1,object2)
}

#' @export
maxima.discphasetype <- function(object1,object2){

  if(class(object1) != "discphasetype"| class(object2) != "discphasetype") stop("Invalid objects! object1 and object2 must be of class 'discphasetype'.")

  initDist1 = object1$initDist
  P_Mat1 = object1$P_Mat
  p <- length(initDist1)
  initDist2 = object2$initDist
  P_Mat2 = object2$P_Mat
  q <- length(initDist2)

  newInitDist <- c(kronecker(initDist1,initDist2),rep(0,p+q))
  newP_Mat <- rbind(cbind(kronecker(P_Mat1,P_Mat2),
                          kronecker(P_Mat1,1-rowSums(P_Mat2)),
                          kronecker(1-rowSums(P_Mat1),P_Mat2)),
                    cbind(matrix(0, ncol = p*q, nrow = p), P_Mat1, matrix(0, ncol = q, nrow = p)),
                    cbind(matrix(0, ncol = p*q, nrow = q), matrix(0, ncol = p, nrow = q), P_Mat2))

  return(discphasetype(initDist = newInitDist, P_Mat = newP_Mat))
}

#' @export
maxima.contphasetype <- function(object1,object2){

  if(class(object1) != "contphasetype"| class(object2) != "contphasetype") stop("Invalid objects! object1 and object2 must be of class 'contphasetype'.")

  initDist1 = object1$initDist
  T_Mat1 = object1$T_Mat
  p <- length(initDist1)
  initDist2 = object2$initDist
  T_Mat2 = object2$T_Mat
  q <- length(initDist2)

  newInitDist <- c(kronecker(initDist1,initDist2),rep(0,p+q))
  firstblockrow <- cbind(kronecker(T_Mat1, diag(nrow = nrow(T_Mat2))) +
                           kronecker(diag(nrow = nrow(T_Mat1)), T_Mat2),
                         kronecker(diag(nrow = nrow(T_Mat1)),rowSums(-T_Mat2)),
                         kronecker(rowSums(-T_Mat1),diag(nrow = nrow(T_Mat2))))

  newT_Mat <- rbind(firstblockrow,
                    cbind(matrix(0, ncol = p*q, nrow = p), T_Mat1, matrix(0, ncol = q, nrow = p)),
                    cbind(matrix(0, ncol = p*q, nrow = q), matrix(0, ncol = p, nrow = q), T_Mat2))

  return(contphasetype(initDist = newInitDist, T_Mat = newT_Mat))
}

