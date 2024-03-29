#' Sum of phase-type distributions
#'
#' The sum of two independent discrete or continuous
#' phase-type distributed variables with initial distributions
#' \code{initDist1} and \code{initDist2} as well as sub-transition/sub-intensity
#' matrices equal to \code{P_Mat1}/\code{T_Mat1} and  \code{P_Mat2}/\code{T_Mat2}.
#'
#' In the discrete case, the sum of two phase-type distributed
#' variables \eqn{tau1 ~ DPH_p(\alpha,S)} and \eqn{tau2 ~ DPH_q(\beta,T)} is
#' again discrete phase-type distributed in the following way
#' \deqn{tau1 + tau2 ~ DPH_{p+q}((\alpha,0),cbind((S, s \beta),(0,T)) ).}
#' In the continuous case, the sum of two phase-type distributed
#' variables \eqn{X ~ PH_p(\alpha,S)} and \eqn{Y ~ PH_q(\beta,T)} is
#' again continuous and phase-type distributed in the following way
#' \deqn{X + Y ~ PH_{p+q}((\alpha,0),cbind((S, s \beta),(0,T)) ).}
#'
#' @param object1,object2 two objects of class \code{discphasetype}
#' or \code{contphasetype} for which the sum should be computed.
#'
#' @return The function \code{phsum} returns an object of type \code{discphasetype}
#' or \code{contphasetype} (depending on the input) holding the phase-type representation
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
#' ## A simple example
#' phsum(T_MRCA$n5,T_Total$n5)
#'
#' ## For n=4, the total length of branches giving rise to
#' ## singletons is phase-type distributed with initial distribution
#' initDist1 <- c(1,0,0)
#' ## and sub-intensity rate matrix
#' T_Mat1 <- matrix(c(-1.5, 1.5, 0,
#'                    0, -1.5, 1,
#'                    0, 0, -1), nrow = 3, byrow = TRUE)
#' ## The total length of branches giving rise to
#' ## double-tons is phase-type distributed with initial distribution
#' initDist2 <- c(1,0)
#' ## and sub-intensity rate matrix
#' T_Mat2 <- matrix(c(-3, 1,
#'                    0, -0.5), nrow = 2, byrow = TRUE)
#' ## Defining two objects of type "contphasetype"
#' T1 <- contphasetype(initDist1, T_Mat1)
#' T2 <- contphasetype(initDist2, T_Mat2)
#' ## Hence, the total length of branches giving rise to
#' ## singletons and doubletons is phase-type distributed
#' ## in the following way
#' phsum(T1,T2)
#' ## (Please compare this distribution with the distribution
#' ## obtained directly from the reward transformation)
#'
#' @export
phsum <- function(object1,object2){

  UseMethod("phsum")
}

#' @export
phsum.default <- function(object1,object2){

  sum(object1,object2)
}

#' @export
phsum.discphasetype <- function(object1,object2){

  if(class(object1) != "discphasetype"| class(object2) != "discphasetype") stop("Invalid objects! object1 and object2 must be of class 'discphasetype'.")

  initDist1 = object1$initDist
  P_Mat1 = object1$P_Mat
  initDist2 = object2$initDist
  P_Mat2 = object2$P_Mat

  newInitDist <- c(initDist1,rep(0,length(initDist2)))
  newP_Mat <- rbind(
    cbind(P_Mat1,(1-rowSums(P_Mat1))%*%t(initDist2)),
    cbind(matrix(0, nrow = length(initDist2), ncol = length(initDist1)), P_Mat2))

  return(discphasetype(initDist = newInitDist, P_Mat = newP_Mat))

}

#' @export
phsum.contphasetype <- function(object1,object2){

  if(class(object1) != "contphasetype"| class(object2) != "contphasetype") stop("Invalid objects! object1 and object2 must be of class 'contphasetype'.")
  initDist1 = object1$initDist
  T_Mat1 = object1$T_Mat
  initDist2 = object2$initDist
  T_Mat2 = object2$T_Mat

  newInitDist <- c(initDist1,rep(0,length(initDist2)))
  newT_Mat <- rbind(
    cbind(T_Mat1, rowSums(-T_Mat1) %*% t(initDist2)),
    cbind(matrix(0, nrow = length(initDist2), ncol = length(initDist1)), T_Mat2))

  return(contphasetype(initDist = newInitDist, T_Mat = newT_Mat))
}
