#' Discretizing a continuous phase-type distribution
#'
#' Discretizes a continuous phase-type distribution with initial distribution
#' \code{initDist} and sub-intensity rate matrix \code{T_Mat}.
#'
#' The relation between continuous and discrete phase-type distributions
#' is given in the following way. If \eqn{T} is the sub-intensity rate matrix of a continuous
#' phase-type distribution with representation \eqn{PH(initDist,T)}, then there exists
#' a constant \eqn{a>0} such that \eqn{P := I + 1/a * T} is a sub-transition
#' probability matrix and \eqn{DPH(initDist, P)} is a representation for a discrete
#' phase-type distribution. This holds for any \eqn{a} larger than the maximum of
#' all diagonal entries in \eqn{T}, as all entries in a sub-transition
#' probability matrix have to be between zero and one.
#' It even holds that for a genealogical model where the total
#' branch length \eqn{\tau ~ PH(initDist, T)} and the mutation rate at the locus is \eqn{\lambda = \theta/2},
#' that the number of segregating sites \eqn{S_{Total}} plus one is discrete phase-type distributed
#' with initial distribution \eqn{initDist} and sub-transition probability matrix
#' \eqn{P = (I-\lambda^{-1} * T)^{-1}}, i.e.
#' \deqn{S + 1 ~ DPH(initDist, P).}
#'
#' @param object a continuous phase-type distributed object of class \code{contphasetype}.
#' @param a a constant that is larger than the maximum of all diagonal
#' entries of the sub-intensity rate matrix.
#' @param lambda the positive mutation rate at the locus.
#'
#' @return Depending on the input, the function returns the discretized phase-type
#' distribution with sub-transition probability matrix equal to either
#' \deqn{P := I + 1/a * T }
#' (if \eqn{a} is provided) or
#' \deqn{P = (I-lambda^{-1} * T)^{-1}}
#' (if \eqn{\lambda} is provided). If both \eqn{a} and \eqn{\lambda} are provided, the function
#' returns both distributions in a list. In all three cases, the returned objects are
#' of type \code{discphasetype}.
#'
#' @source Mogens Bladt and Bo Friis Nielsen (2017):
#' \emph{ Matrix-Exponential Distributions in Applied Probability}.
#' Probability Theory and Stochastic Modelling (Springer), Volume 81.
#'
#'
#' @examples
#' ## For n=4, the total branch length is phase-type
#' ## distributed with initial distribution
#' initDist <- c(1,0,0,0)
#' ## and sub-intensity rate matrix
#' T_Mat <- matrix(c(-1.5, 1.5, 0, 0,
#'                   0, -1, 2/3, 1/3,
#'                   0, 0, -0.5, 0,
#'                   0, 0, 0, -0.5), nrow = 4, byrow = TRUE)
#'
#' TTotal <- contphasetype(initDist,T_Mat)
#'
#' ## Hence, for theta=2, the number of segregating sites plus one is
#' ## discrete phase-type distributed with the same initial
#' ## distribution and sub-transition probability matrix
#' discretization(TTotal, lambda=1)$P_Mat
#'
#' @export
discretization <- function(object, a=NULL, lambda=NULL){

  if(class(object) != "contphasetype") stop("Invalid object! The object has to be of class 'contphasetype'.")

  if(is.null(lambda)){

    if(is.null(a) | a < max(diag(object$T_Mat)) ){stop("Not a valid constant a! a has to be bigger than the maximum of all diagonal entries of the subintensity rate matrix.")}

    res <- discphasetype(initDist = object$initDist,
                  P_Mat = diag(nrow = nrow(object$T_Mat))+1/a*object$T_Mat)

  }else if(is.null(a)){

    if(is.null(lambda) | lambda <= 0){stop("Not a valid mutation rate. lambda has to be positive.")}

    res <- discphasetype(initDist = object$initDist,
                         P_Mat = solve(diag(nrow = nrow(object$T_Mat))-lambda^(-1)*object$T_Mat))


  }else{

    if(lambda <= 0){stop("Not a valid mutation rate. lambda has to be positive.")}
    if(a < max(diag(object$T_Mat)) )stop("Not a valid constant a! a has to be bigger than the maximum of all diagonal entries of the subintensity rate matrix.")

    res <- list(a = discphasetype(initDist = object$initDist,
                                  P_Mat = diag(nrow = nrow(object$T_Mat))+1/a*object$T_Mat),
                lambda = discphasetype(initDist = object$initDist,
                                       P_Mat = solve(diag(nrow = nrow(object$T_Mat))-lambda^(-1)*object$T_Mat)))
  }
  return(res)
}

