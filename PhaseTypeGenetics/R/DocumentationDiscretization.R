#' Discretizing a continuous phase-type disctribution
#'
#' Discretizes a continuous phase-type distribution with initial distribution
#' \code{initDist} and subintensity rate matrix \code{T.mat}.
#'
#' The relation between continuous and discrete phase-type distributions
#' is given in the following way. If \eqn{T} is the subintensity rate matrix of a continuous
#' phase-type distribution with representation \eqn{PH(pi,T)}, then there exists
#' a constant \eqn{a>0} such that \eqn{P := I + 1/a * T} is a subtransition
#' probability matrix and \eqn{DPH(pi, P)} is a representation for a discrete
#' phase-type distribution. This holds for any \eqn{a} larger than the maximum of
#' all diagonal entries in \eqn{T}, as all entries in a subtransition
#' probability matrix have to be between zero and one.
#' It even holds that for a genealogical model where the total
#' branch length \eqn{\tau ~ PH(pi, T)} and the mutation rate at the locus is \eqn{\lambda = \theta/2},
#' the number of segregating sites \eqn{S} plus one is discrete phase-type distributed
#' with inital distribution \eqn{pi} and subtransition probability matrix
#' \eqn{P = (I-\lambda^{-1} * T)^{-1}}, i.e.
#' \deqn{S + 1 ~ DPH(pi, P).}
#'
#' @param object a continuous phase-type distributed object of class \code{contphasetype}.
#' @param a a constant that is larger than the maximum of all diagonal
#' entries of the subintensity rate matrix
#' @param lambda the positive mutation rate at the locus
#'
#' @return Depending on the input, the function returns the discretized phase-type
#' disctribution with subtransition probability matrix equal to either
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
#' ## and subintensity rate matrix
#' T.mat <- matrix(c(-1.5, 1.5, 0, 0,
#'                   0, -1, 2/3, 1/3,
#'                   0, 0, -0.5, 0,
#'                   0, 0, 0, -0.5), nrow = 4, byrow = TRUE)
#'
#' T_Total <- contphasetype(initDist,T.mat)
#'
#' ## Hence, for theta=2, the number of segregating sites plus one is
#' ## discrete phase-type distributed with the same initial
#' ## distribution and subtransition probability matrix
#' discretization(T_Total, lambda=1)$P.mat
#'
#' @export
discretization <- function(object, a=NULL, lambda=NULL){

  if(class(object) != "contphasetype") stop("Invalid object! The object has to be of class 'contphasetype'.")

  if(is.null(lambda)){

    if(is.null(a) | a < max(diag(object$T.mat)) ){stop("Not a valid constant a! a has to be bigger than the maximum of all diagonal entries of the subintensity rate matrix.")}

    res <- discphasetype(initDist = object$initDist,
                  P.mat = diag(nrow = nrow(object$T.mat))+1/a*object$T.mat)

  }else if(is.null(a)){

    if(is.null(lambda) | lambda <= 0){stop("Not a valid mutation rate. lambda has to be positive.")}

    res <- discphasetype(initDist = object$initDist,
                         P.mat = solve(diag(nrow = nrow(object$T.mat))-lambda^(-1)*object$T.mat))


  }else{

    if(lambda <= 0){stop("Not a valid mutation rate. lambda has to be positive.")}
    if(a < max(diag(object$T.mat)) )stop("Not a valid constant a! a has to be bigger than the maximum of all diagonal entries of the subintensity rate matrix.")

    res <- list(a = discphasetype(initDist = object$initDist,
                                  P.mat = diag(nrow = nrow(object$T.mat))+1/a*object$T.mat),
                lambda = discphasetype(initDist = object$initDist,
                                       P.mat = solve(diag(nrow = nrow(object$T.mat))-lambda^(-1)*object$T.mat)))
  }
  return(res)
}

