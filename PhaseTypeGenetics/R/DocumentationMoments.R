#' Statistical moments of phase-type distributions
#'
#' Computing (factorial) moments of a given order of phase-type distributed random variables.
#'
#' In the discrete case \eqn{( \tau ~ DPH(initDist,P) )}, the factorial moments are given by
#' \deqn{E[\tau(\tau-1) ··· (\tau-i+1)] = i! initDist P^(i-1) (I-P)^(-i) e,}
#' where \code{initDist} is the initial distribution and \eqn{P} is the sub-transition probability matrix.
#' For \eqn{\tau ~ PH(initDist, T)}, the \eqn{i}'th-order moment is defined as
#' \deqn{E[\tau^i] = i! initDist (-T)^(-i) e,}
#' where \code{initDist} is again the inital distribution and \eqn{T} is the subintensity rate matrix.
#' In both cases, \eqn{e} is a vector with one in each entry.
#'
#' @param object either a continuous phase-type distributed object of class \code{contphasetype} or
#' a discrete phase-type distributed object of class \code{discphasetype}.
#' @param i a positive number stating the order of the desired moment
#' @param all a logical value indicating whether the function should compute
#' all moments up to the given order. The default is equal to FALSE.
#'
#' @return For \code{all = FALSE}, the function either returns the \eqn{i}'th-order moment (if the object is continuous
#' phase-type distributed) or the \eqn{i}'th factorial moment (if the object is discrete
#' phase-type distributed). In both cases, the length of the output is one.
#' For \code{all = TRUE}, the function computes all (factorial) moments up to the given order,
#' hence the output is a vector of length \eqn{i}.
#'
#' @source Mogens Bladt and Bo Friis Nielsen (2017):
#' \emph{ Matrix-Exponential Distributions in Applied Probability}.
#' Probability Theory and Stochastic Modelling (Springer), Volume 81.
#'
#'
#' @examples
#'
#' ## Using the function moments() to compute the mean
#' ## and variance of a phase-type distribution
#'
#' ## For n=4, the time to the most recent common ancestor is
#' ## phase-type distributed with initital distribution
#' initDist <- c(1,0,0)
#' ## and sub-intensity rate matrix
#' Tmat <- matrix(c(-6,6,0,
#'                    0,-3,3,
#'                    0,0,-1), nrow = 3, ncol = 3, byrow = TRUE)
#' ## Defining an object of type "contphasetype"
#' TMRCA <- contphasetype(initDist, Tmat)
#' ## Computing all moments up to order 2
#' m <- moments(TMRCA, i=2, all = TRUE)
#' ## We get the desired numbers
#' m[1]
#' phmean(TMRCA)
#'
#' m[2] - m[1]^2
#' phvar(TMRCA)
#'
#' ## For theta=2, the number of segregating sites plus one is
#' ## discrete phase-type distributed with initital distribution
#' initDist <- c(1,0,0,0)
#' ## and sub-transition probability matrix
#' Pmat <- matrix(c(0.4, 0.3, 4/30, 2/30,
#'                    0, 0.5, 2/9, 1/9,
#'                    0, 0, 2/3, 0,
#'                    0, 0, 0, 2/3), nrow = 4, ncol = 4, byrow = TRUE)
#' ## Defining an object of type "discphasetype"
#' S_Total <- discphasetype(initDist, Pmat)
#' ## Computing all moments up to order 2
#' m <- moments(S_Total, i=2, all = TRUE)
#' ## We get the desired numbers
#' m[1]
#' phmean(S_Total)
#'
#' m[2] + m[1] - m[1]^2
#' phvar(S_Total)
#'
#' @export
moments <- function(object, i, all = FALSE){

  UseMethod("moments")
}

#' @export
moments.default <- function(object, i, all = FALSE){

  print("The package ‘moments’ can be used to compute sample moments of specified order.")
}

#' @export
#' @importFrom expm %^%
moments.discphasetype <- function(object, i, all = FALSE){

  if(i < 1 ) stop("Invalid order. The number i has to be positive!")
  if(!is.logical(all)) stop(" 'all' must be a logical value")

  initDist = object$initDist
  P.mat = object$P.mat

  res <- NULL
  if(all == FALSE){

    res <- factorial(i)*sum(initDist%*%(P.mat%^%(i-1)%*%
           (solve(diag(1, nrow = nrow(P.mat))-P.mat)%^%(i))))
  }else{
    for (k in 1:i){

      res[k] <- factorial(k)*sum(initDist%*%(P.mat%^%(k-1)%*%
                (solve(diag(1, nrow = nrow(P.mat))-P.mat)%^%(k))))
    }
    names(res) <- 1:i
  }
  return(res)
}

#' @export
#' @importFrom expm %^%
moments.contphasetype <- function(object, i, all = FALSE){

  if(i < 1 ) stop("Invalid order. The number i has to be positive!")
  if(!is.logical(all)) stop(" 'all' must be a logical value")

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
