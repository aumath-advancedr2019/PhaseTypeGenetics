#' Phase-type distributed objects
#'
#' Constructing phase-type distributions
#'
#' \code{discphasetype} creates an object of class "discphasetype", i.e. the
#' object represents a discrete phase-type distribution.
#' \code{contphasetype}creates an object of class "contphasetype", i.e. the
#' object represents a continuous phase-type distribution.
#'
#' @param initDist a vector holding the initial distribution of the phase-type
#' distribution. Hence \code{initDist} should satisfy that all entries sum to one,
#' and that they are all non-negative.
#' @param T.mat A subintensity rate matrix satisfying that \code{det(T.mat)} is not
#' equal to zero. Furthermore, the length of the initial distribution has to be
#' equal to the number of rows and number of columns of \code{P.mat}.
#' @param P.mat A subtransition probability matrix satisfying that all entries
#' are non-negative. Furthermore, \code{P.mat} should be a square matrix, which
#' determinant is different from zero. All rows in \code{P.mat} should sum to a
#' number less than or equal to one, and the length of the initial distribution has to be
#' equal to the number of rows of \code{P.mat}.
#'
#' @return \code{discphasetype} returns a list of class "discphasetype", while
#' \code{contcphasetype} returns a list of class "contphasetype". Both lists
#' hold
#' \describe{
#'   \item{\code{initDist}}{The inital distribution }
#'   \item{\code{P.mat}/\code{T.mat}}{The subtransition/subintensity matrix}
#' }
#'
#' @examples
#'
#'
#' @export

discphasetype <- function(initDist, P.mat){
  P.mat <- as.matrix(P.mat)
  if(sum(initDist) > 1) stop("Not a valid initial distribution, as sum(initDist) > 1")
  if(sum(initDist < 0) > 0) stop("Not a valid initial distribution, as some entries are negative")
  if(nrow(P.mat)!= ncol(P.mat)) stop("The subtransition matrix is not a square matrix")
  if(det(P.mat) == 0) stop("Not a valid subtransition matrix, as det(P.mat) = 0")
  if(sum(P.mat < 0) != 0) stop("Not a valid subtransition matrix, as some entries are negative")
  if(sum(rowSums(P.mat) > 1) > 0) stop("Not a valid subtransition matrix. All rows should sum to a number less than or equal to one")
  if(length(initDist) != nrow(P.mat)) stop("The dimensions of the input should be the same")
  tmp <- list("initDist" = initDist, "P.mat" = P.mat)
  class(tmp) <- "discphasetype"
  return(tmp)
}

#' @rdname add
contphasetype <- function(initDist, T.mat){
  T.mat <- as.matrix(T.mat)
  if(sum(initDist) > 1) stop("Not a valid initial distribution, as sum(initDist) > 1")
  if(sum(initDist < 0) > 0) stop("Not a valid initial distribution, as some entries are negative")
  if(nrow(T.mat)!= ncol(T.mat)) stop("The subintensity matrix is not a square matrix")
  if(det(T.mat)==0) stop("Not a valid subintensity matrix, as det(T.mat)=0")
  if(length(initDist) != nrow(T.mat)) stop("The dimensions of the input should be the same")
  tmp <- list("initDist" = initDist, "T.mat" = T.mat)
  class(tmp) <- "contphasetype"
  return(tmp)
}
