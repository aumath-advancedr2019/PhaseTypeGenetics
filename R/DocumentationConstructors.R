#' Phase-type distributed objects
#'
#' Constructing phase-type representations
#'
#' \code{discphasetype} creates an object of class \emph{"discphasetype"}, i.e. the
#' object represents a discrete phase-type distribution.
#' \code{contphasetype}creates an object of class \emph{"contphasetype"}, i.e. the
#' object represents a continuous phase-type distribution.
#'
#' @param initDist a vector holding the initial distribution of the phase-type
#' distribution. Hence \code{initDist} should satisfy that all entries are non-negative
#' and the sum of the vector must be less than or equal to 1.
#' @param T_Mat A sub-intensity rate matrix. The matrix must be invertible, all diagonal
#' entries must be negative and all non-diagonal entries must be non-negative. Also for each
#' row the sum of the entries must be non-positive.
#' Furthermore, the length of the initial distribution has to be
#' equal to the number of rows and number of columns of \code{T_Mat}.
#' @param P_Mat A sub-transition probability matrix satisfying that all entries
#' are non-negative. Furthermore, \code{P_Mat} should be a square matrix, which
#' determinant is different from zero. All rows in \code{P_Mat} should sum to a
#' number less than or equal to one, and the length of the initial distribution has to be
#' equal to the number of rows of \code{P_Mat}.
#'
#' @return \code{discphasetype} returns a list of class \emph{"discphasetype"}, while
#' \code{contphasetype} returns a list of class \emph{"contphasetype"}. Both lists
#' hold
#' \itemize{
#'   \item \code{initDist} : The initial distribution
#'   \item \code{P_Mat}/\code{T_Mat} : The sub-transition/sub-intensity matrix
#' }
#'
#' @examples
#' ## For n=4, the time to the most recent common ancestor is
#' ## phase-type distributed with
#' ## initial distribution
#' initDist <- c(1,0,0)
#' ## and sub-intensity rate matrix
#' T_Mat <- matrix(c(-6,6,0,
#'                    0,-3,3,
#'                    0,0,-1), nrow = 3, ncol = 3, byrow = TRUE)
#'
#' TMRCA <- contphasetype(initDist, T_Mat)
#'
#' ## For theta=2, the number of segregating sites plus one is
#' ## discrete phase-type distributed with
#' ## initial distribution
#' initDist <- c(1,0,0,0)
#' ## and sub-transition probability matrix
#' P_Mat <- matrix(c(0.4, 0.3, 4/30, 2/30,
#'                    0, 0.5, 2/9, 1/9,
#'                    0, 0, 2/3, 0,
#'                    0, 0, 0, 2/3), nrow = 4, ncol = 4, byrow = TRUE)
#'
#' S_Total <- discphasetype(initDist, P_Mat)
#'
#' @describeIn discphasetype Creating an object of type discphasetype
#' @export
discphasetype <- function(initDist, P_Mat){
  P_Mat <- as.matrix(P_Mat)
  if(sum(initDist) > 1) stop("Not a valid initial distribution, as sum(initDist) > 1")
  if(sum(initDist < 0) > 0) stop("Not a valid initial distribution, as some entries are negative")
  if(nrow(P_Mat)!= ncol(P_Mat)) stop("The subtransition matrix is not a square matrix")
  if(det(P_Mat) == 0) stop("A singular matrix is not a valid subtransition matrix")
  if(sum(P_Mat < 0) > 0) stop("Not a valid subtransition matrix, as some entries are negative")
  if(sum(rowSums(P_Mat) > 1) > 0) stop("Not a valid subtransition matrix. All rows should sum to a number less than or equal to one")
  if(length(initDist) != nrow(P_Mat)) stop("The dimensions of the input should be the same")
  tmp <- list("initDist" = initDist, "P_Mat" = P_Mat)
  class(tmp) <- "discphasetype"
  return(tmp)
}

#' @describeIn discphasetype Creating an object of type contphasetype
#' @export
contphasetype <- function(initDist, T_Mat){
  T_Mat <- as.matrix(T_Mat)
  if(sum(initDist) > 1) stop("Not a valid initial distribution, as sum(initDist) > 1")
  if(sum(initDist < 0) > 0) stop("Not a valid initial distribution, as some entries are negative")
  if(nrow(T_Mat)!= ncol(T_Mat)) stop("The subintensity matrix is not a square matrix")
  if(sum(diag(T_Mat)>=0) > 0) stop("All diagonal entries of the subintensity matrix must be negative")
  if(det(T_Mat)==0) stop("A singular matrix is not a valid subintensity matrix")
  if(sum(round(rowSums(T_Mat),9) > 0) > 0) stop("The sums of the rows must be non-positive")
  if(length(initDist) != nrow(T_Mat)) stop("The dimensions of the input should be the same")
  tmp <- list("initDist" = initDist, "T_Mat" = T_Mat)
  class(tmp) <- "contphasetype"
  return(tmp)
}
