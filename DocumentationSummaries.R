#' Summaries of phase-type distributed objects
#'
#' These \code{summary} functions are \code{\link{methods}} which depend on the \code{\link{class}} of the first argument.
#' For objects of class \code{discphasetype} or \code{contphasetype},
#' \code{summary} prints the inital distribution and the
#' subtransition/subintensity matrix of the phase-type distribution.
#'
#'
#' @param object an object for which a summary is desired. To be able to use
#' the summaries for phase-type distributed objects,the object has to be of
#' class \code{discphasetype} or \code{contphasetype}.
#'
#' @return \code{summary.discphasetype} and \code{summary.contphasetype} print
#' the initial distribution and the subtransition/subintensity matrix of the
#' phase-type distribution. If the sum of all entries in the initial distribution
#' is less than one, it also prints the defect size \code{(1-sum(initDist))}.
#'
#' @seealso \code{\link{summary}}, \code{\link{summary.glm}},
#' \code{\link{summary.lm}}.
#'
#' @examples
#'
#'
#' @export

summary.discphasetype <- function(object){
  cat("A discrete phase-type distribution with initial probability vector \n")
  print(object$initDist)
  cat("and subtransition matrix \n")
  print(object$P.mat)
  if(sum(object$initDist) < 1){
    cat("and defect\n", 1-sum(object$initDist))
  }
}

#' @rdname add
summary.contphasetype <- function(object){
  cat("A continuous phase-type distribution with initial probability vector \n")
  print(object$initDist)
  cat("and subintensity matrix \n")
  print(object$T.mat)
  if(sum(object$initDist) < 1){
    cat("and defect\n", 1-sum(object$initDist))
  }
}
