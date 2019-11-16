#' Summaries of phase-type distributed objects
#'
#' This \code{summary} function is a method that depends on the \code{\link{class}} of the argument.
#' For objects of class \code{discphasetype} or \code{contphasetype},
#' \code{summary} prints the inital distribution and the
#' subtransition/subintensity matrix of the phase-type distribution.
#'
#' @param object an object for which a summary is desired. To be able to use
#' the summaries for phase-type distributed objects,the object has to be of
#' class \code{discphasetype} or \code{contphasetype}.
#'
#' @return \code{summary} prints
#' the initial distribution and the subtransition/subintensity matrix of the
#' phase-type distribution. If the sum of all entries in the initial distribution
#' is less than one, it also prints the defect size \code{(1-sum(initDist))}.
#'
#' @seealso \code{\link{summary}}, \code{\link{summary.glm}},
#' \code{\link{summary.lm}}.
#'
#' @examples
#'
#' summary(T_Total$n10)
#'
#' @export
summary.discphasetype <- function(object){

  cat("A discrete phase-type distribution with initial probability vector \n")

  if(length(object$initDist) < 10){

    print(object$initDist)
    cat("and subtransition matrix \n")
    print(object$P.mat)
  }else{

    print(object$initDist[1:10])
    cat("...\n")
    cat("and subtransition matrix \n")
    print(object$P.mat[1:10,1:10])
    cat("... \n")
    cat("(Showing only the first ten entries)\n")
  }

  if(sum(object$initDist) < 1){
    cat("and defect\n", 1-sum(object$initDist))
  }
}

#' @rdname summary.discphasetype
#' @export
summary.contphasetype <- function(object){

  cat("A continuous phase-type distribution with initial probability vector \n")

  if(length(object$initDist) < 10){

    print(object$initDist)
    cat("and subintensity matrix \n")
    print(object$T.mat)
  }else{

    print(object$initDist[1:10])
    cat("... \n")
    cat("and subintensity matrix \n")
    print(object$T.mat[1:10,1:10])
    cat("... \n")
    cat("(Showing only the first ten entries)\n")
  }
  if(sum(object$initDist) < 1){
    cat("and defect\n", 1-sum(object$initDist))
  }
}
