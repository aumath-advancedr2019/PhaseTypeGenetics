#' Discretizing a continuous phase-type disctribution
#'
#' Discretizes a continuous phase-type distribution with initial distribution
#' \code{initDist} and subintensity rate matrix \code{T.mat}.
#'
#' The relation between continuous and discrete phase-type distributions 
#' is given in the following way. If T is the subintensity rate matrix of a continuous
#' phase-type distribution with representation PH(pi,T), then there exists
#' a constant a>0 such that P := I + 1/a * T is a subtransition 
#' probability matrix and DPH(pi, P) is a representation for a discrete 
#' phase-type distribution. This holds for any a larger than the maximum of 
#' all diagonal entries in T, as all entries in a subtransition 
#' probability matrix have to be between zero and one. 
#' This relation even implies that for a genealogical model where the total
#' branch length tau ~ PH(pi, T) and the mutation rate at the locus is lambda = theta/2,
#' the number of segregating sites S plus one is discrete phase-type distributed
#' with inital disctribution pi and subtransition probability matrix 
#' P = (I-lambda^{-1}*T)^{-1}, i.e. 
#' S + 1 ~ DPH(pi, P).
#' 
#' @param object a continuous phase-type distributed object of class \code{contphasetype}.
#' @param a a constant that is larger than the maximum of all diagonal 
#' entries of the subintensity rate matrix
#' @param lambda the mutation rate at the locus
#'
#' @return Depending on the input, the function returns the discretized phase-type
#' disctribution with subtransition probability matrix equal to either 
#' P := I + 1/a * T (if a is provided) or P = (I-lambda^{-1}*T)^{-1} 
#' (if lambda is provided). If both a and lambda are provided, the function
#' returns both distributions in a list. In all three cases, the returned objects are
#' of type \code{discphasetype}.
#'
#' @source Mogens Bladt and Bo Friis Nielsen (2017):
#' \emph{ Matrix-Exponential Distributions in Applied Probability}.
#' Probability Theory and Stochastic Modelling (Springer), Volume 81.
#' 
#'
#' @examples
#'
#'
#' @export

discretization <- function(object, a=NULL, lambda=NULL){
  
  if(class(object) != "contphasetype") stop("The object has to be of type contphasetype")
  
  if(is.null(lambda)){
    
    if(is.null(a) | a < max(diag(object$T.mat)) ){stop("Not a valid constant a")}
    
    res <- discphasetype(initDist = object$initDist,
                  P.mat = diag(nrow = nrow(object$T.mat))+1/a*object$T.mat)
    
  }else if(is.null(a)){
    
    if(is.null(lambda)){stop("Not a valid mutation rate")}
    
    res <- discphasetype(initDist = object$initDist,
                         P.mat = solve(diag(nrow = nrow(object$T.mat))-lambda^(-1)*object$T.mat))
    
    
  }else{
    
    res <- list(a = discphasetype(initDist = object$initDist,
                                  P.mat = diag(nrow = nrow(object$T.mat))+1/a*object$T.mat),
                lambda = discphasetype(initDist = object$initDist,
                                       P.mat = solve(diag(nrow = nrow(object$T.mat))-lambda^(-1)*object$T.mat)))
  }
  return(res)
}

