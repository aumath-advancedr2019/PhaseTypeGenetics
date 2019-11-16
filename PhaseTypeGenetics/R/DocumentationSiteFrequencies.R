#' Site frequency spectrum
#'
#' Finding the distribution of the site frequencies, the number of
#' segregating sites and the tail statistics using phase-type distributions.
#'
#' This function can be used to compute the distribution of the site frequencies
#' \eqn{\xi_i + 1}, for all \eqn{i} in eqn{{1,...,n-1}}, the total number of segregating sites \eqn{S_{Total} + 1}
#' and the tail statistic \eqn{S_{i+} + 1}.
#' The reason for adding one to the site frequency is that the support for
#' discrete phase-type distributions is on the natural numbers excluding zero.
#' Hence, immediate absorbtion would not be possible. By adding one, we allow
#' the site frequency to be zero.
#' Note that the package does also include the function \code{\link{dSegregatingSites}},
#' which computes the density function of the number of segregating sites
#' for a given sample size \eqn{n}, a mutation parameter \eqn{\theta} and
#' a nonnegative vector of quantiles \eqn{k}.
#'
#' @param n the sample size (>=3)
#' @param lambda the nonnegative mutation rate
#' @param i either the number of the site frequency that should be considered
#' or the number of the first term of the tail statistic. In both cases \eqn{1 <= i <= n-1}.
#' @param nSegSites a logical value indicating whether the function should compute
#' the distribution of the number of segregating sites \eqn{(S_{Total} = \xi_1 + ... + \xi_{n-1})}.
#' If TRUE, any value of \eqn{i} will be ignored. Defaults to FALSE.
#' @param tailStat a logical value indicating whether the function should compute
#' the distribution of the tail statistic \eqn{( S_{i+} = \xi_i +...+ \xi_{n-1})}.
#' If TRUE, \eqn{i} will determine the first term of this statistic. Defaults to FALSE.
#'
#' @return If \code{nSegSites = FALSE} and \code{tailStat= FALSE}, the function
#' returns the distribution of the \eqn{i}'th site frequency \eqn{(\xi_i)} plus one. If \code{nSegSites = TRUE},
#' the function returns the distribution of the total number of segregating sites plus one, and if
#' \code{tailStat= TRUE}, the distribution of the tail statistic (which first term
#' is determined by \eqn{i}) plus one is returned.
#' In all three cases, the returned object is of type \code{discphasetype}.
#'
#' @source Asger Hobolth, Arno Siri-Jégousse, Mogens Bladt (2019):
#' \emph{Phase-type distributions in population genetics}.
#' Theoretical Population Biology, 127, pp. 16-32.
#'
#' @examples
#' SiteFrequencies(n=4, lambda=1, i=2)
#' SiteFrequencies(n=4, lambda=1, nSegSites=TRUE)
#' SiteFrequencies(n=4, lambda=1, i=2, tailStat=TRUE)
#'
#' @seealso \code{\link{dSegregatingSites}}
#'
#' @export
SiteFrequencies <- function(n,lambda,i=NULL, nSegSites=FALSE, tailStat = FALSE){

  if(n < 2) stop("Not a valid sample size. n must be greater than 2.")
  if(lambda < 0) stop("Not a valid mutation rate. lambda must be nonnegative!")
  if(i <1 | i>n-1) stop("i must be between 1 and n-1")
  if(!is.logical(nSegSites)) stop("nSegSites must be a logical value")
  if(!is.logical(tailstat)) stop("tailStat must be a logical value")

  ## For a given number n of samples, we find the state
  ## space and the corresponding rate matrix for the block
  ## counting process in the standard coalescent
  res <- BlockCountProcess(n)

  ## The rate matrix
  Tmat <- res$Rate.mat
  ## and the corresponding inital distribution
  pi.vec <- c(1,replicate(nrow(Tmat)-1,0))

  ## We define an object of type contphasetype
  obj <- contphasetype(pi.vec, Tmat)

  if(nSegSites){

    ## In order to find the distribution for the number
    ## of segregating sites, we need a reward vector that
    ## correpsonds to xi_1+...+_xi_n-1. Hence
    r.vec <- rowSums(res$StateSpace.mat)

    ## As all enties in the reward vector are positive, we
    ## can define the reward-transformed sub-intensity matrix
    ## by multiplying with the matrix that has 1/r(i) on its
    ## diagonal
    Tmat <- diag(1/r.vec)%*%Tmat

    obj <- contphasetype(pi.vec, Tmat)

    ## Now we can compute the distribution of the number of
    ## segregating sites by using the descretization:
    newobj <- discretization(obj, a=NULL, lambda=lambda)

  }else if(tailStat){

    ## In order to find the distribution for the tail
    ## statistic, we need a reward vector that
    ## correpsonds to xi_i+...+_xi_n-1. Hence
    if(length(i:(n-1))==1){
      r.vec <- res$StateSpace.mat[,(n-1)]

    }else{
      r.vec <- rowSums(res$StateSpace.mat[,i:(n-1)])
    }

    ## In this case, some of the enties in the reward vector
    ## are zero. Therefore, we have to use the function
    ## RewTransDistribution in order to get the transformed
    ## distribution
    newobj <- RewTransDistribution(obj, r.vec)

    ## Now we can compute the distribution of the
    ## tail statistic by using the descretization:
    newobj <- discretization(newobj, a=NULL, lambda=lambda)

  }else if(nSegSites == FALSE & tailStat == FALSE & is.null(i)==FALSE){

    if(i <1 | i>n-1) stop("i must be between 1 and n-1")

    ## In order to find the distribution for the site
    ## frequency xi_i, we need a reward vector that
    ## correpsonds to xi_i. Hence
    r.vec <- res$StateSpace.mat[,i]

    ## In this case, some of the enties in the reward vector
    ## are zero. Therefore, we have to use the function
    ## RewTransDistribution in order to get the transformed
    ## distribution
    newobj <- RewTransDistribution(obj, r.vec)

    ## Now we can compute the distribution of the
    ## tail statistic by using the descretization:
    newobj <- discretization(newobj, a=NULL, lambda=lambda)

  }else{

    stop("i is not specified.")
  }
  return(newobj)
}
