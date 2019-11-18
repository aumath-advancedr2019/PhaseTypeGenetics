#' The number of segregating sites
#'
#' The density function of the total number of segregating sites
#'
#' The density of the total number of segregating sites can be obtained
#' by the aid of the block counting process together with the reward
#' transformation and the discretization. For more information on this topic see \code{vignette("PhaseTypeGenetics")} or
#' Hobolth et al. (2019): \emph{Phase-type distributions in population genetics}.
#'
#' @param n the sample size (n >= 1).
#' @param theta the mutation parameter (theta > 0).
#' @param k a non-negative number or a non-negative vector.
#' @param plot a logical value indicating whether the function should
#' plot the density of the total number of segregating sites for the
#' given values of k.
#'
#' @return The function returns the probabilities \eqn{P(S=k)} for all values
#' of \eqn{k}. Hence, the returned object is of the same length as \eqn{k}.
#' If \code{plot=TRUE}, the function also plots the densities as a function of
#' \eqn{k}.
#'
#' @source Mogens Bladt and Bo Friis Nielsen (2017):
#' \emph{ Matrix-Exponential Distributions in Applied Probability}.
#' Probability Theory and Stochastic Modelling (Springer), Volume 81.
#'
#' @source Asger Hobolth, Arno Siri-Jégousse, Mogens Bladt (2019):
#' \emph{Phase-type distributions in population genetics}.
#' Theoretical Population Biology, 127, pp. 16-32.
#'
#' @seealso \code{\link{SiteFrequencies}}, \code{\link{dphasetype}}.
#'
#' @importFrom expm %^%
#'
#' @examples
#'
#' ## Computing the density for a sample of n=5
#' dSegregatingSites(n=5,theta=2,k=5)
#' dSegregatingSites(n=5,theta=2,k=1:20, plot=TRUE)
#'
#' ## We apply the function for different sample sizes
#' ## and theta=2
#' k_vec <- 0:15
#' theta <- 2
#' ## Defining a matrix of results
#' Res_Mat <- dSegregatingSites(n = 1, theta = theta, k = k_vec)
#' ## And Applying the function for all n in {2,...,20}
#' for(n in 2:20){
#'
#' Res_Mat <- cbind(Res_Mat, dSegregatingSites(n = n, theta = theta, k = k_vec))
#' }
#'
#' ## We reproduce Figure 4.1 in John Wakeley (2009):
#' ## "Coalescent Theory: An Introduction",
#' ## Roberts and Company Publishers, Colorado.
#' ## by using the package plot3D.
#' plot3D::hist3D(x=k_vec, y=1:20, z=Res_Mat, col = "grey", border = "black",
#'        xlab = "k", ylab = "n", zlab = "P(S=k)",
#'        main = "The probability function of the number of segregating sites",
#'        sub = expression(paste("The mutation parameter is ", theta,"= 2")),
#'        cex.main = 0.9, colkey = FALSE, zlim = c(0,0.4))
#'
#' @export
dSegregatingSites <- function(n, theta, k, plot =FALSE){

  if(n < 1) stop("Invalid sample size. n has to be greater than or equal to 1.")
  if(n != floor(n)) warning(paste("The proviede sample size n is not a natural number.\n
                   The function will use n= ", floor(n), " instead."))
  n = floor(n)

  if(theta <=0 ) stop("Invalid mutation parameter. Theta must be greater than 0.")
  if(sum(k<0)>0) stop("Invalid vector of quantiles. k has to be nonnegative!")
  if(!is.logical(plot)) stop(" 'plot' must be a logical value")

  if(n==1){

    res <- replicate(length(k),0)

  }else if(n==2){

    ## We define the reward transformed subtransition probability matrix
    P.mat <- theta/(theta+1)
    p.vec <- 1/(theta+1)
    ## We store the results in a matrix
    res <- (P.mat^k)*p.vec

  }else{
    ## For a given number n of samples, we find the state
    ## space and the corresponding rate matrix for the block
    ## counting process in the standard coalescent
    res <- BlockCountProcess(n)

    ## The rate matrix
    Tmat <- res$Rate.mat

    ## and the corresponding inital distribution
    pi.vec <- c(1,replicate(nrow(Tmat)-1,0))

    ## In order to find the distribution for the number
    ## of segregating sites, we need a reward vector that
    ## correpsonds to xi_1+...+_xi_n-1. Hence
    r.vec <- rowSums(res$StateSpace.mat)

    ## As all enties in the reward vector are positive, we
    ## can define the reward-transformed sub-intensity matrix
    ## by multiplying with the matrix that has 1/r(i) on its
    ## diagonal
    Tmat <- diag(1/r.vec)%*%Tmat

    ## Now we can compute the distribution of the number of
    ## segregating sites by using the descretization:
    P.mat <- solve(diag(nrow = nrow(Tmat))-(2/theta)*Tmat)

    res <- NULL

    for (i in k) {

      if(i%%1!=0){

        warning("One or more quantiles are not natural numbers.\n
              The corresponding probabilities are set to 0.")
        res[which(k==i)] <- 0
      }else{

        res[which(k==i)] <- pi.vec%*%(P.mat%^%i)%*%(1-rowSums(P.mat))
      }
    }
  }

  if(plot){

    plot(x=k, y=res, type = "l", col = "darkgrey",
         xlab = "k", ylab = expression(paste("P(", S["Total"], "=k)")),
         main = "The density function of the number of segregating sites",
         sub = paste("The mutation parameter is equal to", theta),
         cex.main = 0.9, ylim = c(0, max(res)*1.2))
  }
  return(res)
}
