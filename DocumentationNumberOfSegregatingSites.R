#' The number of segregating sites
#'
#' The density function of the total number of segregating sites 
#'
#' The density of the total number of segregating sites can be obtained
#' by the aid of the block counting process together with the reward 
#' transformation and the discretization. For more information on this topic see
#' Hobolth et al. (2019): \emph{Phase-type distributions in population genetics}.
#' 
#' @param n the sample size
#' @param theta the mutation parameter
#' @param k a nonnegative number or a nonnegative vector 
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
#' @examples
#'
#' dSegregatingSites(n=5,theta=2,k=5)
#' dSegregatingSites(n=5,theta=2,k=1:20, plot=TRUE)
#' 
#' ## We apply the function for different sample sizes
#' ## and 
#' k.vec <- 0:15
#' theta <- 2
#' 
#' Res.mat <- dSegregatingSites(n = 1, theta = theta, k = k.vec)
#'
#' for(n in 2:20){
#' 
#' Res.mat <- cbind(Res.mat, dSegregatingSites(n = n, theta = theta, k = k.vec))
#' 
#' }
#' 
#' ## Now we reprodue Figure 4.1 in John Wakeley (2008): 
#' ## 'Coalescent Theory: An Indtroduction' by using the package plot3D.
#' library(plot3D)
#' hist3D(x=k.vec, y=1:20, z=Res.mat, col = "grey", border = "black",
#'        xlab = "k", ylab = "n", zlab = "P(S=k)",
#'        main = "The probability function of the number of segregating sites",
#'        sub = expression(paste("The mutation parameter is ", theta,"= 2")),
#'        cex.main = 0.9, colkey = F, zlim = c(0,0.4))
#' 
#' @export
dSegregatingSites <- function(n, theta, k, plot =FALSE){
  
  if(n < 1) stop("Invalid sample size. n has to be positive!")
  if(sum(k<0)>0) stop("Invalid vector of quantiles. k has to be nonnegative!")
  
  if(n==1){
    
    res <- replicate(k,0)
    
  }else if(n==2){
    
    ## We define the reward transformed subtransition probability matrix
    P.mat <- theta/(theta+1)
    p.vec <- 1/(theta+1)
    ## We store the results in a matrix
    ProbSegSites <- (P.mat^k)*p.vec
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
    
    obj <- contphasetype(pi.vec, Tmat)
    
    ## Now we can compute the distribution of the number of
    ## segregating sites by using the descretization:
    newobj <- discretization(obj, a=NULL, lambda=theta/2)
    
    res <- NULL
    for (i in k) {
      
      res[which(k==i)] <- pi.vec%*%(newobj$P.mat%^%i)%*%(1-rowSums(newobj$P.mat))
    }
  }
  
  if(plot){
    
    plot(x=k, y=res, type = "l", col = "grey",
         xlab = "k", ylab = expression(paste("P(", S["Total"], "=k)")),
         main = "The density function of the number of segregating sites",
         sub = paste("The mutation parameter is equal to", theta),
         cex.main = 0.9, ylim = c(0,max(res)*1.2))
  }
  return(res)
}







