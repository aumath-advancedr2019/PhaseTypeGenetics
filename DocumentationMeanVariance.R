#' Mean and variance of phase-type distributions
#'
#' Mean and variance for both the discrete and continuous 
#' phase-type distribution with initial distribution equal to
#' \code{initDist} and subtransition/subintensity matrix equal to
#' \code{P.mat}/\code{T.mat}.
#'
#' In the discrete case, the phase-type distribution has mean
#' \deqn{E[\tau] = initDist (I-P.mat)^{-1} e + 1 - initDist e},
#' where \code{initDist} is the initial distribution, \code{P.mat} is the subtransition
#' probability matrix and \code{e} is the vector having one in each entry. 
#' Furthermore, the variance can be calculated as
#' \deqn{Var[\tau] = E[\tau(\tau-1)] + E[\tau] - E[\tau]^2}, 
#' where 
#' \deqn{E[\tau(\tau-1)] = 2 initDist P.mat (I-P.mat)^{-2} e + 1 - initDist e}.
#' In the continuous case, the phase-type distribution has mean
#' \deqn{E[\tau] = initDist (-T.mat)^{-1} e},
#' where \code{initDist} is the initial distribution and \code{T.mat} is the subintensity
#' rate matrix. Furthermore, the variance can be calculated in the usual way
#' \deqn{Var[tau] = E[tau^2] - E[tau]^2},
#' where 
#' \deqn{E[\tau^2] = 2 initDist (-T.mat)^{-2} e}.
#'
#' @param object an object for which the mean or variance should be computed.
#' To be able to use these function,the object has to be of
#' class \code{discphasetype} or \code{contphasetype}.
#'
#' @return \code{mean} gives the mean and \code{var} gives the 
#' variance of the phase-type distribution. The length of the output is 1.
#'
#' @source Mogens Bladt and Bo Friis Nielsen (2017):
#' \emph{ Matrix-Exponential Distributions in Applied Probability}.
#' Probability Theory and Stochastic Modelling (Springer), Volume 81.
#'
#' @seealso \code{\link{mean}}, \code{\link{var}}.
#'
#' @examples
#'
#' ## We reproduce Figure 3.3 in John Wakeley (2009): 
#' ## "Coalescent Theory: An Introduction", 
#' ## Roberts and Company Publishers, Colorado.
#'
#' ## We define vectors holding the means and variances
#' VecOfMeansMRCA <- replicate(20,0)
#' VecOfVarsMRCA <- replicate(20,0)
#' VecOfMeansTotal <- replicate(20,0)
#' VecOfVarsTotal <- replicate(20,0)
#'
#' ## For n=2, we have that the initial distribution is initDist = 1 and
#' ## the transition probability matrix is T.mat = -1 for T_MRCA and
#' ## T.mat = -1/2 for T_Total,
#' ## hence the mean is given by
#' VecOfMeansMRCA[2] <- 1 ## -pi*T^(-1)* e = -1*(-1)*1 = 1
#' VecOfMeansTotal[2] <- 2 ## -pi*T^(-1)* e = -1*(-1/2)^(-1)*1 = 2
#' ## and the variance is
#' VecOfVarsMRCA[2] <- 1 ## 2-1^2
#' VecOfVarsTotal[2] <- 4 ## 8-2^2
#' ## as variance= E[X^2] - E[X]^2.
#' 
#' # For n=3, we have that the initial distibrution is
#' initDist = c(1,0)
#' ## and the transition probability matrices are
#' T.matMRCA = matrix(c(-3,3,0,-1), nrow = 2, byrow = TRUE)
#' T.matTotal = matrix(c(-2,2,0,-1), nrow = 2, byrow = TRUE)/2
#' ## for T_MRCA and T_Total, respectively. 
#' ## Defining two objects of class "contphasetype"
#' T_MRCA <- contphasetype(initDist, T.matMRCA)
#' T_Total <- contphasetype(initDist, T.matTotal)
#' ## Hence the means are given by
#' VecOfMeansMRCA[3] <- mean(T_MRCA)
#' VecOfMeansTotal[3] <- mean(T_Total)
#' ## and the variances are
#' VecOfVarsMRCA[3] <- var(T_MRCA)
#' VecOfVarsTotal[3] <-var(T_Total)
#' 
#' for (n in 4:20) {
#' 
#'  ## The initial distribution
#'  initDist <- c(1,replicate(n-2,0))
#'  ## The subintensity rate matrix
#'  T.mat <- diag(choose(n:3,2))
#'  T.mat <- cbind(replicate(n-2,0),T.mat)
#'  T.mat <- rbind(T.mat, replicate(n-1,0))
#'  diag(T.mat) <- -choose(n:2,2)
#'  ## Define an object of class "contphasetype"
#'  obj <- contphasetype(initDist,T.mat)
#'  ## Compute the mean and variance 
#'  VecOfMeansMRCA[n] <- mean(obj)
#'  VecOfVarsMRCA[n] <- var(obj)
#'  
#'  ## For T_total, we compute the same numbers
#'  ## The subintensity rate matrix
#'  T.mat <- diag((n-1):2)
#'  T.mat <- cbind(replicate(n-2,0),T.mat)
#'  T.mat <- rbind(T.mat, replicate(n-1,0))
#'  diag(T.mat) <- -((n-1):1)
#'  T.mat <- 1/2*T.mat
#'  ## Define an object of class "contphasetype"
#'  obj <- contphasetype(initDist,T.mat)
#'  ## Compute the mean and variance
#'  VecOfMeansTotal[n] <- mean(obj)
#'  VecOfVarsTotal[n] <- var(obj)
#' }
#'
#' ## Plotting the means
#' plot(x = 1:20, VecOfMeansMRCA, type = "l", main = expression(paste("The dependence of ",E(T[MRCA]),"
#'     and ", E(T[Total]), " on the sample size")), cex.main = 0.9, xlab = "n",
#'     ylab = "Expectation",
#'     xlim = c(0,25), ylim = c(0,8), frame.plot = F)
#' points(x= 1:20, VecOfMeansTotal, type = "l")
#'
#' text(23,tail(VecOfMeansMRCA, n=1),labels = expression(E(T[MRCA])))
#' text(23,tail(VecOfMeansTotal, n=1),labels = expression(E(T[Total])))
#'
#' ## And plotting the variances
#' plot(x = 1:20, VecOfVarsMRCA, type = "l", main = expression(paste("The dependence of ",Var(T[MRCA]), 
#'      " and ", Var(T[Total]), " on the sample size")), cex.main = 0.9, xlab = "n", ylab = "Variance",
#'     xlim = c(0,25), ylim = c(0,7), frame.plot = F)
#' points(x= 1:20, VecOfVarsTotal, type = "l")
#'
#' text(23,tail(VecOfVarsMRCA, n=1),labels = expression(Var(T[MRCA])))
#' text(23,tail(VecOfVarsTotal, n=1),labels = expression(Var(T[Total])))
#'
mean.discphasetype <- function(object){
  
  initDist <- object$initDist
  P.mat <- object$P.mat
  
  return(sum(initDist%*%solve(diag(x=1, nrow = nrow(P.mat))-P.mat)) + 1 - sum(initDist))
  
}

#' @rdname mean.discphasetype
mean.contphasetype <- function(object){
  
  initDist = object$initDist
  T.mat= object$T.mat
  
  return(sum(initDist%*%solve(-T.mat)))
}

#' @rdname mean.discphasetype
var <- function(...){
  
  UseMethod("var")
}

#' @rdname mean.discphasetype
var.default <- function(x,...){
  
  var(x,...)
}

#' @rdname mean.discphasetype
var.discphasetype <- function(object){
  
  initDist = object$initDist
  P.mat = object$P.mat
  defect <- 1 - sum(initDist)
  
  secondMoment <- 2*sum(initDist%*%P.mat%*%solve((diag(x=1, nrow = nrow(P.mat))-P.mat)%^%2))
  firstmoment <- sum(initDist%*%solve(diag(x=1, nrow = nrow(P.mat))-P.mat)) + defect
  return(secondMoment + firstmoment - firstmoment^2)
}

#' @rdname mean.discphasetype
var.contphasetype <- function(object){
  
  initDist = object$initDist
  T.mat = object$T.mat
  
  secondMoment <- 2*sum(initDist%*%solve(T.mat%^%2))
  firstmoment <- sum(initDist%*%solve(-T.mat))
 
  return(secondMoment-firstmoment^2) 
}
