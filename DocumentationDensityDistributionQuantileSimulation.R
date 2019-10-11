#' The phase-type distribution
#'
#' Density, distribution function, quantile function and simulations for
#' the phase-type distribution with initial distribution equal to
#' \code{initDist} and subtransition/subintensity matrix equal to
#' \code{P.mat}/\code{T.mat}.
#'
#' In the discrete case, the phase-type distribution has density
#' f(x) = initDist %*% (P.mat %^% (x-1)) %*% t + (x==1)*(1-sum(initDist)), for integers x>=1 
#' where initDist is the initial distribution, P.mat is the subtransition
#' probability matrix and t = (I-P)e. Furthermore, the distribution function
#' is given by
#' F(x) = 1- initDist %*% (P.mat %^% x) %*% e + (x>=1)*(1-sum(initDist)). If the quantile \code{x} is a
#' real number, the function will round the number down in order to obtain a
#' natural number.
#' In the continuous case, the phase-type distribution has density
#' f(x) = initDist %*% expm(x * T.mat) %*% t, for x>=0
#' where initDist is the initial distribution, T.mat is the subintensity
#' rate matrix and t = -Te. Furthermore, the distribution function
#' is given by
#' F(x) = 1- initDist %*% expm(x * T.mat) %*% e, for x>=0.
#'
#' @param object an object for which the density, distribution function,
#' quantile function or random generation should be computed. To be able to use
#' these function,the object has to be of
#' class \code{discphasetype} or \code{contphasetype}.
#' @param x a positive quantile
#' @param p a probability
#' @param n number of observations
#'
#' @return \code{dphasetype} gives the density, \code{pphasetype}
#' gives the distribution function, \code{qphasetype} gives the quantile
#' function, and \code{rphasetype} simulates from the distribution.
#' The length of the output is 1, except for \code{rphasetype}, which produces
#' an output of length \code{n}.
#'
#' @source Mogens Bladt and Bo Friis Nielsen (2017):
#' \emph{ Matrix-Exponential Distributions in Applied Probability}.
#' Probability Theory and Stochastic Modelling (Springer), Volume 81.
#'
#' @seealso \code{\link{dnorm}}, \code{\link{dt}},
#' \code{\link{dexp}}.
#'
#' @examples
#'
#'
#' @export

dphasetype <- function(...){

  UseMethod("dphasetype")
}

#' @describeIn dphasetype Density of discrete phase-type distributed variable
dphasetype.discphasetype <- function(object,x){

  if(x<1 | !is.integer(x)){
    return(0)
  }else{

  initDist = object$initDist
  P.mat = object$P.mat
  }

 return(sum(initDist%*%(P.mat %^%(x-1))%*%(diag(1, nrow = nrow(P.mat))-P.mat)) + (x==1)*(1-sum(iniDist)))
}

#' @describeIn dphasetype Density of continuous phase-type distributed variable
dphasetype.contphasetype <- function(object, x){

  if(x<0){
    
    return(0)
    
  }else{

  initDist = object$initDist
  T.mat = object$T.mat
  }

  return(-sum(initDist %*% expm(x * T.mat) %*% (diag(1, nrow = nrow(T.mat))-T.mat)))
}

#' @describeIn dphasetype generic
pphasetype <- function(...){

  UseMethod("pphasetype")
}

#' @describeIn dphasetype distribution function of discrete phase-type distributed variable
pphasetype.discphasetype <- function(object,x){

  if(x<0){
    
      return(0)
    
  }else{
    
  x <- floor(x)
  initDist = object$initDist
  P.mat = object$P.mat
  }

  return(1 - sum(initDist%*%(P.mat %^% x)))
}

#' @describeIn dphasetype distribution function of continuous phase-type distributed variable
pphasetype.contphasetype <- function(object, x){

  if(x<0){

    return(0)
    
  }else{

  initDist = object$initDist
  T.mat = object$T.mat
  }

  return(1 - sum(initDist %*% expm(x * T.mat)))
}

#' @describeIn dphasetype generic
qphasetype <- function(...){

  UseMethod("qphasetype")
}

#' @describeIn dphasetype Quantile function of discrete phase-type distributed variable
qphasetype.discphasetype <- function(object, p){

  if( p<0 | p>1 ){

    stop("Not a valid probability")
  }else{

  m <- uniroot(function(y) pphasetype(object = object, x = y)- p, c(0, 400))
  }
  return(round(m$root[1]))
}

#' @describeIn dphasetype Quantile function of continuous phase-type distributed variable
qphasetype.contphasetype <- function(object, p){

  if( p<0 | p>1 ){

    stop("Not a valid probability")
  }else{

  m <- uniroot(function(y) pphasetype(object = object, x = y)-p, c(0, 400))
  }
  return(m$root[1])

}

#' @describeIn dphasetype generic
rphasetype <- function(...){

  UseMethod("rphasetype")
}

#' @describeIn dphasetype Simulating from a discrete phase-type distribution
rphasetype.discphasetype <- function(object, n){

  n=floor(abs(n))

  ## Extracting the initial distribution
  ## and the subtransition matrix
  initDist = object$initDist
  P.mat = object$P.mat

  # Calculate the number of transient states
  p <- length(initDist)

  # Calculate the defect
  defect <- 1 - sum(initDist)

  # Initialize the vector with 1's because if the Markov Chain
  # is absorbed immediately then tau would be 1
  tau <- rep(1,n)

  # If the defect is positive immediate absorption needs to be a possibility
  # for that we generate n uniform(0,1) variables
  if(defect < 1){
    u <- runif(n)
  }

  # We make the rows of the transition matrix corresponding to
  # the transient states from the subtransition matrix
  TransMat <- cbind(P.mat, 1-rowSums(P.mat))

  # Now we simulate the Markov Chain until absorption in 'p+1'
  for(i in 1:n){
    if(u[i] <= defect){
      next()
    }
    else{
      initState <- sample(p, size = 1, prob = initDist)
      x <- initState
      while(x != (p + 1)){

        tau[i] <- tau[i] + 1
        x <- sample(p + 1, size = 1, prob = TransMat[x,] )
      }
    }
  }
  return(tau)
}

#' @describeIn dphasetype Simulating from a continuous phase-type distribution
rphasetype.contphasetype <- function(object,n){

  n=floor(abs(n))

  ## Extracting the initial distribution
  ## and the subtransition matrix
  initDist = object$initDist
  T.mat = object$T.mat

  # Calculate the number of transient states
  p <- length(initDist)

  # Calculate the defect
  defect <- 1 - sum(initDist)

  # Initialize the vector with 0's because if the Markov Chain is absorbed immediately
  # then tau would be 0
  tau <- rep(0,n)

  # If the defect is positive immediate absorption needs to be a possibility
  # for that we generate n uniform(0,1) variables
  if(defect < 1){
    u <- runif(n)
  }

  # We make the rows of the intensity matrix corresponding to the transient states from
  # the subintensity matrix
  IntenseMat <- cbind(T.mat, 1-rowSums(T.mat))

  # Now we simulate the Markov Jump Process until absorption in 'p+1'
  for(i in 1:n){
    if(u[i] <= defect){
      next()
    }
    else{
      initState <- sample(p, size = 1, prob = initDist)
      x <- initState
      while(x != p + 1){

        tau[i] <- tau[i] + rexp(1, rate = -IntenseMat[x,x])
        x <- sample((1:(p+1))[-x], size = 1, prob = IntenseMat[x,-x] )

      }
    }
  }
  return(tau)
}

