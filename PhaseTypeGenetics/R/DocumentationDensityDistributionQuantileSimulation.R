#' The phase-type distribution
#'
#' Density, distribution function, quantile function and simulations for
#' the phase-type distribution with initial distribution equal to
#' \code{initDist} and subtransition/subintensity matrix equal to
#' \code{P.mat}/\code{T.mat}.
#'
#' In the discrete case, the phase-type distribution has density
#' \deqn{f(x) = initDist (P.mat ^ (x-1))t + (x=1)(1-sum(initDist)), }
#' for integers \eqn{x \ge 1}, where \code{initDist} is the initial distribution, \code{P.mat} is the subtransition
#' probability matrix and \code{t = (I-P)e}. Furthermore, the distribution function
#' is given by
#' \deqn{F(x) = 1- initDist (P.mat ^ x) e + (x \ge 1)(1-sum(initDist)). }
#' If the quantile \eqn{x} is a
#' real number, the function will round the number down in order to obtain a
#' natural number.
#' In the continuous case, the phase-type distribution has density
#' \deqn{f(x) = initDist expm(x T.mat) t, for x \ge 0 , }
#' where \code{initDist} is the initial distribution, \code{T.mat} is the subintensity
#' rate matrix and \code{t = -Te}. Furthermore, the distribution function
#' is given by
#' \deqn{F(x) = 1- initDist expm(x T.mat) e, for x \ge 0. }
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
#' an output of length \eqn{n}.
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
#' ## We reproduce Figure 3.4 in John Wakeley (2009):
#' ## "Coalescent Theory: An Introduction",
#' ## Roberts and Company Publishers, Colorado.
#'
#' x.vec <- seq(0,4, by=0.1)
#' dist <- matrix(nrow = 6, ncol = length(x.vec))
#' for(x in x.vec){
#'
#'  dist[2,which(x.vec ==x)] <- dphasetype(T_MRCA$n5,x)
#'  dist[3,which(x.vec ==x)] <- dphasetype(T_MRCA$n10,x)
#'  dist[4,which(x.vec ==x)] <- dphasetype(T_MRCA$n20,x)
#'  dist[5,which(x.vec ==x)] <- dphasetype(T_MRCA$n50,x)
#'  dist[6,which(x.vec ==x)] <- dphasetype(T_MRCA$n100,x)
#' }
#'
#' ## For n=2, the initial distribution is initDist = 1 and the transition probability
#' ## matrix is T.mat = -1, hence the distribution is given by
#' dist[1,] <- exp(-x.vec)
#'
#' plot(x.vec, dist[1,], type = "l", main = expression(paste("The distribution of ", T["MRCA"],
#'      " for n=2,5,10,20,50,100")), cex.main = 0.9, xlab = "x", ylab = expression(f[T[MRCA]](x)),
#'      xlim = c(0,4), ylim = c(0,1), frame.plot = F)
#' points(x.vec, dist[2,], type = "l")
#' points(x.vec, dist[3,], type = "l")
#' points(x.vec, dist[4,], type = "l")
#' points(x.vec, dist[5,], type = "l")
#' points(x.vec, dist[6,], type = "l")
#'
#'
#' x.vec <- seq(0,15, by=0.1)
#' dist <- matrix(,nrow = 6, ncol = length(x.vec))
#' for(x in x.vec){
#'
#'  dist[2,which(x.vec ==x)] <- dphasetype(T_Total$n5,x)
#'  dist[3,which(x.vec ==x)] <- dphasetype(T_Total$n10,x)
#'  dist[4,which(x.vec ==x)] <- dphasetype(T_Total$n20,x)
#'  dist[5,which(x.vec ==x)] <- dphasetype(T_Total$n50,x)
#'  dist[6,which(x.vec ==x)] <- dphasetype(T_Total$n100,x)
#' }
#'
#' ## For n=2, the initial distribution is initDist = 1 and the transition probability
#' ## matrix is T.mat = -1/2, hence the distribution is given by
#' dist[1,] <- exp(-x.vec/2)/2
#'
#' plot(x.vec, dist[1,], type = "l", main = expression(paste("The distribution of ", T["Total"],
#'     " for n=2,5,10,20,50,100")), cex.main = 0.9, xlab = "x", ylab = expression(f[T[Total]](x)),
#'     xlim = c(0,15), ylim = c(0,0.5), frame.plot = F)
#' points(x.vec, dist[2,], type = "l")
#' points(x.vec, dist[3,], type = "l")
#' points(x.vec, dist[4,], type = "l")
#' points(x.vec, dist[5,], type = "l")
#' points(x.vec, dist[6,], type = "l")
#'
#' ## Simulating ten total branch lengths
#' ## for a sample of size 5
#' rphasetype(T_Total$n5, n=10)
#'
#' @export
dphasetype <- function(object,x){

  UseMethod("dphasetype")
}

#' @export
dphasetype.discphasetype <- function(object,x){

  if(x<1 | !is.integer(x)){
    return(0)
  }else{

  initDist = object$initDist
  P.mat = object$P.mat
  }

 return(sum(initDist%*%(P.mat %^%(x-1))%*%(diag(1, nrow = nrow(P.mat))-P.mat)) + (x==1)*(1-sum(iniDist)))
}

#' @export
dphasetype.contphasetype <- function(object, x){

  if(x<0){

    return(0)

  }else{

  initDist = object$initDist
  T.mat = object$T.mat
  }

  return(-sum(initDist %*% expm(x * T.mat) %*% T.mat))
}

#' @rdname dphasetype
#' @export
pphasetype <- function(object,x){

  UseMethod("pphasetype")
}

#' @export
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

#' @export
pphasetype.contphasetype <- function(object, x){

  if(x<0){

    return(0)

  }else{

  initDist = object$initDist
  T.mat = object$T.mat
  }

  return(1 - sum(initDist %*% expm(x * T.mat)))
}

#' @rdname dphasetype
#' @export
qphasetype <- function(object, p){

  UseMethod("qphasetype")
}

#' @export
qphasetype.discphasetype <- function(object, p){

  if( p<0 | p>1 ){

    stop("Not a valid probability")
  }else{

  m <- uniroot(function(y) pphasetype(object = object, x = y)- p, c(0, 400))
  }
  return(round(m$root[1]))
}

#' @export
qphasetype.contphasetype <- function(object, p){

  if( p<0 | p>1 ){

    stop("Not a valid probability")
  }else{

  m <- uniroot(function(y) pphasetype(object = object, x = y)-p, c(0, 400))
  }
  return(m$root[1])

}

#' @rdname dphasetype
#' @export
rphasetype <- function(object, n){

  UseMethod("rphasetype")
}

#' @export
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

#' @export
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

