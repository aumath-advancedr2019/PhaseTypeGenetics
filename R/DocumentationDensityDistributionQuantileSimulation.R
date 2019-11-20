#' The phase-type distribution
#'
#' Density, distribution function, quantile function and simulations for
#' the phase-type distribution with initial distribution equal to
#' \code{initDist} and sub-transition/sub-intensity matrix equal to
#' \code{P_Mat}/\code{T_Mat}.
#'
#' In the discrete case, the phase-type distribution has density
#' \deqn{f(x) = initDist (P_Mat ^ (x-1))t + (x=1)(1-sum(initDist)), }
#' for integers \eqn{x \ge 1}, where \code{initDist} is the initial distribution, \code{P_Mat} is the sub-transition
#' probability matrix and \code{t = (I-P)e}. Furthermore, the distribution function
#' is given by
#' \deqn{F(x) = 1- initDist (P_Mat ^ x) e + (x \ge 1)(1-sum(initDist)). }
#' If the quantile \eqn{x} is a
#' real number, the function will round the number down in order to obtain a
#' natural number.
#' In the continuous case, the phase-type distribution has density
#' \deqn{f(x) = initDist expm(x T_Mat) t, for x \ge 0 , }
#' where \code{initDist} is the initial distribution, \code{T_Mat} is the sub-intensity
#' rate matrix and \code{t = -Te}. Furthermore, the distribution function
#' is given by
#' \deqn{F(x) = 1- initDist expm(x T_Mat) e, for x \ge 0. }
#'
#' @param object an object for which the density, distribution function,
#' quantile function or random generation should be computed. To be able to use
#' these function,the object has to be of
#' class \code{discphasetype} or \code{contphasetype}.
#' @param x a single non-negative quantile or a vector of non-negative quantiles.
#' @param p a single probability or a vector of probabilities.
#' @param n the number of observations (n>=1).
#'
#' @return \code{dphasetype} gives the density, \code{pphasetype}
#' gives the distribution function, \code{qphasetype} gives the quantile
#' function, and \code{rphasetype} simulates from the distribution.
#' The length of the output is equal to the length of the input, i.e. for
#' \code{dphasetype} and \code{pphasetype} the length of the output is equal to
#' the length of the quantile vector \code{x}, for \code{qphasetype} the output is of
#' the same length as the input vector \code{p}, and \code{rphasetype} produces
#' an output of length \eqn{n}.
#'
#' @source Mogens Bladt and Bo Friis Nielsen (2017):
#' \emph{ Matrix-Exponential Distributions in Applied Probability}.
#' Probability Theory and Stochastic Modelling (Springer), Volume 81.
#'
#' @seealso \code{\link{dnorm}}, \code{\link{dt}},
#' \code{\link{dexp}}.
#'
#' @importFrom expm %^%
#'
#' @examples
#'
#' ## We reproduce Figure 3.4 in John Wakeley (2009):
#' ## "Coalescent Theory: An Introduction",
#' ## Roberts and Company Publishers, Colorado.
#'
#' x <- seq(0,4, by=0.1)
#' dist <- matrix(nrow = 6, ncol = length(x))
#'
#' dist[2,] <- dphasetype(T_MRCA$n5, x)
#' dist[3,] <- dphasetype(T_MRCA$n10, x)
#' dist[4,] <- dphasetype(T_MRCA$n20, x)
#' dist[5,] <- dphasetype(T_MRCA$n50, x)
#' dist[6,] <- dphasetype(T_MRCA$n100, x)
#'
#' ## For n=2, the initial distribution is initDist = 1 and the sub-transition probability
#' ## matrix is T_Mat = -1, hence the distribution is given by
#' dist[1,] <- exp(-x)
#'
#' plot(x, dist[1,], type = "l", main = expression(paste("The distribution of ", T["MRCA"],
#'      " for n=2,5,10,20,50,100")), cex.main = 0.9, xlab = "x", ylab = expression(f[T[MRCA]](x)),
#'      xlim = c(0,4), ylim = c(0,1), frame.plot = FALSE)
#' points(x, dist[2,], type = "l")
#' points(x, dist[3,], type = "l")
#' points(x, dist[4,], type = "l")
#' points(x, dist[5,], type = "l")
#' points(x, dist[6,], type = "l")
#'
#'
#' x <- seq(0,15, by=0.1)
#' dist <- matrix(nrow = 6, ncol = length(x))
#'
#' dist[2,] <- dphasetype(T_Total$n5,x)
#' dist[3,] <- dphasetype(T_Total$n10,x)
#' dist[4,] <- dphasetype(T_Total$n20,x)
#' dist[5,] <- dphasetype(T_Total$n50,x)
#' dist[6,] <- dphasetype(T_Total$n100,x)
#'
#' ## For n=2, the initial distribution is initDist = 1 and the sub-transition probability
#' ## matrix is T_Mat = -1/2, hence the distribution is given by
#' dist[1,] <- exp(-x/2)/2
#'
#' plot(x, dist[1,], type = "l", main = expression(paste("The distribution of ", T["Total"],
#'     " for n=2,5,10,20,50,100")), cex.main = 0.9, xlab = "x", ylab = expression(f[T[Total]](x)),
#'     xlim = c(0,15), ylim = c(0,0.5), frame.plot = FALSE)
#' points(x, dist[2,], type = "l")
#' points(x, dist[3,], type = "l")
#' points(x, dist[4,], type = "l")
#' points(x, dist[5,], type = "l")
#' points(x, dist[6,], type = "l")
#'
#' ## Simulating ten total branch lengths
#' ## for a sample of size 5
#' rphasetype(T_Total$n5, n=10)
#'
#' @export
dphasetype <- function(object, x){

  UseMethod("dphasetype")
}

#' @export
#' @importFrom expm %^%
dphasetype.discphasetype <- function(object, x){

  res <- NULL
  for(l in x){

    if(l<1 | l%%1!=0){

      warning("One or more quantiles are less than 1 or not natural numbers.\n
              The corresponding probabilities are set to 0.")
      res[which(x==l)] <- 0
    }else{

      res[which(x==l)] <- sum(object$initDist%*%(object$P_Mat %^%(l-1))%*%(diag(1, nrow = nrow(object$P_Mat))-object$P_Mat)) +
                          (l==1)*(1-sum(object$initDist))
    }
  }
 return(res)
}

#' @export
dphasetype.contphasetype <- function(object, x){

  res <- NULL
  for(l in x){

    if(l<0){

      warning("One or more quantiles are negative. The corresponding probabilities are set to 0.")
      res[which(x==l)] <- 0

    }else{

      res[which(x==l)] <- -sum(object$initDist %*% expm::expm(l * object$T_Mat) %*% object$T_Mat)
    }
  }
  return(res)
}

#' @rdname dphasetype
#' @export
pphasetype <- function(object, x){

  UseMethod("pphasetype")
}

#' @export
#' @importFrom expm %^%
pphasetype.discphasetype <- function(object, x){

  res <- NULL
  for(l in x){

    if(l<1 | l%%1!=0){

      warning("One or more quantiles are less than 1 or not natural numbers.\n
              The corresponding probabilities are set to 0.")
      res[which(x==l)] <- 0

    }else{

      res[which(x==l)] <- 1 - sum(object$initDist%*%(object$P_Mat %^% l))
    }
  }
  return(res)
}

#' @export
pphasetype.contphasetype <- function(object, x){

  res <- NULL
  for(l in x){

    if(l<0){

      warning("One or more quantiles are negative. The corresponding probabilities are set to 0.")
      res[which(x==l)] <- 0

    }else{

      res[which(x==l)] <- 1 - sum(object$initDist %*% expm::expm(l * object$T_Mat))
    }
  }
  return(res)
}

#' @rdname dphasetype
#' @export
qphasetype <- function(object, p){

  UseMethod("qphasetype")
}

#' @export
qphasetype.discphasetype <- function(object, p){

  if( sum(p<0) >0 | sum(p>1) >0 ){

    stop("Not a valid probability vector. One or more entries are less than 0 or bigger than 1.")

    }else{

    res <- NULL
    for(l in p){

      res[which(p==l)] <- stats::uniroot(function(y) pphasetype(object = object, x = y)- l, c(0, 400))$root[1]
    }
  }
  return(round(res))
}

#' @export
qphasetype.contphasetype <- function(object, p){

  if( sum(p<0) >0 | sum(p>1) >0 ){

    stop("Not a valid probability vector. One or more entries are less than 0 or bigger than 1.")
  }else{

    res <- NULL
    for(l in p){

      res[which(p==l)] <- stats::uniroot(function(y) pphasetype(object = object, x = y)-l, c(0, 400))$root[1]
    }
  }
  return(res)
}

#' @rdname dphasetype
#' @export
rphasetype <- function(object, n){

  UseMethod("rphasetype")
}

#' @export
rphasetype.discphasetype <- function(object, n){

  if(n != floor(abs(n))){

    warning(paste("The proviede number n is either negative or not a natural number.\n
                   The function will generate ", floor(abs(n)), "simulations instead."))
  }

  n=floor(abs(n))

  ## Extracting the initial distribution
  ## and the subtransition matrix
  initDist = object$initDist
  P_Mat = object$P_Mat

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
    u <- stats::runif(n)
  }

  # We make the rows of the transition matrix corresponding to
  # the transient states from the subtransition matrix
  TransMat <- cbind(P_Mat, 1-rowSums(P_Mat))

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
rphasetype.contphasetype <- function(object, n){

  if(n != floor(abs(n))){

  warning(paste("The proviede number n is either negative or not a natural number.\n
                 The function will generate ", floor(abs(n)), "simulations instead."))
  }
  n=floor(abs(n))

  ## Extracting the initial distribution
  ## and the subtransition matrix
  initDist = object$initDist
  T_Mat = object$T_Mat

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
    u <- stats::runif(n)
  }

  # We make the rows of the intensity matrix corresponding to the transient states from
  # the subintensity matrix
  IntenseMat <- cbind(T_Mat, -rowSums(T_Mat))

  # Now we simulate the Markov Jump Process until absorption in 'p+1'
  for(i in 1:n){
    if(u[i] <= defect){
      next()
    }
    else{
      initState <- sample(p, size = 1, prob = initDist)
      x <- initState
      while(x != p + 1){

        tau[i] <- tau[i] + stats::rexp(1, rate = -IntenseMat[x,x])
        x <- sample((1:(p+1))[-x], size = 1, prob = IntenseMat[x,-x] )

      }
    }
  }
  return(tau)
}

