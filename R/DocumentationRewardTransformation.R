#' Reward transformation
#'
#' Reward transforming a continuous phase-type distribution with initial distribution
#' \code{initDist} and sub-intensity rate matrix \code{T_Mat}.
#'
#' It is possible to assign a linear (non-negative) reward to each of the transient states
#' of the Markov jump process underlying the phase-type distribution under consideration.
#' More precisely, a reward is earned in each transient state and that reward is
#' proportional to the time spent in the state. By assigning rewards, the original
#' phase-type distribution is transformed into another continuous phase-type distribution.
#' More details on the reward transformation can be found in \code{vignette("PhaseTypeGenetics")} or in Mogens Bladt and Bo Friis Nielsen (2017):
#' \emph{ Matrix-Exponential Distributions in Applied Probability}.
#'
#' @param object a continuous phase-type distributed object of class \code{contphasetype}.
#' @param rewards a non-negative reward vector. The length of the reward vector
#' should be equal to the number of rows of the sub-intensity rate matrix \code{T_Mat}.
#'
#' @return The function returns the reward transformed phase-type distribution,
#' which is again an object of type \code{contphasetype}.
#'
#' @source Mogens Bladt and Bo Friis Nielsen (2017):
#' \emph{ Matrix-Exponential Distributions in Applied Probability}.
#' Probability Theory and Stochastic Modelling (Springer), Volume 81.
#'
#'
#' @examples
#'
#' ## The time to the most recent common ancestor (T_MRCA)
#' ## is phase-type distributed with initial distribution
#' initDist <- c(1,0,0,0)
#' ## and sub-intensity rate matrix
#' T_Mat <- matrix(c(-6,6,0,0,
#'                   0,-3,2,1,
#'                   0,0,-1,0,
#'                   0,0,0,-1), nrow = 4, byrow = TRUE)
#' ## Defining an object of type "contphasetype"
#' TMRCA <- contphasetype(initDist, T_Mat)
#' ## In order to obtain the distribution of the total
#' ## length of all branches giving rise to singletons,
#' ## we have to give the following rewards to the
#' ## different states
#' r.vec <- c(4,2,1,0)
#' ## Hence,
#' RewTransDistribution(TMRCA,r.vec)
#'
#' @export

RewTransDistribution <- function(object, rewards){

  if(class(object) != "contphasetype") stop("Invalid object! The object has to be of type 'contphasetype'.")

  if(sum(rewards < 0) > 0 | sum(rewards)==0 | length(rewards)!=nrow(object$T_Mat)){

    stop("The reward vector has to be nonnegative and contain at least one positve number.\n
         Furthermore, the length of the reward vector has to be equal to the number of rows of the sub-intensity rate matrix.")

  }else{

    #We extract the initial distribution and subintensity matrix from the object
    initDist <- object$initDist
    T_Mat <- object$T_Mat

    ## We define the sets S^+ and S^0
    S.plus <- which(rewards>0)
    S.zero <- which(rewards==0)
    ## as well as the matrix Q
    Q.mat <- -T_Mat/diag(T_Mat)
    diag(Q.mat) <- 0

    ## We rearrange the states according to the
    ## sets S^+ and S^0:
    Q.matpp <- Q.mat[S.plus,S.plus]
    Q.matp0 <- Q.mat[S.plus,S.zero]
    Q.mat0p <- Q.mat[S.zero,S.plus]
    Q.mat00 <- Q.mat[S.zero,S.zero]

    ## Now we can define the transition matrix
    ## for the new markov chain
    if(length(S.zero)==1){

      P.mat <- Q.matpp+(1-Q.mat00)^(-1)*Q.matp0%*%t(Q.mat0p)

    }else if(length(S.zero)==0){
      P.mat <- Q.matpp

    }else{
      P.mat <- Q.matpp+Q.matp0%*%solve(diag(1,nrow = nrow(Q.mat00))-
                                         Q.mat00)%*%Q.mat0p
    }

    ## We define the new exit vector
    p.vec <- 1-rowSums(P.mat)

    ## We also split the original initial distribution
    ## into pi = (pi^+,pi^0)
    pi.vecp <- initDist[S.plus]
    pi.vec0 <- initDist[S.zero]

    ## Then, the initial distribution of the new
    ## Markov chain is given by
    if(length(S.zero)==1){

      newInitDist <- pi.vecp + (1-Q.mat00)*pi.vec0%*%t(Q.mat0p)

    }else if(length(S.zero)==0){

      newInitDist <- pi.vecp

    }else{
      newInitDist <- pi.vecp + pi.vec0%*%solve(diag(1,nrow = nrow(Q.mat00))
                                               -Q.mat00)%*%Q.mat0p
    }

    ## Now we can define the subintensity matrix
    ## newT_Mat as
    if(length(S.plus)==1){

      newT_Mat <- -T_Mat[S.plus,S.plus]/rewards[S.plus]*P.mat
      exitrate <- -T_Mat[S.plus,S.plus]/rewards[S.plus]*p.vec
    }else{

      newT_Mat <- -diag(T_Mat[S.plus,S.plus])/rewards[S.plus]*P.mat
      exitrate <- -diag(T_Mat[S.plus,S.plus])/rewards[S.plus]*p.vec
    }


    diag(newT_Mat) <- -rowSums(newT_Mat)-exitrate

  }
  return(contphasetype(initDist = newInitDist, T_Mat = newT_Mat))
}
