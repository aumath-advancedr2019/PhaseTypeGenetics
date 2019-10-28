#' Block counting process
#'
#' Computing the state space and the corresponding rate matrix for 
#' the block counting process for a given sample size \eqn{n} in 
#' the standard coalescent model.
#' 
#' For a given sample size \eqn{n}, one can have one or more possible 
#' coalescent trees. Each coalescent event in these trees correspond 
#' to a state of the block counting process. Furthermore, each state is
#' represented by a \eqn{(n-1)}- dimensional row vector, where each entry \eqn{i}
#' corresponds to the the number of bracnhes giving rise to \eqn{i} 
#' descendents. Hence, state 1 is always a vector of the form \eqn{(n,0,0,...,0)},
#' and state 2 is always given by the vector \eqn{(n-2,1,0,...,0)}.
#' 
#' @param n the sample size
#' 
#' @return The function returns a list containing the subintensity 
#' rate matrix \code{Rate.mat} and the state space matrix \code{StateSpace.mat}. 
#' In the latter, each row corresponds to a state and each state is a
#' \eqn{(n-1)}-dimensional row vector. 
#'
#' @examples
#' a <- BlockCountProcess(4)
#' a$Rate.mat
#' a$StateSpace.mat
#'
#' @export
BlockCountProcess <- function(n){
  ##----------------------------------------------------
  ## Possible states
  ##----------------------------------------------------
  ## Size of the state space (number of states)
  d <- P(n)
  ## Set of partitions of [n]
  Partition.mat <- parts(n)
  ## Rewriting the partitions as (a1,...,an)/
  ## Definition of the state space matrix
  StateSpace.mat <- t(apply(Partition.mat, 2,tabulate, nbins = n))
  
  ## Reordering
  StateSpace.mat <- StateSpace.mat[order(rowSums(StateSpace.mat),decreasing=TRUE),]
  ## Because of this ordering we can't 'go back', i.e.
  ## below the diagonal the entries are always zero
  
  ##----------------------------------------------------
  ## Intensity matrix
  ##----------------------------------------------------
  RateM <- matrix(0, ncol=d, nrow=d)
  ## Algorithm for finding rates between states
  for (i in 1:(d-1)){
    for (j in (i+1):d){
      # cat(i," state i",StSpM[i,])
      # cat(" ",j," state j",StSpM[j,])
      cvec <- StateSpace.mat[i,]-StateSpace.mat[j,]
      # cat(" cvec",cvec)
      ## Two branches are merged, i.e. removed from state i 
      check1 <- sum(cvec[cvec>0])==2
      # cat(" check1",check1)
      ## One new branch is created, i.e. added in state from j
      check2 <- sum(cvec[cvec<0])==-1
      # cat(" check2",check2)
      if (check1 & check2){
        ## Size(s) of the block(s) and the corresponding rates
        tmp <- StateSpace.mat[i,which(cvec>0)]
        RateM[i,j] <- ifelse(length(tmp)==1,tmp*(tmp-1)/2,prod(tmp))
      }
      
    }
  }
  ## Diagonal part of the rate matrix
  for (i in 1:d){
    RateM[i,i] <- -sum(RateM[i,])
  }
  return(list(Rate.mat = RateM[-nrow(RateM), -ncol(RateM)],
              StateSpace.mat = StateSpace.mat[-nrow(StateSpace.mat), 
                                              -ncol(StateSpace.mat)]))
}