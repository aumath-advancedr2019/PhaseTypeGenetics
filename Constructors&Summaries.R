# Constructor for continuous phase-type distribution
cphasetype <- function(initDist, T.mat){
  if((sum(initDist) > 1) | (sum(initDist < 0) > 0)) stop("Not a valid initial distribution")
  if(det(T.mat)==0) stop("Not a valid subintensity matrix")
  if(length(initDist) != nrow(T.mat)) stop("The dimensions of the input should be the same")
  tmp <- list("initDist" = initDist, "T.mat" = T.mat)
  class(tmp) <- "cphasetype"
  tmp
}

# Constructor for discrete phase-type distribution
dphasetype <- function(initDist, T.mat){
  if((sum(initDist) > 1) | (sum(initDist < 0) > 0)) stop("Not a valid initial distribution")
  if((det(T.mat) == 0) | (sum(T.mat < 0) == 0) | (sum(rowSums(T.mat) > 1) > 0)) {
    stop("Not a valid subtransition matrix") 
  }
  if(length(initDist) != nrow(T.mat)) stop("The dimensions of the input should be the same")
  tmp <- list("initDist" = initDist, "T.mat" = T.mat)
  class(tmp) <- "dphasetype"
  tmp
}

# Summary for continuous phase-type distribution
summary.cphasetype <- function(x){
  cat("A continuous phase-type distribution with initial probability vector \n")
  print(x$initDist)
  cat("and subintensity matrix \n")
  print(x$T.mat)
  if(sum(x$initDist) < 1){
    cat("and defect\n", 1-sum(x$initDist))
  }
}

# Summary for discrete phase-type distribution
summary.dphasetype <- function(x){
  cat("A discrete phase-type distribution with initial probability vector \n")
  print(x$initDist)
  cat("and subtransition matrix \n")
  print(x$T.mat)
  if(sum(x$initDist) < 1){
    cat("and defect\n", 1-sum(x$initDist))
  }
}

