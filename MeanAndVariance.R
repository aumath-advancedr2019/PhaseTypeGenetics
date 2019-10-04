
mean.cphasetype <- function(initDist, Tmat){
  
  return(LaplacePhaseType(initDist, Tmat, i=1))
  
}

var <- function(...){
  
  UseMethod("var")
}

var.default <- function(x,...){
  
  var(x,...)
}

var.cphasetype <- function(){
  
  return(LaplacePhaseType(initDist = q, Tmat = TmatMRCA, i=2)-
    LaplacePhaseType(initDist = q, Tmat = TmatMRCA, i=1)^2)
}

