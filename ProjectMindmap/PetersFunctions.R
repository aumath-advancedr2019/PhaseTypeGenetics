# The distribution function for PH(pivec, Tm) evaluated at q
pphasetype <- function(q, Tm, pivec) 1 - sum(pivec %*% expm(q*Tm))

# The probability function for DPH(initprob, subtrans) at k
ddiscphastyp <- function(initprob, subtrans, k){
  tmpnr <- length(initprob)
  tmpmat <- diag(1,tmpnr)-subtrans
  sum(initprob %*% (subtrans %^% (k-1)) %*% tmpmat) + (1 - sum(initprob))*(k==1)
}

# The density function for PH(pivec, Tm) evaluated at x
dcphasetype <- function(x, Tm, pivec) -sum(pivec %*% expm(x*Tm) %*% Tm)

# The p-quantile of PH(pivec, Tm)
#(Note that that the input '20' in uniroot is based on the plots of the distribution function 
#for Ttotal and Tmrca if other T matrices this should perhaps be another value)
qphasetype <- function(p, Tm, pivec){
  uniroot(function(x) pphasetype(x, Tm, pivec)-p, c(0, 20))$root[1]
}

# Reward transformation
rewardtransformparm <- function(rewards, initprob, subintensemat){
  d <- sum(rewards > 0)
  p <- length(initprob)
  if(p==1){
    list("newinitprob" = initprob, "newsubintensitymatrix" = as.matrix(subintensemat/rewards), "defect" = 0)
  }
  else{
    qmat <- matrix(rep(0,p^2), ncol = p)
    for(i in 1:p){
      
      for(j in (1:p)[-i]){
        qmat[i,j] <- -subintensemat[i,j]/subintensemat[i,i]
      }
    }
    ##
    # If all rewards are stricly postive everything is simpler
    ##
    if(d == p){
      pmat <- qmat
      alphavec <- initprob
    }
    else{
      qplusplus <- qmat[which(rewards > 0),which(rewards > 0)]
      qpluszero <- qmat[which(rewards > 0),which(rewards == 0)]
      qzeroplus <- qmat[which(rewards == 0),which(rewards > 0)]
      qzerozero <- qmat[which(rewards == 0 ),which(rewards == 0)]
      pmat <- qplusplus + qpluszero %*% solve(diag(1, nrow = p-d)-qzerozero) %*% qzeroplus
      piplus <- initprob[which(rewards > 0)]
      pizero <- initprob[which(rewards == 0)]
      alphavec <- piplus + pizero %*% solve(diag(1, nrow = p-d)-qzerozero) %*% qzeroplus 
      subintensemat <- as.matrix(subintensemat[which(rewards > 0), which(rewards >0)])
      rewards <- rewards[which(rewards > 0)]
    }
    pvec <- 1 - rowSums(pmat)
    Tstarmat <- matrix(rep(0,d^2), ncol = d)
    tstarvec <- rep(0,d)
    for(i in 1:d){
      for(j in (1:d)[-i]){
        Tstarmat[i,j] <- -subintensemat[i,i]/rewards[i]*pmat[i,j]
      }
      tstarvec[i] <- -subintensemat[i,i]/rewards[i]*pvec[i]
      Tstarmat[i,i] <- -sum(Tstarmat[i,])-tstarvec[i]
    }
    list("newinitprob" = alphavec, "newsubintensitymatrix" = Tstarmat, "defect" = 1 - sum(alphavec))
  }
}