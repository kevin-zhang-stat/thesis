
#####################3
require(snowfall)
sfInit(parallel=TRUE, cpus=4, type="SOCK",socketHosts=rep("localhost",4))
source("fun.r",encoding="utf-8")
simk_power <- matrix(nrow=3,ncol=1000)
simdealt_power <- matrix(nrow=3,ncol=1000)
simI_power <- matrix(nrow=3,ncol=1000)
#######参数设计
C <- list()
C[[1]] <- matrix(c(1,0,0,0,0,1),nrow=2)
C[[2]] <- matrix(c(0,0,1,0,0,1),nrow=2)
exam2 <- matrix(nrow=25,ncol=3)
exam2[,1] <- c(0,0,1,0,0,1,1,0,1,0,0,0,0,0,1,1,0,1,0,1,0,0,1,1,1)
exam2[,2] <- c(0,1,1,0,0,0,1,0,0,0,0,1,1,0,1,0,0,1,0,1,0,1,0,1,1)
exam2[,3] <- c(4,4,3,2,7,12,8,1,9,5,6,4,2,5,9,19,7,10,5,8,14,8,14,9,17)
t_hat <- list()
t_hat[[1]] <- exam2[exam2[,2]==0,][,-2]
t_hat[[2]] <- exam2[exam2[,2]==1,][,-2]
x <- list()
z <- list()
t <- list()
for (j in 1:2){
  x[[j]] <- t_hat[[j]][t_hat[[j]][,1]==0,][,-1]
  z[[j]] <- t_hat[[j]][t_hat[[j]][,1]==1,][,-1]
  t[[j]] <- t_hat[[j]][,-1]
}




J <- 2
n <- rep(0,J)
n0 <-rep(0,J)
n1 <- rep(0,J)
rou <- rep(0,J)
for (j in 1:J){
  n[j] <- length(t[[j]])
  n0[j] <- length(x[[j]])
  n1[j] <- n[j]-n0[j]
  rou[j] <- n1[j]/n0[j]
}
N <- sum(n)
result <- optim(c(rep(1,J),1),logistic.lik,Z=z,y=t,rou=rou,method="BFGS")
alphahead <- result$par[1:J]
betahead <- result$par[J+1]
  
  
  ####抽样总体#
  ptx <- list()
  ptz <- list()
  t.unique <- list()
  for (j in 1:J){
    t.sort <- sort(t[[j]]) 
    tt <- as.numeric(!duplicated(t.sort))
    repeats <- tapply(t.sort, t.sort, length)
    repeats <- rep(repeats, repeats)
    ttt <- tt*repeats
    p1 <- 1/((1+rou[j]*exp(alphahead[j]+betahead*t.sort))*n0[j])
    p <- p1*ttt
    ###H分布
    ph <- exp(alphahead[j]+betahead*t.sort)/((1+rou[j]*exp(alphahead[j]+betahead*t.sort))*n0[j])
    phh <- ph*ttt
    t.unique[[j]] <- unique(t.sort)
    ptx[[j]] <- setdiff(as.numeric(p),0)
    ptz[[j]] <- setdiff(as.numeric(phh),0)
  }
  
  
  #######统计量#########  
  
  S <- matrix(nrow=J+1,ncol=J+1)
  xi <- (n-n0)/N
  S[1:J,1:J] <- diag(xi*mapply(Aky,0,Inf,n,n0,alphahead,betahead,t))
  S[1:J,J+1] <- xi*mapply(Aky,1,Inf,n,n0,alphahead,betahead,t)
  S[J+1,1:J] <- t(S[1:J,J+1])
  S[J+1,J+1] <- sum(xi*mapply(Aky,2,Inf,n,n0,alphahead,betahead,t))
  S_tilde <- list()
  S_tilde[[1]] <- S
  S_tilde[[2]] <- S
  
  K <- mapply(KK,x,t,N,n,n0,alphahead,betahead,C,S_tilde)
  K <- sum(n*K)/sum(n)
  
  Dealt <- mapply(dealt,x,t,N,n,n0,alphahead,betahead)
  Dealt <- sum(n*Dealt)/sum(n)
  
  
  I <- mapply(kernel,t,N,n,n0,alphahead,betahead)
  I <- sum(n*I)/sum(n)
  #######bootstrap模拟
  
  #######再抽样########
  sfExportAll()
  boot <- function(idx){
    x <- list()
    z <- list()
    t <- list()
    for (j in 1:J){
      x[[j]] <-sample(t.unique[[j]], size = n0[j], replace = TRUE, prob = ptx[[j]])
      z[[j]] <-sample(t.unique[[j]], size = n1[j], replace = TRUE, prob = ptz[[j]])
      t[[j]] <- c(x[[j]],z[[j]])
    }
    ####极大似然估计####
    result <- optim(c(rep(1,J),1),logistic.lik,Z=z,y=t,rou=rou,method="BFGS")
    alphahead <- result$par[1:J]
    betahead <- result$par[J+1]
    #######统计量#######
    
    S <- matrix(nrow=J+1,ncol=J+1)
    S[1:J,1:J] <- diag(xi*mapply(Aky,0,Inf,n,n0,alphahead,betahead,t))
    S[1:J,J+1] <- xi*mapply(Aky,1,Inf,n,n0,alphahead,betahead,t)
    S[J+1,1:J] <- t(S[1:J,J+1])
    S[J+1,J+1] <- sum(xi*mapply(Aky,2,Inf,n,n0,alphahead,betahead,t))
    S_tilde <- list()
    S_tilde[[1]] <- S
    S_tilde[[2]] <- S
    K_boot <- mapply(KK,x,t,N,n,n0,alphahead,betahead,C,S_tilde)
    K_boot <- sum(n*K_boot)/sum(n)
    
    Dealt_boot <- mapply(dealt,x,t,N,n,n0,alphahead,betahead)
    Dealt_boot <- sum(n*Dealt_boot)/sum(n)
    
    
    I_boot <- mapply(kernel,t,N,n,n0,alphahead,betahead)
    I_boot <- sum(n*I_boot)/sum(n)
    
    return(c(K_boot,Dealt_boot,I_boot))
  }
  
  Result <- sfSapply(1:1000,boot)
  sfStop ()
  simk <- Result[1,]
  simdealt <- Result[2,]
  simI <- Result[3,]
  simk <- simk[!is.na(simk)]
  pdealt <- mean(simdealt>=Dealt)
  pI <- mean(simI >= I)
  pk <- mean(simk >= K)





