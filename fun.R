J <- 2
L <- 4
l <- 1/L
pl <- c(rep(0,L)+l*c(1:L))


logistic.lik <- function(theta,Z,y,rou){   #经验似然函数
  alpha1 <- theta[1]
  alpha2 <- theta[2]
  beta <- theta[3]
  logl <- (sum(alpha1+Z[[1]]*beta)+sum(alpha2+Z[[2]]*beta)
           -sum(log(1+rou[1]*exp(alpha1+y[[1]]*beta)))-sum(log(1+rou[2]*exp(alpha2+y[[2]]*beta))))
  return(-logl)
}


#############ks统计量部分

Aky <- function(k,u,n,n0,alpha,beta,t){    #A函数
  rou <- (n-n0)/n0
  t.sort <- sort(t)
  u <- c(rep(u,n))
  expp <- exp(alpha+beta*t.sort)
  dd <- as.numeric(t.sort<=u)
  Aky <- sum(expp*t.sort^k*dd/(1+rou*expp)^2)/n0
  return(Aky)
}




#########卡方统计量部分
KK <- function (x,y,N,n,n0,alpha,beta,C,S){#卡方统计量
  rou <- (n-n0)/n0
  xi <- (n-n0)/N
  maxt <- max(y)
  y.sort <- sort(y)
  rou <- (n-n0)/n0
  dsmell <- quantile(y,  probs =pl , type=1)  
  
  
  d <- c(min(y)-1,dsmell)
  B <- matrix(nrow=2,ncol=L-1)
  D <- c(1:(L-1)) 
  for (l in 1:(L-1)){     
    B[,l] <- c(Aky(0,d[l+1],n,n0,alpha,beta,y)
               -Aky(0,d[l],n,n0,alpha,beta,y),
               Aky(1,d[l+1],n,n0,alpha,beta,y)
               -Aky(1,d[l],n,n0,alpha,beta,y))  
    D[l] <- Aky(0,d[l+1],n,n0,alpha,beta,y)-Aky(0,d[l],n,n0,alpha,beta,y)
  }
  
  D <- diag(D)  
  ####Qtilde计算
  Dsmell <- matrix(c(rep(as.numeric(dsmell),each=n)),nrow=n)
  p1 <- 1/((1+rou*exp(alpha+beta*y.sort))*n0)
  tt <- matrix(c(rep(y.sort,L)),nrow=n)
  p1 <- matrix(c(rep(p1,L)),nrow=n)
  dd <- matrix(as.numeric(tt<=Dsmell),nrow=n)
  q0 <- apply(p1*dd,2,sum)
  q00 <- c(0,q0[1:L-1])
  q_tilde <- q0-q00
  ######以下是qhead的计算方法
  #x <- y[1:n0]
  Fn <- ecdf(x)
  qhead0 <- Fn(dsmell)
  qhead00 <-c(0,qhead0[1:L-1])
  q_head <- qhead0-qhead00 
  
  Q_tilde <- sqrt(N)*(q_head-q_tilde)[1:L-1]  
  
  RJ <- rou^2*(D/xi-t(B)%*%C%*%solve(S)%*%t(C)%*%B)
  
  K <- tryCatch({
    t(Q_tilde)%*%solve(RJ)%*%Q_tilde
  },
                error = function(err){
                  return(NA)
                }
  )
  
  return(K)
}


dealt <- function(x,y,N,n,n0,alpha,beta){   #dealt统计量
  rou <- (n-n0)/n0
  y.sort <- sort(y) 
  yy <- as.numeric(!duplicated(y.sort))
  repeats <- tapply(y.sort, y.sort, length)
  repeats <- rep(repeats, repeats)
  yyy <- yy*repeats
  p1 <- 1/((1+rou*exp(alpha+beta*y.sort))*n0)
  p <- p1*yyy
  G <- as.numeric(cumsum(p))
  
 # x <- y[1:n0]
  Fn <- ecdf(x)
  Ghead <- Fn(y.sort)
  return(max(sqrt(N)*abs(Ghead-G)))
}

#########kernel部分

gauss <- function(x) 1/sqrt(4*pi) * exp(-(x^2)/4) 
kernel <- function(t,N,n,n0,alpha,beta){
  n1 <- n-n0
  rou <- n1/n0
  y <- c(rep(0,n0),rep(1,n1))
  r <- (y - rou*exp(alpha+beta*t)/(1+rou*exp(alpha+beta*t)))
  x1 <- matrix(rep(t,n),nrow=n)
  x2 <- t(x1)
  Q <- N/(n0*n0) * gauss(x1-x2)
  I <- t(r)%*%Q%*%r
}
