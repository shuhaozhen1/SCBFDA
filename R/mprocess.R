#' This part considers multi-dimensional functional data generation.
#'
#'
#'
#'

library(L1pack)
###
### Kernel functions
EpaK <- function(x){
  y <- ifelse (x^2 > 1, 0 ,  3/4*(1-x^2))
  return(y)
}
###
#### ### Use Fourier basis functions to generate random process
Laplace.process <- function(x , tvar = 1 , k=100){
  num <- 1:k
  coefx <- as.vector(rmLaplace(n=1,Scatter = diag(sqrt(tvar/(2^(2: (k+1)) ) ) ) ) )
  basis <- sapply(num,function(n,x){
    (1- (-1)^n)/2 * cos(2* pi * n * x) + (1- (-1)^n)/2 * sin(2* pi * n * x)
  },x=x)
  if(length(x) == 1) {
    valueatx <- sum( basis * coefx )
  } else {
    valueatx <- colSums(t(basis) * coefx)
  }
  return(valueatx)
}
####

### Transfer matrix
p <- 3
transfermatrix <- matrix(0,p,p)
for(i in 1:p){
  for(j in 1:p){
    transfermatrix[i,j] <- (p-abs(i-j))/p
  }
}


### data generation for one trajectory
datafori <- function(m,sigma=0.1){
  mcandidate <- (m-2):(m+2)
  mi <- sample(mcandidate,1, replace = T)
  tpoints <- sort(runif(mi))
  rawx <- matrix(replicate(n=p, Laplace.process(tpoints)), nrow = mi)

  y0 <- transfermatrix %*% t(rawx)

  y <- y0 + rnorm(mi,0,sigma)
  return(cbind(tpoints,t(y),rep.int(mi,mi)))
}

dataforn <- function(n,m,sigma=0.1) {
  tdata <- replicate(n,datafori(m,sigma))
  return(tdata)
}

### data

ty <- dataforn(100,3)
tdata <- do.call(rbind, ty)


### localpoly

localpt <- function(data,t,d=1,h){
  Tnmt <- sapply(0:d, function(data,t,d){
    (data[,1]-t)^d
  }, data=tdata,t=t)
  Wnmt <- diag(EpaK((tdata[,1]-t)/h) / h / tdata[,p+2])
  Ynmt <- tdata[,2:(p+1)]

  estd <- (solve(t(Tnmt) %*% Wnmt %*% Tnmt)) %*% (t(Tnmt) %*% Wnmt %*% Ynmt)
  muhat <- estd[1,]
  return(muhat)
}



localp <- function(data, t, d=1, h){
  hat <- sapply(t,localpt,data=data,d=d, h=h)
  return(hat)
}

localp(t=seq(0,1,length.out=11), data=tdata, h=0.1)


