#'
#'
#'
#'
#'
#'
#'
#'
#'


library(foreach)
library(doParallel)
### Kernel functions
EpaK <- function(x){
  y <- ifelse (x^2 > 1, 0 ,  3/4*(1-x^2))
  return(y)
}
###
tpoints <- seq(0,1,length.out=31)

mean1 <- function(x){
  sin(pi*x)
}

mean2 <- function(x){
  2*sin(pi*x)
}

mean3 <- function(x){
  -sin(pi*x)
}

mean4 <- function(x){
  cos(pi*x)
}

mean5 <- function(x){
  -cos(pi*x)
}

meanlist <- list(mean1,mean2,mean3,mean4,mean5)









## KL-representation of random process
rpatx <- function(x , tvar = 1 , k=100, distribution='gaussian'){
  var <- tvar/(2^(1:k))
  if (distribution == 'gaussian') {
    coefx <- rnorm(k,mean=0, sd = sqrt(var) )
  } else if (distribution == 'exponential') {
    coefx <- rexp(k, rate = var^(-1/2))- var^(1/2)
  }
  basis <- sapply(1:k,function(n,x){
    (1- (-1)^(n))/2 * cos(2* pi * n * x) + (1- (-1)^(n+1))/2 * sin(2* pi * n * x)
  },x=x)
  if(length(x) == 1) {
    valueatx <- sum( basis * coefx )
  } else {
    valueatx <- colSums(t(basis) * coefx)
  }
  return(valueatx)
}

### p-dimensional data generation



## data generation for one trajectory
datafori <- function(m,p,sigma=0.1, mean, distribution='gaussian'){
  mcandidate <- (m-2):(m+2)
  mi <- sample(mcandidate,1, replace = T)
  tpoints <- sort(runif(mi))
  rawx <- matrix(replicate(n=p, rpatx(tpoints,distribution = distribution)), nrow = mi)

  ### Transfer matrix
  transfermatrix <- matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      transfermatrix[i,j] <- (p-abs(i-j))/p
    }
  }
  ###

  y0 <- transfermatrix %*% t(rawx)

  mu <- matrix(sapply(mean, function(x){x(tpoints)}), nrow=mi)
  y <- y0 + t(mu) + rnorm(mi,0,sigma)
  return(cbind(tpoints,t(y),rep.int(mi,mi)))
}

dataforn <- function(n,m,p=5,mean,sigma=0.1, distribution='gaussian') {
  tdata <- replicate(n,datafori(m=m,p=p,mean=mean,sigma=sigma,distribution=distribution))
  return(tdata)
}

### data simulation



### localpoly

localpt <- function(data,t,d=1,h){
  p <- ncol(data)-2
  Tnmt <- sapply(0:d, function(data,t,d){
    (data[,1]-t)^d
  }, data=data,t=t)
  Wnmt <- diag(EpaK((data[,1]-t)/h) / h / data[,p+2])
  Ynmt <- data[,2:(p+1)]

  estd <- (solve(t(Tnmt) %*% Wnmt %*% Tnmt)) %*% (t(Tnmt) %*% Wnmt %*% Ynmt)
  muhat <- estd[1,]
  return(muhat)
}



localp <- function(data, t, d=1, h){
  hat <- sapply(t,localpt,data=data,d=d, h=h)
  return(hat)
}

### get centered data
# center <- function(data,totaldata,h){
#   p <- ncol(data)-2
#   est <- localp(t=data[,1], data=totaldata, h=h)
#   data[,2:(p+1)] <- data[,2:(p+1)] - t(est)
#   return(data)
# }
#
#
# ### individual smoothing
# ismootht <- function(t,data,h){
#   p <- ncol(data)-2
#   Wmt <- EpaK((data[,1]-t)/h) / h / data[,p+2]
#   hatit <- matrix(Wmt * data[,2:(p+1)], ncol= p)
#   hatit <- colSums(hatit)
#   return(hatit)
# }
#
#
# ismooth <- function(tpoints,data,h){
#   totalsmooth <- sapply(tpoints,ismootht,data,h=h)
#   return(totalsmooth)
# }

# gaussbs <- function(data){
#   bsdata <- lapply(data,function(x){
#     rnorm(1)*x
#   })
#   bs <- Reduce('+',bsdata)/sqrt(n)
#   return(bs)
# }



gaussbs <- function(data,mis,p,t,h){
  rawdata <- data
  rawdata[,2:(p+1)] <-  rep(rnorm(length(mis)), mis) * data[,2:(p+1)]
  gshat <- localp(rawdata,t=t, h=h)
  return(sqrt(length(mis))*gshat)
}



###
truemean <- t(sapply(meanlist,function(x){x(tpoints)}))


n <- 100
m <- 5
sigma <- 0.1

registerDoParallel(detectCores())

results <- foreach (rp=1:50, .combine=cbind, .errorhandling = 'remove') %dopar% {
  ty <- dataforn(n,m,mean=meanlist,sigma = sigma)

  mis <- sapply(ty, nrow)

  matrixformty <- do.call(rbind,ty)
  p <- ncol(matrixformty)-2

  tdata <- matrixformty[,1]
  mdata <- matrixformty[,p+2]



  hcv1 <-1.2 * bw.bcv(matrixformty[,1])

  hcv <- hcv1/1.5

  muhat <- localp(matrixformty, t = tpoints, h=hcv1)

  hatij <- localp(matrixformty, t = matrixformty[,1], h=hcv1)

  centerdatay <- matrixformty[,2:(p+1)] - t(hatij)
  centerdata <- cbind(tdata,centerdatay,mdata)

  gs <- replicate(400,gaussbs(centerdata, mis=mis, p=p, t= tpoints, h=hcv),simplify = F)

  sdest <- sqrt(Reduce('+',lapply(gs, function(x){x^2}))/n)

  gssd <- lapply(gs, function(x,sd){x/sd},sd=sdest)

  supdata <- sapply(gssd, function(x){
    max(abs(x))
  })

  qhat <- quantile(supdata,0.95)

  bandup <- muhat + qhat * sdest/sqrt(n)
  bandlow <- muhat - qhat* sdest/sqrt(n)

  type1 <- 1 - floor(mean((bandup >= truemean & bandlow <= truemean)))

  bandwidth <- mean(bandup-bandlow)

  c(type1, bandwidth)
}


stopImplicitCluster()

rowMeans(results)
