### Kernel functions
EpaK <- function(x){
  y <- ifelse (x^2 > 1, 0 ,  3/4*(1-x^2))
  return(y)
}
###
tpoints <- seq(0,1,length.out=21)

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
datafori <- function(m,p = 5,sigma=0.1, mean,distribution='gaussian'){
  mcandidate <- (m-2):(m+2)
  mi <- sample(mcandidate,1, replace = T)
  tpoints <- sort(runif(mi))
  rawx <- matrix(replicate(n=p, rpatx(tpoints,distribution = distribution )), nrow = mi)

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

dataforn <- function(n,m,p=5,mean,sigma=0.1,distribution='gaussian') {
  tdata <- replicate(n,datafori(m=m,p=p,mean=mean,sigma=sigma,distribution = distribution))
  return(tdata)
}

###

a <- dataforn(n=100,m=5,mean=meanlist,distribution = 'exponential')

