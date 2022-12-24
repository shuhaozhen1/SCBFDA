#' This part considers multi-dimensional functional data generation.
#'
#'
#'
#'


library(L1pack)
### Use Fourier basis functions to generate random process

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


p <- 3
transfermatrix <- matrix(0,p,p)
for(i in 1:p){
  for(j in 1:p){
    transfermatrix[i,j] <- (p-abs(i-j))/p
  }
}

mi <- 110

re <- 0.01

tpoints <- sort(runif(mi))

rawx <- matrix(replicate(n=p, Laplace.process(tpoints)), nrow = mi)


y0 <- transfermatrix %*% t(rawx)

y <- y0 + rnorm(mi,0,re)



