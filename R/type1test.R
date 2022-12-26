library(foreach)
library(doParallel)

####
n <- 40
m <- 5
###
mean1 <- function(x){
  sin(pi*x)
}

mean2 <- function(x){
  2*sin(pi*x)
}

mean3 <- function(x){
  -sin(pi*x)
}

meanlist <- list(mean1,mean2,mean3)

####

####
### Transfer matrix
p <- 3
transfermatrix <- matrix(0,p,p)
for(i in 1:p){
  for(j in 1:p){
    transfermatrix[i,j] <- (p-abs(i-j))/p
  }
}
####
ttoest <- seq(0,1,length.out=21)




registerDoParallel(detectCores()-1)
foreach (rp=1:10, .combine=cbind, .packages = 'L1pack',.errorhandling = 'remove') %dopar% {
  ty <- dataforn(n,m)
  tdata <- do.call(rbind, ty)

  # hcv <- bw.bcv(tdata[,1])
  hcv <- bw.ucv(tdata[,1])

  muhat <- localp(tdata, t=ttoest,h = hcv)

  centerty <- lapply(ty, center, totaldata =tdata, h=hcv)

  centerdata <-  do.call(rbind, centerty)


  bstime <- 1000
  bsresult <- replicate(bstime,gausshat(centerty,h=hcv),simplify = F)
  bsresult2 <- lapply(bsresult,function(x){x^2})

  sdest <- sqrt(Reduce('+',bsresult2)/bstime)

  sdbsresult <- lapply(bsresult,function(x,sd){return(x/sd)}, sd=sdest)

  qhat <- quantile( sapply(sdbsresult,function(x){max(abs(x))}) ,probs = 0.95)

  scbup <- muhat + qhat * sdest / sqrt(n)
  scblow <- muhat - qhat * sdest / sqrt(n)

  type1 <- 1 - floor(mean((scbup >= t(sapply(meanlist,function(x){x(ttoest)})) ) & (t(sapply(meanlist,function(x){x(ttoest)})) >= scblow )))

  c(type1, mean(scbup-scblow))
}

stopImplicitCluster()

