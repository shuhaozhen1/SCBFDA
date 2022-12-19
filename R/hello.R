data <- irreg.fd(mu=1,X=kl.process(distribution = 'LAPLACE'),n=100,m=10)

t_est <- seq(0,1,length.out=21)

EpaK <- function(x){
  y <- ifelse (x^2 > 1, 0 ,  3/4*(1-x^2))
  return(y)
}


locallineart <- function(data, t, h= 0.1 ) {
  S0 <- 0
  S1 <- 0
  S2 <- 0
  R0 <- 0
  R1 <- 0

  n <- length(data$t)
  for (i in 1: n) {
    mi <- length(data$t[[i]])
    for (j in 1:mi) {
      S0 <- S0 + mi^{-1}* EpaK((data$t[[i]][j]-t)/h)/h
      S1 <- S1 + mi^{-1}* EpaK((data$t[[i]][j]-t)/h)/h * ((data$t[[i]][j]-t)/h)
      S2 <- S2 + mi^{-1}* EpaK((data$t[[i]][j]-t)/h)/h * ((data$t[[i]][j]-t)/h)^2

      R0 <- R0 + mi^{-1}* EpaK((data$t[[i]][j]-t)/h)/h * data$y[[i]][j]
      R1 <- R1 + mi^{-1}* EpaK((data$t[[i]][j]-t)/h)/h * data$y[[i]][j] * ((data$t[[i]][j]-t)/h)
    }
  }

  hatt <- (R0*S2-R1*S1)/(S0*S2-S1^2)
  return(hatt)
}

locallinear <- function(t_est, data, h) {
  hat <- sapply(t_est, locallineart, data = data, h = h)
  return(hat)
}


h_c <- 10^seq(-2,0,length.out=21)

vc <- function(data, h) {
  cvsqr <- 0
  for(i in 1:length(data$t)){
    newdata <- list(t = data$t[-i],y = data$y[-i])
    oldt <- data$t[[i]]
    oldy <- data$y[[i]]
    cvest <- locallinear(oldt, data = newdata, h= h)

    cvsqr <- cvsqr + mean((oldy-cvest)^2)
  }
  return(cvsqr/length(data$t))
}

cvk <- function(data, h, k = 5) {
  n <- length(data$t)
  nk <- floor(n/k)
  cvsqr <- 0
  for(i in 1:nk){
    out <- (k*(i-1)+1):k*i
    newdata <- list(t = data$t[-out],y = data$y[-out])
    oldt <- unlist(data$t[out])
    oldy <- unlist(data$y[out])
    cvest <- locallinear(oldt, data = newdata, h= h)

    cvsqr <- cvsqr + mean((oldy-cvest)^2)
  }
  return(cvsqr/length(data$t))
}

hcvk <- h_c[which.min(sapply(h_c, cvk, data=data))]

hat <- locallinear(t_est, data= data, h=hcvk)

plot(hat)


