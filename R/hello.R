### data without random error
data <- irreg.fd(mu=1,X=kl.process(distribution = 'LAPLACE'),n=100,m=5, sig=0)

### generate different Yp and add random error
random <- function(data) {
  yu <- rnorm(length(unlist(data)), 0, 0.1)
  return(yu)
}

unifadd <- function(data) {
  yu <- rnorm(2,2)*(unlist(data))
  return(yu)
}

gaussadd <- function(data) {
  ynorm<- rnorm(1,1)*(unlist(data))
  return(ynorm)
}

re <- lapply(data$t,random)
yu <- lapply(data$y,unifadd)
ynorm <- lapply(data$y,gaussadd)

### data is completed
tdata <- list()
for (i in 1:length(data$t)){
  tdata[[i]] <- cbind(data$t[[i]],as.numeric(data$y[[i]])+re[[i]],as.numeric(yu[[i]])+re[[i]],as.numeric(ynorm[[i]])+re[[i]])
}

rm(data)

### Kernel functions
EpaK <- function(x){
  y <- ifelse (x^2 > 1, 0 ,  3/4*(1-x^2))
  return(y)
}

### time points to be estimated
t_est <- seq(0,1,length.out=21)

### local linear estimator

SRi <- function(data, t , h=0.1){
  tij <- data[,1]
  s0ij <- (1/length(tij)) * EpaK((tij-t)/h)/h
  s1ij <- s0ij * ((tij-t)/h)
  s2ij <- s0ij * ((tij-t)/h)^2
  s0repeat <- replicate(3,s0ij)
  s1repeat <- replicate(3,s1ij)
  r0ij <- s0repeat * data[,-1]
  r1ij <- s1repeat * data[,-1]
  return(c(sum(s0ij),sum(s1ij),sum(s2ij),colSums(r0ij),colSums(r1ij)))
}

locallineart <- function(data, t, h = 0.1) {
  totalsr <- sapply(tdata, SRi, t = t , h= h)
  sr <- rowMeans(totalsr)

  s0 <- sr[1]
  s1 <- sr[2]
  s2 <- sr[3]
  r0p1 <- sr[4]
  r0p2 <- sr[5]
  r0p3 <- sr[6]
  r1p1 <- sr[7]
  r1p2 <- sr[8]
  r1p3 <- sr[9]

  hatp1 <- (r0p1 * s2 - r1p1 * s1)/(s0 * s2 - s1^2)
  hatp2 <- (r0p2 * s2 - r1p2 * s1)/(s0 * s2 - s1^2)
  hatp3 <- (r0p3 * s2 - r1p3 * s1)/(s0 * s2 - s1^2)
  hat <- c(hatp1,hatp2,hatp3)

  return(hat)
}

locallinear <- function(data , t_est , h=0.1) {
  hat <- sapply(t_est, locallineart, data = data, h=h)
  return(hat)
}

### bandwidth choice vc and vc k-fold
h_c <- 10^seq(-2,0,length.out=21)


cvk <- function(data, h, k = 5) {
  n <- length(data)
  nk <- floor(n/k)
  cvsqr <- 0
  for(i in 1:nk){
    out <- (k*(i-1)+1):k*i
    newdata <- data[-out]
    oldt <- unlist(sapply(data[out], function(x){x[,1]}))
    oldy <- unlist(sapply(data[out], function(x){x[,-1]}))
    cvest <- c(t(locallinear(oldt, data = newdata, h= h)))

    cvsqr <- cvsqr + mean((oldy-cvest)^2)
  }
  return(cvsqr/n)
}


hcvk <- h_c[which.min(sapply(h_c, cvk, data=tdata))]

hat <- locallinear(t_est, data= tdata, h=hcvk)

plot(hat[3,])



