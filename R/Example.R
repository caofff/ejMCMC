library(deSolve)
require(stats)
library(MASS)
library(hetGP)
DE_model <- function(t, state, parameters) {
  ## assign parameter values to parameter variables
  ## the function 'unname()' just removes the name of the parameter - this is unnecessary, it just cleans things up a bit
  kappa1 <- unname(parameters['kappa1'])
  kappa2 <- unname(parameters['kappa2'])

  ## assign state variables to names
  X1 <- unname(state['X1'])
  X2 <- unname(state['X2'])

  ## compute the rates of change
  dX1dt <- 72/(36+X2) - kappa1
  dX2dt <- kappa2*X1 - 1

  ## return as a list object, with the first element of the list being the derivatives. The order of derivatives must be the same as the order in the initial condition vector!
  return(list(c(dX1dt, dX2dt)))
}
times <- seq(0, 60, 0.5)
initial <- c(X1=7, X2=-10)
parameters <- c(kappa1 = 2, kappa2 = 1 )
## Simulate the system
y_fun <- function(theta,seed=0){
  names(theta) <- NULL
  output <- ode(y=initial, times=times, func=DE_model, parms=c(kappa1 = theta[1], kappa2 = theta[2] ))
  re <- matrix(0,nrow(output),2)
  set.seed(seed)
  re[,1] <- output[,'X1']+rnorm(nrow(output),0,1)
  re[,2] <- output[,'X2']+rnorm(nrow(output),0,3)
  return(c(re[,1],re[,2]))
}
Resample <- function(W,N){
  n.re <- rep(0,N)
  u <- (runif(1)+0:(N-1))/N
  Psum <- cumsum(W)
  i=1
  for(j in 1:N){
    while(i<N+1&Psum[j]>u[i]){
      i=i+1
      n.re[j] <- n.re[j]+1
    }
  }
  id <- rep(1:N,n.re)
  return(id)
}
y_obs <- y_fun(parameters,122)
Dis <- function(theta){
  theta <- as.vector(theta)
  names(theta) <- NULL
  pa<- c(kappa1 = theta[1], kappa2 = theta[2])
  output <- ode(y=initial, times=times, func=DE_model, parms=pa)
  y <- matrix(0,121,2)
  y[,1] <- output[,'X1']+rnorm(121,0,1)
  y[,2] <- output[,'X2']+rnorm(121,0,3)
  re <- mean(sqrt(apply((y-y_obs)^2,1,sum)))
  return(re)
}
theta1_lower <- 1.8;theta1_upper=2.2
theta2_lower <- 0.8;theta2_upper=1.2
rprior <- function(N){
  theta1 = runif(N,theta1_lower,theta1_upper)
  theta2 = runif(N,theta2_lower,theta2_upper)
  re <- cbind(theta1,theta2)
  return(re)
}
dprior <- function(theta){
  if(is.vector(theta)){
    re <- (theta[1]>theta1_lower)*(theta[2]>theta2_lower)*
      (theta[1]<theta1_upper)*(theta[2]<theta2_upper)
  }else if(dim(theta)[1]==1){
    re <- (theta[1]>theta1_lower)*(theta[2]>theta2_lower)*
      (theta[1]<theta1_upper)*(theta[2]<theta2_upper)
  }else{
    re <- (theta[1,]>theta1_lower)*(theta[2,]>theta2_lower)*
      (theta[1,]<theta1_upper)*(theta[2,]<theta2_upper)
  }
  return(re)
}

#Generte Train Data set
theta0 <- rprior(100)
dis0 <- matrix(apply(theta0,1,Dis),ncol=1)
Train_Data <- Train_SMC(100,3000,theta0,dis0,Dis,dprior)
#Trian the discrepancy model and obtain function h()
hom1 <- mleHomGP(as.matrix(Train_Data[,1:2]), Train_Data[,3], covtype = "Gaussian")
h_fun <- function(theta){
  if(is.vector(theta)){
    theta <- matrix(theta,nrow = 1)
  }else{
    theta <- as.matrix(theta)
  }
  predictions <- predict(x = theta, object =  hom1)
  re <- predictions$mean+qnorm(0.05)*sqrt(predictions$sd2+predictions$nugs)
  return(re)
}
#ejMCMC
Re2 <- ejMCMC(10,4.5,c(2,1),4.4,Dis,dprior,diag(c(0.05,0.05)),h_fun)
#ejASMC
theta1 <- rprior(3000)
dis1 <- matrix(apply(theta1,1,Dis),ncol=1)
Re2 <- ejASMC(3000,4,theta1,dis1,Dis,dprior,h_fun)
