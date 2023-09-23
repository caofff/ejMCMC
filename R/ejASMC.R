

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
K_fun <- function(x,h,type='Uniform'){
  Class <- c('Uniform','Tringular','Epanechnikov','Quartic','Triweight','Tricube')
  re <- rep(0,length(x))
  switch (which(Class==type),
          #Uniform
          re[abs(x/h)<1] <- 1,
          #Tringular
          re[abs(x/h)<1] <- 1-abs(x[abs(x/h)<1]/h),
          #Epanechnikov
          re[abs(x/h)<1] <- 1-(x[abs(x/h)<1]/h)^2,
          #Quartic
          re[abs(x/h)<1] <- (1-(x[abs(x/h)<1]/h)^2)^2,
          #Triweight
          re[abs(x/h)<1] <- (1-(x[abs(x/h)<1]/h)^2)^3,
          #Tricube
          re[abs(x/h)<1] <- (1-abs(x[abs(x/h)<1]/h)^3)^3
  )
  return(re)
}
#' Early rejection adaptive sequential Monte Carlo
#'
#' @param N Number of parameters at each iteration
#' @param eps.tag  Target threshold epsilon
#' @param theta0 Initial parameters
#' @param dis0 Initial discrepancies
#' @param Dis Discrepancy function, the input of which is parameter theta
#' @param dprior The prior density
#' @param h_fun The prediction function  h
#' @param Kernel The type of kernel function in the definition of ABC, Kernel=Uniform','Tringular','Epanechnikov','Quartic','Triweight' or 'Tricube'.
#' @param gamma1 The scale parameter to update threshold
#' @param gamma2 gamma2*N is the effective sample size
#'
#' @return the parameters from the posterior distribution
#' @export
#'
#' @examples
ejASMC <- function(N,eps.tag,theta0,dis0,Dis,dprior,h_fun,Kernel='Uniform',gamma1=0.8,gamma2=0.5){
  t=1
  Theta1 <- theta0
  dis1 <- dis0
  dis_hat_old <- dis1
  Weight <- data.frame(matrix(NA,nrow=N,ncol=1))
  Weight[,t] <- 1/N
  Theta_Re <- list(Theta1)
  Dis_Re <- as.data.frame(dis1)
  Eps <- c(10^6)
  while(Eps[t]>eps.tag){
    t=t+1
    U_Dis <- unique(Dis_Re[Weight[,t-1]>0,t-1])
    Eps[t] <- quantile(U_Dis,gamma1)
    if(Eps[t]<eps.tag){
      Eps[t]==eps.tag
    }
    id.alive <- which(Weight[,t-1]>0)
    Weight[,t] <- 0
    Weight[id.alive,t] <- Weight[id.alive,t-1]*K_fun(Dis_Re[id.alive,t-1],Eps[t],type=Kernel)/
      K_fun(Dis_Re[id.alive,t-1],Eps[t-1],type=Kernel)
    Weight[,t] <- Weight[,t]/sum(Weight[,t])
    if(1/sum(Weight[,t]^2)<gamma2*N+1){
      id <- Resample(Weight[,t],N)
      Theta_Re[[t-1]] <- Theta_Re[[t-1]][id,]
      Dis_Re[,t-1] <- Dis_Re[id,t-1]
      dis_hat_old <- dis_hat_old[id]
      Weight[,t] <- 1/N
    }
    pro_var <-  var(log(Theta_Re[[t-1]][Weight[,t]>0,]))
    Theta_Re[[t]] <- Theta_Re[[t-1]]
    Dis_Re[,t] <- Dis_Re[,t-1]
    #proposal(w>0)
    id.alive <- which(Weight[,t]>0)
    Theta_temp <- Theta_Re[[t-1]]
    Theta_temp[id.alive,] <- exp(t(apply(log(Theta_Re[[t-1]][id.alive,]),1,function(x) mvrnorm(1,x,pro_var))))
    w <- runif(N)
    #ordinary ejMCMC
    alpha_breve <- apply(Theta_temp,1,dprior)/(apply(Theta_Re[[t-1]],1,dprior)*
                                                 K_fun(Dis_Re[,t-1],Eps[t],type=Kernel))
    id1 <- which(w<alpha_breve)
    #ejMCMC
    dis_hat <- dis_hat_old
    dis_hat[id1] <- h_fun(Theta_temp[id1,])
    alpha_tilde <- alpha_breve*K_fun(dis_hat,Eps[t],type=Kernel)
    id2 <- which(w<alpha_tilde)
    # Simulation
    dis_temp <- Dis_Re[,t-1]
    dis_temp[id2] <- apply(Theta_temp[id2,],1,Dis)
    alpha_hat <- alpha_breve*K_fun(dis_temp,Eps[t],type=Kernel)
    id3 <- which(w<alpha_hat)
    Theta_Re[[t]][id3,] <- Theta_temp[id3,]
    Dis_Re[id3,t] <- dis_temp[id3]
    dis_hat_old[id3] <- dis_hat[id3]
  }
  id <- Resample(Weight[,t],N)
  Theta_re <- Theta_Re[[t]][id,]
  return(Theta_re)
}
