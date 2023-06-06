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


#' @param N Number of parameters at each iteration.
#' @param N_sim  The number of simulations.
#' @param theta0 Initial parameters.
#' @param dis0 Initial discrepancies.
#' @param Dis Discrepancy function, the input of which is parameter theta.
#' @param dprior The prior density.
#' @param h_fun The prediction function  h.
#' @param Kernel The type of kernel function in the definition of ABC, Kernel=Uniform','Tringular','Epanechnikov','Quartic','Triweight' or 'Tricube'.
#' @param gamma1 The scale parameter to update threshold.
#' @param gamma2 gamma2*N is the effective sample size.
#' @param eps.tag  Target threshold epsilon.
#'
#' @return Train data set, the last column of which is discrepancies.
#' @export
#'
#' @examples
Train_SMC <- function(N,N_sim,theta0,dis0,Dis,dprior,Kernel='Uniform',gamma1=0.5,gamma2=0.5,eps.tag=0){
  t=1
  Theta1 <- theta0
  dis1 <- dis0
  Weight <- data.frame(matrix(NA,nrow=N,ncol=1))
  Weight[,t] <- 1/N
  Theta_Re <- list(Theta1)
  Dis_Re <- as.data.frame(dis1)
  Eps <- c(10^6)
  num_sim=N
  Theta_Train <- Theta1
  Dis_Train <- dis1
  while(num_sim<N_sim){
    t=t+1
    U_Dis <- unique(Dis_Re[Weight[,t-1]>0,t-1])
    Eps[t] <- sort(U_Dis)[min(gamma1*N,length(U_Dis))]
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
      Weight[,t] <- 1/N
    }
    pro_var <-  var(log(Theta_Re[[t-1]][Weight[,t]>0,]))
    Theta_Re[[t]] <- Theta_Re[[t-1]]
    Dis_Re[,t] <- Dis_Re[,t-1]
    #proposal(w>0)
    id.alive <- which(Weight[,t]>0)
    Theta_temp <- exp(t(apply(log(Theta_Re[[t-1]][id.alive,]),1,function(x) mvrnorm(1,x,pro_var))))
    w <- runif(length(id.alive))
    #ordinary ejMCMC
    alpha1 <- apply(Theta_temp,1,dprior)/
      (apply(Theta_Re[[t-1]][id.alive,],1,dprior)*K_fun(Dis_Re[id.alive,t-1],Eps[t],type=Kernel))
    id1 <- which(w<alpha1)
    # Simulation
    num_sim=num_sim+length(id1)
    dis_temp <- apply(Theta_temp[id1,],1,Dis)
    Theta_Train <- rbind(Theta_Train,Theta_temp[id1,])
    Dis_Train <- c(Dis_Train,dis_temp)
    alpha <- alpha1[id1]*K_fun(dis_temp,Eps[t],type=Kernel)
    id2 <- id1[w[id1]<alpha]
    Theta_Re[[t]][id.alive[id2],] <- Theta_temp[id2,]
    Dis_Re[id.alive[id2],t] <- dis_temp[w[id1]<alpha]
  }
  Train_Data_SMC <- cbind(Theta_Train,Dis_Train)
  Train_Data_SMC <- as.data.frame(Train_Data_SMC)
  colnames(Train_Data_SMC) <- c(paste('theta',1:ncol(Theta_Train)),'dis')
  return(Train_Data_SMC)
}

