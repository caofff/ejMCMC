library(MASS)
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

#' ejMCMC function
#'
#' @param N Total number of MCMC iterations.
#' @param eps ABC threshold.
#' @param theta0 the initial parameter.
#' @param dis0 the initial discrepancy.
#' @param Dis the discrepancy function, the input and output of which are parameters(theta) and discrepancy between the synthetic data and observed data.
#' @param dprior the prior density function.
#' @param Sigma_prop the covariance matrix of Gaussian proposal distribution.
#' @param h_fun the prediction function h().
#' @param Kernel The type of kernel function in the definition of ABC, Kernel=Uniform','Tringular','Epanechnikov','Quartic','Triweight' or 'Tricube'.
#'
#' @return the number of predictions, the number of simulations, the number of accepted parameters, the parameters from the ABC posterior distribution.
#' @export
#'
#' @examples
ejMCMC <- function(N,eps,theta0,dis0,Dis,dprior,Sigma_prop,h_fun,Kernel='Uniform'){
  num_pre=0
  p = length(theta0)
  Theta_re <- as.data.frame(matrix(0,N,p))
  colnames(Theta_re) <- paste0('theta',1:p)
  t=1;w_old=0;
  dis_old <- dis0
  dis_hat_old <- dis_old
  k_old <- min(K_fun(dis_old,eps,type=Kernel),K_fun(dis_hat_old,eps,type=Kernel))
  Theta_re[t,] <- theta0
  num_pre <- 0
  num_simu <- 0
  num_acc <- 0
  while(t<N){
    t=t+1
    theta_temp <- mvrnorm(1,t(Theta_re[t-1,]),Sigma_prop)
    alpha_breve <- dprior(theta_temp)/
      (dprior(Theta_re[t-1,])*k_old)
    w <- runif(1)
    if(w>alpha_breve){
      Theta_re[t,] <- Theta_re[t-1,]
    }else{
      dis_hat_temp <- h_fun(theta_temp)
      alpha_tilde <- alpha_breve*K_fun(dis_hat_temp,eps,type=Kernel)
      num_pre=num_pre+1
      if(w>alpha_tilde){
        Theta_re[t,] <- Theta_re[t-1,]
      }else{
        dis_temp <- Dis(theta_temp)
        num_simu <- num_simu+1
        k_temp <- min(K_fun(dis_temp,eps,type=Kernel),K_fun(dis_hat_temp,eps,type=Kernel))
        alpha_hat <- alpha_breve*k_temp
        if(w<alpha_hat){
          num_acc = num_acc+1
          Theta_re[t,] <- theta_temp
          dis_old=dis_temp
          k_old=k_temp
        }else{
          Theta_re[t,] <- Theta_re[t-1,]
        }
      }
    }
  }
  Re <- list(num_pre=num_pre,num_simu=num_simu,num_acc=num_acc,Theta_re=Theta_re)
  return(Re)
}
