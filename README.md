# ejMCMC

Summary
-------
ejMCMC is an R package of a recently proposed approximate Bayesian computing method that targets for problems with intractable or missing likelihood function. It is early rejection Markov chain Monte Carlo (ejMCMC) sampler based on Gaussian processes to accelerate inference speed. 

Installation
------------
install_github(``caofff/ejMCMC'')

Quick Start
------------

function Trian_SMC()
-----
A function to sequentially generate training data using the ABC SMC algorithm with OejMCMC as proposals.

Input Arguments
-----
- N: Number of parameters at each iteration.
- N_sim:  The number of simulations.
- theta0: Initial parameters.
- dis0: Initial discrepancies.
- Dis: Discrepancy function, the input of which is parameter theta.
- dprior: The prior density.
- h_fun: The prediction function  h.
- Kernel: The type of kernel function in the definition of ABC, Kernel=Uniform','Tringular','Epanechnikov','Quartic','Triweight' or 'Tricube'.
- gamma1: The scale parameter to update threshold.
- gamma2: gamma2*N is the effective sample size.
- eps.tag:  Target threshold epsilon.

Output
-----
- Train data set, the last column of which is discrepancies.


function ejMCMC()
-----
An early rejection ABC Markov chain Monte Carlo

Input Arguments
-----
- N: Total number of MCMC iterations.
- eps: ABC threshold.
- theta0: The initial parameter.
- dis0: The initial discrepancy.
- Dis: The discrepancy function, the input and output of which are parameters(theta) and discrepancy between the synthetic data and observed data.
- dprior: The prior density function.
- Sigma_prop: The covariance matrix of Gaussian proposal distribution.
- h_fun: The prediction function h().
- Kernel: The type of kernel function in the definition of ABC, Kernel=Uniform','Tringular','Epanechnikov','Quartic','Triweight' or 'Tricube'.

Output
- num_pre: The number of predictions. 
- num_simu: The number of simulations.
- num_acc: The number of accepted parameters.
- Theta_re: The ABC posterior distribution of parameters.



function ejASMC()
-----
An early rejection adaptive sequential Monte Carlo

Input Arguments
-----
- N: Number of parameters at each iteration
- eps.tag:  Target threshold epsilon
- theta0: Initial parameters
- dis0: Initial discrepancies
- Dis: Discrepancy function, the input of which is parameter theta
- dprior: The prior density
- h_fun: The prediction function  h
- Kernel: The type of kernel function in the definition of ABC, Kernel=Uniform','Tringular','Epanechnikov','Quartic','Triweight' or 'Tricube'.
- gamma1: The scale parameter to update threshold
- gamma2: gamma2*N is the effective sample size

Output
-----
- The parameters from the posterior distribution.

Demo
-----
We refer users to 'R/Example.R' for  an example.


