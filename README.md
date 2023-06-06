# ejMCMC

Summary
-------
ejMCMC is an R package of a recently proposed approximate Bayesian computing method that targets for problems with intractable or missing likelihood function. It is early rejection Markov chain Monte Carlo (ejMCMC) sampler based on Gaussian processes to accelerate inference speed. 

Installation
------------
install_github(``caofff/ejMCMC'')

Quick Start
------------


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
-----
- num_pre: The number of predictions. 
- num_simu: The number of simulations.
- num_acc: The number of accepted parameters.
- Theta_re: The ABC posterior distribution of parameters.

Demo
-----
We refer users to 'R/Example.R' for  an example. 


