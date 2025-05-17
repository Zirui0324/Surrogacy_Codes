library(tidyverse)
library(rje)
library(caret)
library(ks)
library(kernlab)

N=1000
n=1000
K=5
givenestimator = 'glm' # c('glm', 'svmLinear2', 'rf', 'nnet)

sim_function <- function(N, n, K, givenestimator, NBoot) {
  
  # haven't run this function for simulation, could just use the codes in: for (i in 1:N) {......}
  # specify outputs here
  # specify outputs here
  # specify outputs here
  # specify outputs here

  for (i in 1:N) {
    
    DP <- GD(n) 
    Index <- Split2(K, n) # sample split
    dimS <- ncol(DP$S0) # dim(theta) = dim(S)
    dimX <- 3
    gamma0 <- matrix(1e-6, dimX+1, 1) # gamma = (000) would case error
    
    # array for storing B, C, D from the K folds
    NuisanceFit_folds <- vector("list", K)
    TargetFit_folds   <- vector("list", K)
    P_hat_folds <- array(NA, dim = c(dimS, dimS, K))
    Q_hat_folds <- matrix(NA, nrow=dimS, ncol=K)
    V <- numeric(K)
    gamma_i <- matrix(NA, nrow=dimS, ncol=K)
    ha_folds <- array(NA, dim = c(dimX+1, dimX+1, K))
    xi_folds <- matrix(NA, nrow=dimX+1, ncol=K)
    #B_boo_folds <- array(NA, dim = c(dimS, NBoot, K))
    #C_boo_folds <- array(NA, dim = c(dimS, dimS, NBoot, K))
    #D_boo_folds <- array(NA, dim = c(NBoot, K))
    
    for (j in 1:K){
      
      idNj = (Index!=j) # for training, (1-1/K)*n dp
      idj = (Index==j) # for estimation, (1/K)*n dp
      
      # training set
      In <- lapply(DP, function(component) {
        if (is.matrix(component)) {
          component[idNj, , drop = FALSE]
        } else {
          component[idNj]
        }})
      
      # estimation set
      Out <- lapply(DP, function(component) {
        if (is.matrix(component)) {
          component[idj, , drop = FALSE]
        } else {
          component[idj]
        }})
      
      # target set
      Target <- DP$Xt
      
      # estimation
      result_j <- tryCatch(
        estimate(In, Out, Target, givenestimator, NBoot),
        error = function(e) {
          list(
            NuisanceFit = list(),
            TargetFit = list()
            # se_theta= rep(NaN, dimS),
            # se_tau  = NaN
          )
        }
      )
      
      # point estimate:
      NuisanceFit_folds[[j]]  <- result_j$NuisanceFit
      TargetFit_folds[[j]]    <- result_j$TargetFit
      
    }
    
    ###----------------------------------------- Iteration -----------------------------------------###
    
    tol        <- 0.1 # adjust later
    max_iter   <- 100
    prev_avg_V <- Inf
    gamma      <- gamma0 
    
    for (iter in 1:max_iter) {
      
      # P, Q:
      for (j in 1:K) {
        P_hat_folds[,,j]  <- calculation_pq(NuisanceFit_folds[[j]], TargetFit_folds[[j]], gamma)$P # can also move out
        Q_hat_folds[,j]   <- calculation_pq(NuisanceFit_folds[[j]], TargetFit_folds[[j]], gamma)$Q
      }
      
      # Beta: (average P, Q over folds)
      avg_P    <- apply(P_hat_folds, c(1, 2), mean, na.rm = TRUE)
      avg_Q    <- rowMeans(Q_hat_folds, na.rm = TRUE)
      beta_opt <- solve(avg_P) %*% avg_Q
      
      # V, ha, xi:
      for (j in 1:K){
        V[j] <- Vab(NuisanceFit_folds[[j]], TargetFit_folds[[j]])
        ha_folds[,,j] <- new_gamma(NuisanceFit_folds[[j]], TargetFit_folds[[j]], gamma, beta_opt)$ha
        xi_folds[,j]  <- new_gamma(NuisanceFit_folds[[j]], TargetFit_folds[[j]], gamma, beta_opt)$xi
      }
      
      avg_V <- mean(V) # evaluate the convergence of avg_V
      avg_ha   <- apply(ha_folds, c(1, 2), mean, na.rm = TRUE)
      avg_xi   <- rowMeans(xi_folds, na.rm = TRUE)
      gamma_new <- gamma - solve(avg_ha) %*% avg_xi
      
      cat(sprintf("iter %d : avg_V = %.5f  |ΔV| = %.5f\n",
                  iter, avg_V, abs(avg_V - prev_avg_V)))
      
      
      if (abs(avg_V - prev_avg_V) < tol) {
        message("Converged!")
        break
      }
      
      ## update for next iteration
      prev_avg_V <- avg_V
      gamma      <- gamma_new
    }

    
}}
