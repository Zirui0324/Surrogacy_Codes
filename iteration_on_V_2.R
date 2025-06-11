
# incorporating gamma convergence in before each V
# V generally converges within 3 steps (tol = 1e-5)

iteration_V <- function(NuisanceFit_folds, TargetFit_folds) {
  
  # storage:
  P_hat_folds <- array(NA, dim = c(dimS, dimS, K))
  Q_hat_folds <- matrix(NA, nrow=dimS, ncol=K)
  V           <- numeric(K)
  tau_numr    <- numeric(K)
  tau_denom   <- numeric(K)
  
  # initiate:
  gamma <- gamma0
  prev_avg_V  <- Inf
  
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
    #print(beta_opt) # debug
    
    # Gamma iteration:
    gamma_new <- iteration_gamma(NuisanceFit_folds, TargetFit_folds, beta_opt, gamma)
    
    # V, ha, xi:
    for (j in 1:K){
      V[j] <- Vab(NuisanceFit_folds[[j]], TargetFit_folds[[j]], gamma_new, beta_opt)
    }

    avg_V <- mean(V) # then evaluate the convergence of avg_V
    delta_V <- abs(avg_V - prev_avg_V)
    
    # ░░░░░░░░░░░░░░░░░░░░░░░░░ display block: ΔV, γ, β ░░░░░░░░░░░░░░░░░░░░░░░░░░
    #cat(sprintf("iter %d : avg_V = %.5f  ΔV = %.5f\n",
    #            iter, avg_V, avg_V - prev_avg_V))
    #
    #cat(sprintf(
    #  "\niter %3d | |ΔV| = %.5f\n", iter, delta_V))
    #
    #cat("γ")
    #print(round(gamma_new, 6))                              # keeps 4×1 layout
    #
    #cat("β")
    #print(round(beta_opt, 6))                              # keeps 2×1 layout
    #cat("----------------------------------------------------------------\n")
    # ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
    
    if (abs(avg_V - prev_avg_V) < tol) {
      #message("Converged!")
      gamma <- gamma_new
      break
    }
    
    ## update for next iteration
    prev_avg_V <- avg_V
    
  }
  
  for (j in 1:K){
    tau_numr[j] <- Tau(NuisanceFit_folds[[j]], TargetFit_folds[[j]], beta_opt)$numr
    tau_denom[j] <- Tau(NuisanceFit_folds[[j]], TargetFit_folds[[j]], beta_opt)$denom
  }
  
  tau <- 1-mean(tau_numr)/mean(tau_denom)
  return(list(
    tau = tau,
    beta_opt = beta_opt,
    gamma = gamma
         ))
}
