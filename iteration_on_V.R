iteration <- function(NuisanceFit_folds, TargetFit_folds) {
  
  tol        <- 0.0001 # adjust later
  max_iter   <- 100
  prev_avg_V <- Inf
  gamma      <- gamma0 
  step_size <- 1 # set at = 0.1/0.2 for gradient update
  
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
    
    if (abs(avg_V - prev_avg_V) < tol) {
      #message("Converged!")
      break
    }
    
    ## update for next iteration
    prev_avg_V <- avg_V
    gamma      <- gamma_new
  }
  
  for (j in 1:K){
    tau_numr[j] <- Tau(NuisanceFit_folds[[j]], TargetFit_folds[[j]], beta_opt)$numr
    tau_denom[j] <- Tau(NuisanceFit_folds[[j]], TargetFit_folds[[j]], beta_opt)$denom
  }
  
  tau <- 1-mean(tau_numr)/mean(tau_denom)
  tau
}