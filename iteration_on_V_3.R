
# incorporating gamma convergence in before each V
# V generally converges within 3 steps (tol = 1e-5)

iteration_V <- function(NuisanceFit_folds, TargetFit_folds) {
  
  gamma_reached_max <- FALSE # indicator variable for hit_max_gamma in iteration_gamma()
  no_drop <- 3 # plateau length: stops if does not decrease for 3 consecutive steps
  
  # storage:
  best_V      <- Inf # stores the smallest possible so far
  best_gamma  <- gamma0  
  best_beta   <- matrix(NA, nrow = dimS, ncol = 1)
  
  P_hat_folds <- array(NA, dim = c(dimS, dimS, K))
  Q_hat_folds <- matrix(NA, nrow=dimS, ncol=K)
  V           <- numeric(K)
  tau_numr    <- numeric(K)
  tau_denom   <- numeric(K)
  
  # initiate:
  gamma <- gamma0
  prev_avg_V  <- Inf
  plateau_cnt <- 0
  
  for (iter in 1:max_iter_V) {
    
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
    iteration_gamma_result = iteration_gamma(NuisanceFit_folds, TargetFit_folds, beta_opt, gamma)
    gamma_new <- iteration_gamma_result$gamma
    gamma_reached_max <- gamma_reached_max || iteration_gamma_result$hit_max_gamma # once any gamma iteration hits max, stored as TRUE
    
    # V, ha, xi:
    for (j in 1:K){
      V[j] <- Vab(NuisanceFit_folds[[j]], TargetFit_folds[[j]], gamma_new, beta_opt)
    }

    avg_V <- mean(V) # then evaluate the convergence of avg_V
    
    # ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
    #cat(sprintf(
    #  "\niter %3d | V = %.6f  ΔV = %.6f  plateau_cnt = %d\n",
    #  iter, avg_V, avg_V - prev_avg_V, plateau_cnt))
    #
    #cat("gamma =")
    #print(round(gamma_new, 6), quote = FALSE)     # 4 × 1 layout
    #
    #cat("beta  =")
    #print(round(beta_opt, 6), quote = FALSE)      # 2 × 1 layout
    #cat("------------------------------------------------------------\n")
    # ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
    
    # store best positive V so far
    if (avg_V > 0 && avg_V < best_V) {
      best_V     <- avg_V             
      best_gamma <- gamma_new         
      best_beta  <- beta_opt          
    }
    
    # stop immediately if V ≤ 0
    if (avg_V <= 0) {
      message("V became non-positive; rolling back to the previous best V.")
      break # once break, output the best_XXXs
    }
    
    # stop if non-decrease for 3 consecutive rounds:

    if (avg_V >= prev_avg_V - .Machine$double.eps) { # if not decreasing
      plateau_cnt <- plateau_cnt + 1 # considering only consecutive
    } else {
      plateau_cnt <- 0
    }
    
    delta_V <- abs(avg_V - prev_avg_V)
    
    if (plateau_cnt >= no_drop) {
      break # do NOT touch best_XXX → roll back to global minimum
    }
    
    if (delta_V < tol) {
      best_V     <- avg_V
      best_gamma <- gamma_new
      best_beta  <- beta_opt
      break
    }
    
    prev_avg_V <- avg_V
    gamma      <- gamma_new
    
  }

  beta_opt <- best_beta # used in tau calculation later
  
  hit_V_max <- as.integer(iter == max_iter_V) # V only does one set of iteration, so compare with max_iter_V once (differentiate from gamma)
  hit_gamma_max <- as.integer(gamma_reached_max) 
  
  for (j in 1:K){
    tau_numr[j] <- Tau(NuisanceFit_folds[[j]], TargetFit_folds[[j]], beta_opt)$numr
    tau_denom[j] <- Tau(NuisanceFit_folds[[j]], TargetFit_folds[[j]], beta_opt)$denom
  }
  
  tau <- 1-mean(tau_numr)/mean(tau_denom)
  return(list(
    tau = tau,
    beta_opt = beta_opt,
    gamma = gamma,
    hit_gamma_max = hit_gamma_max,
    hit_V_max     = hit_V_max
         ))
  
}
