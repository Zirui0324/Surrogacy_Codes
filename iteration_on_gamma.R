
# Gamma generally converges within 4 steps (tol = 1e-5)

iteration_gamma <- function(NuisanceFit_folds, TargetFit_folds, beta_opt, gamma){
  
  hit_max_gamma <- FALSE # indicator variable, indicates whether max_iter_gamma is reached before stopping
  
  ha_folds    <- array(NA, dim = c(dimX+1, dimX+1, K))
  xi_folds    <- matrix(NA, nrow=dimX+1, ncol=K)
  
  for (iter in 1:max_iter_gamma) {
    
    for (j in 1:K){
      ha_folds[,,j] <- new_gamma(NuisanceFit_folds[[j]], TargetFit_folds[[j]], gamma, beta_opt)$ha
      xi_folds[,j]  <- new_gamma(NuisanceFit_folds[[j]], TargetFit_folds[[j]], gamma, beta_opt)$xi
    }
    avg_ha    <- apply(ha_folds, c(1, 2), mean, na.rm = TRUE)
    avg_xi    <- rowMeans(xi_folds, na.rm = TRUE)
    gamma_new <- gamma - step_size*(solve(avg_ha) %*% avg_xi)
    
    delta_gamma <- max(abs(gamma_new - gamma), na.rm = TRUE)
    
    # ░░░░░░░░░░░░░░░░░░░░░░░░░ display block: γ  ░░░░░░░░░░░░░░░░░░░░░░░░░░
    #cat(sprintf(
    #  "\niter %3d | |Δγ| = %.5f\n", iter, delta_gamma))
    #
    #cat("γ")
    #print(round(gamma_new, 6))                              # keeps 4×1 layout
    #
    #cat("----------------------------------------------------------------\n")
    # ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
    
    if (delta_gamma < tol) {
      gamma <- gamma_new
      #message("Converged!")
      break
    }
    gamma <- gamma_new
  }
  
  if (iter == max_iter_gamma) 
    hit_max_gamma <- TRUE
  
  list(gamma = gamma, hit_max_gamma = hit_max_gamma)
  
}







