
sim_function <- function(n, K, givenestimator, NBoot) {
    
    DP <- GD(n, a, b, rho_y, rho_s1, rho_s2)
    Index <- Split2(K, n) # sample split
    dimS <- ncol(DP$S0) # dim(theta) = dim(S)
    
    NuisanceFit_folds <- vector("list", K)
    TargetFit_folds   <- vector("list", K)
    # Bootstrap storage:
    NuisanceBoot_folds <- lapply(seq_len(NBoot), function(...) vector("list", K))
    TargetBoot_folds <- lapply(seq_len(NBoot), function(...) vector("list", K))
    tau_boot = numeric(NBoot)
    
    dimS = 2
    dimX = 3
    
    ###----------------------------------------- Nuisance Estimate -----------------------------------------###
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
          )
        }
      )
      
      # point estimate:
      NuisanceFit_folds[[j]]  <- result_j$NuisanceFit
      TargetFit_folds[[j]]    <- result_j$TargetFit
      
    }
    
    ###----------------------------------------- Iteration -----------------------------------------###
    
    # true tau:
    
    iteration_result = iteration_V(NuisanceFit_folds, TargetFit_folds)
  
    tau = iteration_result$tau
    beta_opt = iteration_result$beta_opt
    gamma = iteration_result$gamma
    hit_gamma_max = iteration_result$hit_gamma_max
    hit_V_max = iteration_result$hit_V_max
    
    message(sprintf("[%s]  true_tau obtained — start bootstrapping",
                    format(Sys.time(), "%m‑%d %H:%M")))
    
    # se(tau) through bootstrapping:
    for (i in 1:NBoot) {
      
      for (j in 1:K) {
        NuisanceBoot_folds[[i]][[j]] = resample_data(NuisanceFit_folds[[j]], TargetFit_folds[[j]])$NuisanceBoot
        TargetBoot_folds[[i]][[j]] = resample_data(NuisanceFit_folds[[j]], TargetFit_folds[[j]])$TargetBoot
      }
      
      tau_boot[i] <- tryCatch( # store as NA if errors occur (doesn't converge/haissen)
        iteration_V(NuisanceBoot_folds[[i]], TargetBoot_folds[[i]])$tau,
        error = function(e) NA_real_ 
      )
      
    }
    
    se_tau <- sd(tau_boot[is.finite(tau_boot)], na.rm = TRUE)
    
    ###----------------------------------------- Output -----------------------------------------###
    
    return(list(
      tau = tau,
      se_tau = se_tau,
      beta = beta_opt,
      gamma = gamma,
      hit_gamma_max = hit_gamma_max,
      hit_V_max = hit_V_max
    ))
}
