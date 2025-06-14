DP <- GD(n) 
Index <- Split2(K, n) # sample split
dimS <- ncol(DP$S0) # dim(beta) = dim(S)
dimX <- 3
gamma0 <- matrix(0, dimX+1, 1)

# array for storing B, C, D from the K folds
NuisanceFit_folds <- vector("list", K)
TargetFit_folds   <- vector("list", K)
# Bootstrap storage:
NuisanceBoot_folds <- lapply(seq_len(NBoot), function(...) vector("list", K))
TargetBoot_folds <- lapply(seq_len(NBoot), function(...) vector("list", K))
tau_boot = numeric(NBoot)

P_hat_folds <- array(NA, dim = c(dimS, dimS, K))
Q_hat_folds <- matrix(NA, nrow=dimS, ncol=K)
V <- numeric(K)
gamma_i <- matrix(NA, nrow=dimS, ncol=K)
ha_folds <- array(NA, dim = c(dimX+1, dimX+1, K))
xi_folds <- matrix(NA, nrow=dimX+1, ncol=K)
tau_numr <- numeric(K)
tau_denom <- numeric(K)

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

# true tau:
tau = iteration(NuisanceFit_folds, TargetFit_folds)

# se(tau) through bootstrapping:
for (i in 1:NBoot) {
  
  for (j in 1:K) {
    NuisanceBoot_folds[[i]][[j]] = resample_data(NuisanceFit_folds[[j]], TargetFit_folds[[j]])$NuisanceBoot
    TargetBoot_folds[[i]][[j]] = resample_data(NuisanceFit_folds[[j]], TargetFit_folds[[j]])$TargetBoot
  }
  
  tau_boot[i] = iteration(NuisanceBoot_folds[[i]], TargetBoot_folds[[i]])

}

ase_tau = sd(tau_boot)

