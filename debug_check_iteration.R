
# ░░░░░░░░░░░░░
i = 10 # ░░░░░░
# ░░░░░░░░░░░░░

dimS <- 2
dimX <- 3
gamma0 <- matrix(0, dimX+1, 1)

set.seed(200427 + i)
tol        <- 0.00001
max_iter   <- 100
step_size  <- 1

DP <- GD(n) 
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
###-------------------------------------------- V --------------------------------------------###

iteration_V(NuisanceFit_folds, TargetFit_folds)

###------------------------------------------ gamma ------------------------------------------###
## storage:
#P_hat_folds <- array(NA, dim = c(dimS, dimS, K))
#Q_hat_folds <- matrix(NA, nrow=dimS, ncol=K)
## initiate:
#gamma <- gamma0
#
#  # P, Q:
#  for (j in 1:K) {
#    P_hat_folds[,,j]  <- calculation_pq(NuisanceFit_folds[[j]], TargetFit_folds[[j]], gamma)$P # can also move out
#    Q_hat_folds[,j]   <- calculation_pq(NuisanceFit_folds[[j]], TargetFit_folds[[j]], gamma)$Q
#  }
#  
#  # Beta: (average P, Q over folds)
#  avg_P    <- apply(P_hat_folds, c(1, 2), mean, na.rm = TRUE)
#  avg_Q    <- rowMeans(Q_hat_folds, na.rm = TRUE)
#  beta_opt <- solve(avg_P) %*% avg_Q
#  
#  # Gamma iteration:
#  gamma_new <- iteration_gamma(NuisanceFit_folds, TargetFit_folds, beta_opt, gamma)






