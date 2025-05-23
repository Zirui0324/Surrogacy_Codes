
calculation_pq <- function(NuisanceFit, TargetFit, gamma) {
  
  Xnew=NuisanceFit$Xnew
  Anew=NuisanceFit$Anew
  omega=NuisanceFit$omega
  Y0new=NuisanceFit$Y0new
  Y1new=NuisanceFit$Y1new
  S0new=NuisanceFit$S0new
  S1new=NuisanceFit$S1new
  #m0_pred=NuisanceFit$m0_pred
  #m1_pred=NuisanceFit$m1_pred
  #r0_pred=NuisanceFit$r0_pred
  #r1_pred=NuisanceFit$r1_pred
  #b0_pred=NuisanceFit$b0_pred
  #b1_pred=NuisanceFit$b1_pred
  c0_pred=NuisanceFit$c0_pred
  c1_pred=NuisanceFit$c1_pred
  d0_pred=NuisanceFit$d0_pred
  d1_pred=NuisanceFit$d1_pred
  b0_pred_e=NuisanceFit$b0_pred_e
  b1_pred_e=NuisanceFit$b1_pred_e
  c0_pred_e=NuisanceFit$c0_pred_e
  c1_pred_e=NuisanceFit$c1_pred_e
  d0_pred_e=NuisanceFit$d0_pred_e
  d1_pred_e=NuisanceFit$d1_pred_e
  alpha_k = exp(cbind(X0 = 1, Xnew) %*% gamma)
  
  Xt = TargetFit$Xt
  #m0_t=TargetFit$m0_t
  #m1_t=TargetFit$m1_t
  #r0_t=TargetFit$r0_t
  #r1_t=TargetFit$r1_t
  #b0_t=TargetFit$b0_t
  #b1_t=TargetFit$b1_t
  c0_t=TargetFit$c0_t
  c1_t=TargetFit$c1_t
  d0_t=TargetFit$d0_t
  d1_t=TargetFit$d1_t
  #b0_t_e=TargetFit$b0_t_e
  #b1_t_e=TargetFit$b1_t_e
  c0_t_e=TargetFit$c0_t_e
  c1_t_e=TargetFit$c1_t_e
  d0_t_e=TargetFit$d0_t_e
  d1_t_e=TargetFit$d1_t_e
  alpha_t = exp(cbind(X0 = 1, Xt) %*% gamma)
  
  dims <- ncol(TargetFit$r0_t)
  
  # unconditional debiased:
  D0_e <- matrix(colMeans(d0_t_e), nrow = dims, ncol = dims) + # target # D0, D1: stored as p*p matrix; other nuisance as numeric
    matrix(colMeans(((t(apply(S0new-r0_pred, 1, function(v) as.vector(outer(v,v)))) - d0_pred_e)*omega)[Anew == 0,,drop = FALSE]),
           nrow = dims, ncol = dims) # source k-th fold
  
  D1_e <- matrix(colMeans(d1_t_e), nrow = dims, ncol = dims) +
    matrix(colMeans(((t(apply(S1new-r1_pred, 1, function(v) as.vector(outer(v,v)))) - d1_pred_e)*omega)[Anew == 1,,drop = FALSE]),
           nrow = dims, ncol = dims)
  
  C0_e <- colMeans(c0_t_e) +
    colMeans((((S0new-r0_pred)*as.vector(Y0new-m0_pred) - c0_pred_e)*omega)[Anew == 0,,drop = FALSE])
  C1_e <- colMeans(c1_t_e) +
    colMeans((((S1new-r1_pred)*as.vector(Y1new-m1_pred) - c1_pred_e)*omega)[Anew == 1,,drop = FALSE])
  
  # conditional debiased:
  ad1_e_target <- as.vector(alpha_t)*d1_t_e
  ad1_e_source <- as.vector(alpha_k)*omega*(t(apply(S1new-r1_pred, 1, function(v) as.vector(outer(v,v)))) - d1_pred_e)
  ad0_e_target <- as.vector(1/alpha_t)*d0_t_e
  ad0_e_source <- as.vector(1/alpha_k)*omega*(t(apply(S0new-r0_pred, 1, function(v) as.vector(outer(v,v)))) - d0_pred_e)
  
  ac1_e_target <- as.vector(alpha_t)*c1_t_e
  ac1_e_source <- as.vector(alpha_k)*omega*((S1new-r1_pred)*as.vector(Y1new-m1_pred)-c1_pred_e)
  ac0_e_target <- as.vector(1/alpha_t)*c0_t
  ac0_e_source <- as.vector(1/alpha_k)*omega*((S0new-r0_pred)*as.vector(Y0new-m0_pred)-c0_pred_e)
  
  # P hat:
  P <- D1_e + D0_e +
    # + a*d1_e:
    matrix(colMeans(ad1_e_target) + colMeans(ad1_e_source[Anew == 1,,drop = FALSE]),
           nrow = dims, ncol = dims
    ) +
    # + 1/a*d0_e:
    matrix(colMeans(ad0_e_target) + colMeans(ad0_e_source[Anew == 0,,drop = FALSE]),
           nrow = dims, ncol = dims
    ) 
  
  # Q hat (calculated as t(Q)):
  QT <- C1_e + C0_e +
    # + a*c1_e
    colMeans(ac1_e_target) + colMeans(ac1_e_source[Anew == 1,,drop = FALSE]) +
    # + 1/a*c0_e
    colMeans(ac0_e_target) + colMeans(ac0_e_source[Anew == 0,,drop = FALSE])
  Q <- t(t(QT))
  
  # solve(P) %*% Q
  
  # Return them
  list(P = P, Q = Q)
}
