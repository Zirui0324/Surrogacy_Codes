
calculation_pq <- function(NuisanceFit, TargetFit, gamma) {
  
  Xnew=NuisanceFit$Xnew
  Anew=NuisanceFit$Anew
  omega=NuisanceFit$omega
  Y0new=NuisanceFit$Y0new
  Y1new=NuisanceFit$Y1new
  S0new=NuisanceFit$S0new
  S1new=NuisanceFit$S1new
  m0_pred=NuisanceFit$m0_pred
  m1_pred=NuisanceFit$m1_pred
  r0_pred=NuisanceFit$r0_pred
  r1_pred=NuisanceFit$r1_pred
  #b0_pred=NuisanceFit$b0_pred
  #b1_pred=NuisanceFit$b1_pred
  c0_pred=NuisanceFit$c0_pred
  c1_pred=NuisanceFit$c1_pred
  d0_pred=NuisanceFit$d0_pred
  d1_pred=NuisanceFit$d1_pred
  #b0_pred_e=NuisanceFit$b0_pred_e
  #b1_pred_e=NuisanceFit$b1_pred_e
  c0_pred_e=NuisanceFit$c0_pred_e
  c1_pred_e=NuisanceFit$c1_pred_e
  d0_pred_e=NuisanceFit$d0_pred_e
  d1_pred_e=NuisanceFit$d1_pred_e
  alpha_k = exp(cbind(X0 = 1, Xnew) %*% gamma)
  
  Xt = TargetFit$Xt
  m0_t=TargetFit$m0_t
  m1_t=TargetFit$m1_t
  r0_t=TargetFit$r0_t
  r1_t=TargetFit$r1_t
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
  M0 <- mean(m0_t) + mean(((Y0new - m0_pred)*omega)[Anew == 0,drop = FALSE])
  M1 <- mean(m1_t) + mean(((Y1new - m1_pred)*omega)[Anew == 1,drop = FALSE])
  
  R0 <- colMeans(r0_t) +
    colMeans(((S0new - r0_pred)*omega)[Anew == 0,,drop = FALSE])
  R1 <- colMeans(r1_t) +
    colMeans(((S1new - r1_pred)*omega)[Anew == 1,,drop = FALSE])
  
  #B0 <- mean(b0_t) + mean((((Y0new)^2-b0_pred)*omega)[Anew == 0,drop = FALSE])
  #B1 <- mean(b1_t) + mean((((Y1new)^2-b1_pred)*omega)[Anew == 1,drop = FALSE])
  
  C0 <- colMeans(c0_t) +
    colMeans(((S0new*as.vector(Y0new) - c0_pred)*omega)[Anew == 0,,drop = FALSE])
  C1 <- colMeans(c1_t) +
    colMeans(((S1new*as.vector(Y1new) - c1_pred)*omega)[Anew == 1,,drop = FALSE])
  
  D0 <- matrix(colMeans(d0_t), nrow = dims, ncol = dims) + # target # D0, D1: stored as p*p matrix; other nuisance as numeric
    matrix(colMeans(((t(apply(S0new, 1, function(v) as.vector(outer(v,v)))) - d0_pred)*omega)[Anew == 0,,drop = FALSE]),
           nrow = dims, ncol = dims) # source k-th fold
  D1 <- matrix(colMeans(d1_t), nrow = dims, ncol = dims) +
    matrix(colMeans(((t(apply(S1new, 1, function(v) as.vector(outer(v,v)))) - d1_pred)*omega)[Anew == 1,,drop = FALSE]),
           nrow = dims, ncol = dims)
  
  #C0_e <- colMeans(c0_t_e) +
  #  colMeans((((S0new-r0_pred)*as.vector(Y0new-m0_pred) - c0_pred_e)*omega)[Anew == 0,,drop = FALSE])
  #C1_e <- colMeans(c1_t_e) +
  #  colMeans((((S1new-r1_pred)*as.vector(Y1new-m1_pred) - c1_pred_e)*omega)[Anew == 1,,drop = FALSE])
  #
  #D0_e <- matrix(colMeans(d0_t_e), nrow = dims, ncol = dims) + # target # D0, D1: stored as p*p matrix; other nuisance as numeric
  #  matrix(colMeans(((t(apply(S0new-r0_pred, 1, function(v) as.vector(outer(v,v)))) - d0_pred_e)*omega)[Anew == 0,,drop = FALSE]),
  #         nrow = dims, ncol = dims) # source k-th fold
  #
  #D1_e <- matrix(colMeans(d1_t_e), nrow = dims, ncol = dims) +
  #  matrix(colMeans(((t(apply(S1new-r1_pred, 1, function(v) as.vector(outer(v,v)))) - d1_pred_e)*omega)[Anew == 1,,drop = FALSE]),
  #         nrow = dims, ncol = dims)
  
  # conditional debiased: 
  r1r0_target <- t(sapply(seq_len(nrow(r1_t)), function(i) {as.vector(outer(r1_t[i, ], r0_t[i, ]))}))
  r1r0_source <- colMeans((omega*t(sapply(seq_len(nrow(S1new)), function(i)
  {as.vector(outer((S1new-r1_pred)[i,], r0_pred[i,]))})))[Anew == 1,,drop = FALSE]) + 
    colMeans((omega*t(sapply(seq_len(nrow(S0new)), function(i) 
    {as.vector(outer(r1_pred[i,], (S0new-r0_pred)[i,]))})))[Anew == 0,,drop = FALSE])
  
  r0r1_target <- t(sapply(seq_len(nrow(r0_t)), function(i) {as.vector(outer(r0_t[i, ], r1_t[i, ]))}))
  r0r1_source <- colMeans((omega*t(sapply(seq_len(nrow(S0new)), function(i)
  {as.vector(outer((S0new-r0_pred)[i,], r1_pred[i,]))})))[Anew == 0,,drop = FALSE]) + 
    colMeans((omega*t(sapply(seq_len(nrow(S1new)), function(i) 
    {as.vector(outer(r0_pred[i,], (S1new-r1_pred)[i,]))})))[Anew == 1,,drop = FALSE])
  
  m1r0_target <- r0_t*as.vector(m1_t)
  m1r0_source <- colMeans((omega*r0_pred*as.vector(Y1new-m1_pred))[Anew == 1,,drop = FALSE]) + colMeans((omega*(S0new-r0_pred)*as.vector(m1_pred))[Anew == 0,,drop = FALSE])
  m0r1_target <- r1_t*as.vector(m0_t)
  m0r1_source <- colMeans((omega*r1_pred*as.vector(Y0new-m0_pred))[Anew == 0,,drop = FALSE]) + colMeans((omega*(S1new-r1_pred)*as.vector(m0_pred))[Anew == 1,,drop = FALSE])
  
  # conditional debiased - error terms:
  ad0_e_target <- as.vector(1/alpha_t)*d0_t_e
  ad0_e_source <- as.vector(1/alpha_k)*omega*(t(apply(S0new-r0_pred, 1, function(v) as.vector(outer(v,v)))) - d0_pred_e)
  ad1_e_target <- as.vector(alpha_t)*d1_t_e
  ad1_e_source <- as.vector(alpha_k)*omega*(t(apply(S1new-r1_pred, 1, function(v) as.vector(outer(v,v)))) - d1_pred_e)
  
  ac0_e_target <- as.vector(1/alpha_t)*c0_t_e
  ac0_e_source <- as.vector(1/alpha_k)*omega*((S0new-r0_pred)*as.vector(Y0new-m0_pred)-c0_pred_e)
  ac1_e_target <- as.vector(alpha_t)*c1_t_e
  ac1_e_source <- as.vector(alpha_k)*omega*((S1new-r1_pred)*as.vector(Y1new-m1_pred)-c1_pred_e)
  
  # P hat:
  P <- D1 - outer(R1, R1) + D0 - outer(R0, R0) + outer(R1, R0) + outer(R0, R1) + 
    # + a*d1_e:
    matrix(colMeans(ad1_e_target) + colMeans(ad1_e_source[Anew == 1,,drop = FALSE]),
           nrow = dims, ncol = dims
    ) +
    # + 1/a*d0_e:
    matrix(colMeans(ad0_e_target) + colMeans(ad0_e_source[Anew == 0,,drop = FALSE]),
           nrow = dims, ncol = dims
    ) +
    # - r1r0
    (-1)*(
      matrix(colMeans(r1r0_target) + r1r0_source,
             nrow = dims, ncol = dims)
    ) +
    # - r0r1
    (-1)*(
      matrix(colMeans(r0r1_target) + r0r1_source,
             nrow = dims, ncol = dims)
    )
  
  # Q hat (calculated as t(Q)):
  QT <- C1 - M1*R1 + C0 - M0*R0 + M1*R0 + M0*R1 +
    # + a*c1_e
    colMeans(ac1_e_target) + colMeans(ac1_e_source[Anew == 1,,drop = FALSE]) +
    # + 1/a*c0_e
    colMeans(ac0_e_target) + colMeans(ac0_e_source[Anew == 0,,drop = FALSE]) +
    # -m1r0
    (-1)*(colMeans(m1r0_target) + m1r0_source) +
    # -m0r1
    (-1)*(colMeans(m0r1_target) + m0r1_source)
  
  
  Q <- t(t(QT))
  
  # solve(P) %*% Q
  
  # Return them
  list(P = P, Q = Q)
}
