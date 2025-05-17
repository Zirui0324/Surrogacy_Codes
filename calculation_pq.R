
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
  b0_pred=NuisanceFit$b0_pred
  b1_pred=NuisanceFit$b1_pred
  c0_pred=NuisanceFit$c0_pred
  c1_pred=NuisanceFit$c1_pred
  d0_pred=NuisanceFit$d0_pred
  d1_pred=NuisanceFit$d1_pred
  alpha_k = cbind(X0 = 1, Xnew) %*% gamma
  
  Xt = TargetFit$Xt
  m0_t=TargetFit$m0_t
  m1_t=TargetFit$m1_t
  r0_t=TargetFit$r0_t
  r1_t=TargetFit$r1_t
  b0_t=TargetFit$b0_t
  b1_t=TargetFit$b1_t
  c0_t=TargetFit$c0_t
  c1_t=TargetFit$c1_t
  d0_t=TargetFit$d0_t
  d1_t=TargetFit$d1_t
  alpha_t = cbind(X0 = 1, Xt) %*% gamma
  
  dims <- ncol(TargetFit$r0_t)
  
  # unconditional debiased:
  D0 <- matrix(colMeans(d0_t), nrow = dims, ncol = dims) + # target # D0, D1: stored as p*p matrix; other nuisance as numeric
    matrix(colMeans(((t(apply(S0new, 1, function(v) as.vector(outer(v,v)))) - d0_pred)*omega)[Anew == 0,,drop = FALSE]),
           nrow = dims, ncol = dims) # source k-th fold
  D1 <- matrix(colMeans(d1_t), nrow = dims, ncol = dims) +
    matrix(colMeans(((t(apply(S1new, 1, function(v) as.vector(outer(v,v)))) - d1_pred)*omega)[Anew == 1,,drop = FALSE]),
           nrow = dims, ncol = dims)
  R0 <- colMeans(r0_t) +
    colMeans(((S0new - r0_pred)*omega)[Anew == 0,,drop = FALSE])
  R1 <- colMeans(r1_t) +
    colMeans(((S1new - r1_pred)*omega)[Anew == 1,,drop = FALSE])
  C0 <- colMeans(c0_t) +
    colMeans(((S0new*as.vector(Y0new) - c0_pred)*omega)[Anew == 0,,drop = FALSE])
  C1 <- colMeans(c1_t) +
    colMeans(((S1new*as.vector(Y1new) - c1_pred)*omega)[Anew == 1,,drop = FALSE])
  M0 <- mean(m0_t) + mean(((Y0new - m0_pred)*omega)[Anew == 0,drop = FALSE])
  M1 <- mean(m1_t) + mean(((Y1new - m1_pred)*omega)[Anew == 1,drop = FALSE])
  B0 <- mean(b0_t) + mean((((Y0new)^2-b0_pred)*omega)[Anew == 0,drop = FALSE])
  B1 <- mean(b1_t) + mean((((Y1new)^2-b1_pred)*omega)[Anew == 1,drop = FALSE])
  
  # conditional debiased:
  ad1_target <- as.vector(alpha_t)*d1_t
  ad1_source <- (as.vector(alpha_k)*omega*(t(apply(S1new, 1, function(v) as.vector(outer(v,v)))) - d1_pred))
  ad0_target <- as.vector(1/alpha_t)*d0_t
  ad0_source <- (as.vector(1/alpha_k)*omega*(t(apply(S0new, 1, function(v) as.vector(outer(v,v)))) - d0_pred))
  
  ar1sq_target <- as.vector(alpha_t)*t(apply(r1_t, 1, function(s) as.vector(outer(s,s))))
  ar1sq_source <- as.vector(alpha_k)*omega*(
    t(sapply(seq_len(nrow(S1new)), function(i) {as.vector(outer(S1new[i, ]-r1_pred[i, ], r1_pred[i, ]))})) +
      t(sapply(seq_len(nrow(S1new)), function(i) {as.vector(outer(r1_pred[i, ], S1new[i, ]-r1_pred[i, ]))})))
  ar0sq_target <- as.vector(1/alpha_t)*t(apply(r0_t, 1, function(s) as.vector(outer(s,s))))
  ar0sq_source <- as.vector(1/alpha_k)*omega*(
    t(sapply(seq_len(nrow(S0new)), function(i) {as.vector(outer(S0new[i, ]-r0_pred[i, ], r0_pred[i, ]))})) +
      t(sapply(seq_len(nrow(S0new)), function(i) {as.vector(outer(r0_pred[i, ], S0new[i, ]-r0_pred[i, ]))})))
  
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
  
  ac1_target <- as.vector(alpha_t)*c1_t
  ac1_source <- (as.vector(alpha_k)*omega*(S1new*as.vector(Y1new)-c1_pred))
  ac0_target <- as.vector(1/alpha_t)*c0_t
  ac0_source <- as.vector(1/alpha_k)*omega*(S0new*as.vector(Y0new)-c0_pred)
  
  am1r1_target <- as.vector(alpha_t)*r1_t*as.vector(m1_t)
  am1r1_source <- as.vector(alpha_k)*omega*(r1_pred*as.vector(Y1new-m1_pred) + (S1new-r1_pred)*as.vector(m1_pred))
  am0r0_target <- as.vector(1/alpha_t)*r0_t*as.vector(m0_t)
  am0r0_source <- as.vector(1/alpha_k)*omega*(r0_pred*as.vector(Y0new-m0_pred) + (S0new-r0_pred)*as.vector(m0_pred))
  
  m1r0_target <- r0_t*as.vector(m1_t)
  m1r0_source <- colMeans((omega*r0_pred*as.vector(Y1new-m1_pred))[Anew == 1,,drop = FALSE]) + colMeans((omega*(S0new-r0_pred)*as.vector(m1_pred))[Anew == 0,,drop = FALSE])
  m0r1_target <- r1_t*as.vector(m0_t)
  m0r1_source <- colMeans((omega*r1_pred*as.vector(Y0new-m0_pred))[Anew == 0,,drop = FALSE]) + colMeans((omega*(S1new-r1_pred)*as.vector(m0_pred))[Anew == 1,,drop = FALSE])
  
  m1m0_target <- m1_t*m0_t
  m1m0_source <- mean((omega*(Y1new-m1_pred)*m0_pred)[Anew == 1,,drop = FALSE]) + mean((omega*(Y0new-m0_pred)*m1_pred)[Anew == 0,,drop = FALSE])
  
  ab1_target <- alpha_t*b1_t
  ab1_source <- (alpha_k*omega*(Y1new^2-b1_pred))
  ab0_target <- (1/alpha_t)*b0_t
  ab0_source <- ((1/alpha_k)*omega*(Y0new^2-b0_pred))
  
  am1sq_target <- alpha_t*m1_t^2
  am1sq_source <- (alpha_k*omega*(Y1new-m1_pred)*m1_pred)
  am0sq_target <- (1/alpha_t)*m0_t^2
  am0sq_source <- ((1/alpha_k)*omega*(Y0new-m0_pred)*m0_pred)
  
  # P hat:
  P <- D1 - outer(R1, R1) + D0 - outer(R0, R0) + outer(R1, R0) + outer(R0, R1) + # unconditional
    # + a*d1:
    matrix(colMeans(ad1_target) + colMeans(ad1_source[Anew == 1,,drop = FALSE]),
           nrow = dims, ncol = dims
    ) +
    # + 1/a*d0:
    matrix(colMeans(ad0_target) + colMeans(ad0_source[Anew == 0,,drop = FALSE]),
           nrow = dims, ncol = dims
    ) -
    # - a*r1r1:
    matrix(colMeans(ar1sq_target) + colMeans(ar1sq_source[Anew == 1,,drop = FALSE]),
           nrow = dims, ncol = dims
    ) -
    # - 1/a*r0r0:
    matrix(colMeans(ar0sq_target) + colMeans(ar0sq_source[Anew == 0,,drop = FALSE]),
           nrow = dims, ncol = dims
    ) -
    # -r1r0
    matrix(colMeans(r1r0_target) + r1r0_source,
           nrow = dims, ncol = dims
    ) -
    # -r0r1
    matrix(colMeans(r0r1_target) + r0r1_source,
           nrow = dims, ncol = dims
    )
  
  # Q hat (calculated as t(Q)):
  QT <- C1 - M1*R1 + C0 - M0*R0 + M1*R0 + M0*R1 +
    # + a*c1
    colMeans(ac1_target) + colMeans(ac1_source[Anew == 1,,drop = FALSE]) +
    # + 1/a*c0
    colMeans(ac0_target) + colMeans(ac0_source[Anew == 0,,drop = FALSE]) -
    # - a*m1r1
    (colMeans(am1r1_target) + colMeans(am1r1_source[Anew == 1,,drop = FALSE])
    ) -
    # -1/a*m0r0
    (colMeans(am0r0_target) + colMeans(am0r0_source[Anew == 0,,drop = FALSE])
    ) +
    # m1r0
    colMeans(m1r0_target) + m1r0_source
  +
    # m0r1
    colMeans(m0r1_target) + m0r1_source
  
  Q <- t(t(QT))
  
  # solve(P) %*% Q
  
  # Return them
  list(P = P, Q = Q)
}
