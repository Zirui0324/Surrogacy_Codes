Tau <- function(NuisanceFit, TargetFit, beta_opt) {
  
  Mean <- function(x) {
    mean(x[is.finite(x)])
  }
  
  dims <- ncol(TargetFit$r0_t)
  
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
  b0_pred_e=NuisanceFit$b0_pred_e
  b1_pred_e=NuisanceFit$b1_pred_e
  c0_pred_e=NuisanceFit$c0_pred_e
  c1_pred_e=NuisanceFit$c1_pred_e
  d0_pred_e=NuisanceFit$d0_pred_e
  d1_pred_e=NuisanceFit$d1_pred_e
  # variance term:
  var_0 = b0_pred_e - 2*c0_pred_e %*% beta_opt + 
          d0_pred_e %*% c(beta_opt[1]^2,beta_opt[1]*beta_opt[2],beta_opt[1]*beta_opt[2],beta_opt[2]^2)
  var_1 = b1_pred_e - 2*c1_pred_e %*% beta_opt + 
    d1_pred_e %*% c(beta_opt[1]^2,beta_opt[1]*beta_opt[2],beta_opt[1]*beta_opt[2],beta_opt[2]^2)
  

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
  b0_t_e=TargetFit$b0_t_e
  b1_t_e=TargetFit$b1_t_e
  c0_t_e=TargetFit$c0_t_e
  c1_t_e=TargetFit$c1_t_e
  d0_t_e=TargetFit$d0_t_e
  d1_t_e=TargetFit$d1_t_e
  # variance term:
  var_0_t = b0_t_e - 2*c0_t_e %*% beta_opt + 
    d0_t_e %*% c(beta_opt[1]^2,beta_opt[1]*beta_opt[2],beta_opt[1]*beta_opt[2],beta_opt[2]^2)
  var_1_t = b1_t_e - 2*c1_t_e %*% beta_opt + 
    d1_t_e %*% c(beta_opt[1]^2,beta_opt[1]*beta_opt[2],beta_opt[1]*beta_opt[2],beta_opt[2]^2)
  
  # unconditional debiased:
  M0 <- mean(m0_t) + mean(((Y0new - m0_pred)*omega)[Anew == 0,drop = FALSE])
  M1 <- mean(m1_t) + mean(((Y1new - m1_pred)*omega)[Anew == 1,drop = FALSE])
  
  R0 <- colMeans(r0_t) +
    colMeans(((S0new - r0_pred)*omega)[Anew == 0,,drop = FALSE])
  R1 <- colMeans(r1_t) +
    colMeans(((S1new - r1_pred)*omega)[Anew == 1,,drop = FALSE])
  
  B0 <- mean(b0_t) + mean((((Y0new)^2-b0_pred)*omega)[Anew == 0,drop = FALSE])
  B1 <- mean(b1_t) + mean((((Y1new)^2-b1_pred)*omega)[Anew == 1,drop = FALSE])
  
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
  
  # conditional debiased:
  m1m0_target <- m1_t*m0_t
  m1m0_source <- mean((omega*(Y1new-m1_pred)*m0_pred)[Anew == 1,drop = FALSE]) + mean((omega*(Y0new-m0_pred)*m1_pred)[Anew == 0,drop = FALSE])
  
  m1r0_target <- r0_t*as.vector(m1_t)
  m1r0_source <- colMeans((omega*r0_pred*as.vector(Y1new-m1_pred))[Anew == 1,,drop = FALSE]) + colMeans((omega*(S0new-r0_pred)*as.vector(m1_pred))[Anew == 0,,drop = FALSE])
  m0r1_target <- r1_t*as.vector(m0_t)
  m0r1_source <- colMeans((omega*r1_pred*as.vector(Y0new-m0_pred))[Anew == 0,,drop = FALSE]) + colMeans((omega*(S1new-r1_pred)*as.vector(m0_pred))[Anew == 1,,drop = FALSE])
  
  r1r0_target <- t(sapply(seq_len(nrow(r1_t)), function(i) {as.vector(outer(r1_t[i, ], r0_t[i, ]))}))
  r1r0_source <- colMeans((omega*t(sapply(seq_len(nrow(S1new)), function(i)
  {as.vector(outer((S1new-r1_pred)[i,], r0_pred[i,]))})))[Anew == 1,,drop = FALSE]) + 
    colMeans((omega*t(sapply(seq_len(nrow(S0new)), function(i) 
    {as.vector(outer(r1_pred[i,], (S0new-r0_pred)[i,]))})))[Anew == 0,,drop = FALSE])
      
  var_e_target <- sqrt(var_1_t*var_0_t)
  var_e_source_0 <- omega*sqrt(var_1/var_0)*((Y0new - S0new %*% beta_opt - m0_pred + r0_pred %*% beta_opt)^2 - var_0)
  var_e_source_1 <- omega*sqrt(var_0/var_1)*((Y1new - S1new %*% beta_opt - m1_pred + r1_pred %*% beta_opt)^2 - var_1)
  
  var_target <- sqrt(b1_t_e)*sqrt(b0_t_e)
  var_source_0 <- omega*sqrt(b1_pred_e/b0_pred_e)*((Y0new-m0_pred)^2 - b0_pred_e)
  var_source_1 <- omega*sqrt(b0_pred_e/b1_pred_e)*((Y1new-m1_pred)^2 - b1_pred_e)
  
  # CALCULATION:
  numr <- B1 - 2*C1 %*% beta_opt +  t(D1 %*% beta_opt) %*% beta_opt - M1^2 + 2*M1*t(R1) %*% beta_opt - (t(R1) %*% beta_opt)^2 +
    B0 - 2*C0 %*% beta_opt +  t(D0 %*% beta_opt) %*% beta_opt - M0^2 + 2*M0*t(R0) %*% beta_opt - (t(R0) %*% beta_opt)^2 +
    2*M1*M0 - 2*M1*t(R0) %*% beta_opt - 2*M0*t(R1) %*% beta_opt + 2*(R1 %*% beta_opt)*(R0 %*% beta_opt) +
    (-2)*(
      mean(m1m0_target) + m1m0_source - 
        (colMeans(m1r0_target) + m1r0_source) %*% beta_opt -
        (colMeans(m0r1_target) + m0r1_source) %*% beta_opt +
        t(matrix(colMeans(r1r0_target) + r1r0_source,
                 nrow = dims, ncol = dims
        ) %*% beta_opt) %*% beta_opt
    ) +
    2*(
      Mean(var_e_target) +
        0.5*(Mean(var_e_source_0[Anew == 0,drop = FALSE])) +
        0.5*(Mean(var_e_source_1[Anew == 1,drop = FALSE]))
    )
  
  denom <- B1 - M1^2 + B0 - M0^2 + 2*M1*M0 +
    # -2*m1m0
    (-2)*(mean(m1m0_target) + m1m0_source) +
    # -2*√√
    (-2)*(
      Mean(var_target) +
        0.5*Mean(var_source_0[Anew == 0,drop = FALSE]) + # drop NaN and Inf
        0.5*Mean(var_source_1[Anew == 1,drop = FALSE])
    )
  
  list(
    numr = numr,
    denom = denom)
  
}




