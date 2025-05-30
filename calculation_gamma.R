
dimX <- 3

new_gamma <- function(NuisanceFit, TargetFit, gamma, beta_opt) {
  
  Xnew=NuisanceFit$Xnew
  Xnew_i = cbind(X0 = 1, Xnew)
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
  
  Xt = TargetFit$Xt
  Xt_i = cbind(X0 = 1, Xt)
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
  
  # c1 debias:
  c1_e_d <- (as.vector(Y1new - m1_pred))*(S1new-r1_pred) - c1_pred_e
  c0_e_d <- (as.vector(Y0new - m0_pred))*(S0new-r0_pred) - c0_pred_e
  d1_e_d <- t(apply(S1new-r1_pred, 1, function(v) as.vector(v %o% v))) - d1_pred_e
  d0_e_d <- t(apply(S0new-r0_pred, 1, function(v) as.vector(v %o% v))) - d0_pred_e
  
  # calculation in parts
  b1_target = as.vector(exp(Xt_i %*% gamma)*b1_t_e)
  b1_source = as.vector(exp(Xnew_i %*% gamma)*omega*((Y1new-m1_pred)^2 - b1_pred_e))
  b0_target = as.vector(exp(-Xt_i %*% gamma)*b0_t_e)
  b0_source = as.vector(exp(-Xnew_i %*% gamma)*omega*((Y0new-m0_pred)^2 - b0_pred_e))
  c1_target = as.vector(exp(Xt_i %*% gamma)*c1_t_e %*% beta_opt)
  c1_source = as.vector(exp(Xnew_i %*% gamma)*omega*c1_e_d %*% beta_opt)
  c0_target = as.vector(exp(-Xt_i %*% gamma)*c0_t_e %*% beta_opt)
  c0_source = as.vector(exp(-Xnew_i %*% gamma)*omega*c0_e_d %*% beta_opt)
  d1_target = as.vector(exp(Xt_i %*% gamma)*(d1_t_e %*% c(beta_opt[1]^2,beta_opt[1]*beta_opt[2],beta_opt[1]*beta_opt[2],beta_opt[2]^2)))
  d1_source = as.vector(exp(Xnew_i %*% gamma)*omega*(d1_e_d %*% c(beta_opt[1]^2,beta_opt[1]*beta_opt[2],beta_opt[1]*beta_opt[2],beta_opt[2]^2)))
  d0_target = as.vector(exp(-Xt_i %*% gamma)*(d0_t_e %*% c(beta_opt[1]^2,beta_opt[1]*beta_opt[2],beta_opt[1]*beta_opt[2],beta_opt[2]^2)))
  d0_source = as.vector(exp(-Xnew_i %*% gamma)*omega*(d0_e_d %*% c(beta_opt[1]^2,beta_opt[1]*beta_opt[2],beta_opt[1]*beta_opt[2],beta_opt[2]^2)))
  
  # xi_debias:
  xi <- 0 +
    # b1:
    (colMeans(b1_target*Xt_i) + colMeans((b1_source*Xnew_i)[Anew == 1,,drop = FALSE]) 
    ) +
    # b0:
    -(colMeans(b0_target*Xt_i) + colMeans((b0_source*Xnew_i)[Anew == 0,,drop = FALSE])
    ) +
    # c1:
    (-2)*(colMeans(c1_target*Xt_i) + colMeans((c1_source*Xnew_i)[Anew == 1,,drop = FALSE])
    ) +
    # c0:
    2*(colMeans(c0_target*Xt_i) + colMeans((c0_source*Xnew_i)[Anew == 0,,drop = FALSE])
    ) + 
    # d1:
    (colMeans(d1_target*Xt_i) + colMeans((d1_source*Xnew_i)[Anew == 1,,drop = FALSE])
    ) +
    # d0:
    -(colMeans(d0_target*Xt_i) + colMeans((d0_source*Xnew_i)[Anew == 0,,drop = FALSE])
    )
  
  # h_debias:
  ha <- 0 +
    # b1:
    (t((b1_target*Xt_i)) %*% Xt_i / nrow(Xt_i) + 
       t((b1_source*Xnew_i)[Anew == 1,,drop = FALSE])%*%(Xnew_i[Anew == 1,,drop = FALSE]) / nrow(Xnew_i[Anew == 1,,drop = FALSE])
    ) + # source
    # b0:
    (t((b0_target*Xt_i)) %*% Xt_i / nrow(Xt_i) + 
        t((b0_source*Xnew_i)[Anew == 0,,drop = FALSE])%*%(Xnew_i[Anew == 0,,drop = FALSE]) / nrow(Xnew_i[Anew == 0,,drop = FALSE])
    ) +
    # c1:
    (-2)*(
      t((c1_target*Xt_i)) %*% Xt_i / nrow(Xt_i) + 
        t((c1_source*Xnew_i)[Anew == 1,,drop = FALSE])%*%(Xnew_i[Anew == 1,,drop = FALSE]) / nrow(Xnew_i[Anew == 1,,drop = FALSE])
    ) +
    # c0:
    (-2)*(
      t((c0_target*Xt_i)) %*% Xt_i / nrow(Xt_i) + 
        t((c0_source*Xnew_i)[Anew == 0,,drop = FALSE])%*%(Xnew_i[Anew == 0,,drop = FALSE]) / nrow(Xnew_i[Anew == 0,,drop = FALSE])
    ) + 
    # d1:
    (t((d1_target*Xt_i)) %*% Xt_i / nrow(Xt_i) + 
       t((d1_source*Xnew_i)[Anew == 1,,drop = FALSE])%*%(Xnew_i[Anew == 1,,drop = FALSE]) / nrow(Xnew_i[Anew == 1,,drop = FALSE])
    ) +
    # d0:
    (t((d0_target*Xt_i)) %*% Xt_i / nrow(Xt_i) + 
        t((d0_source*Xnew_i)[Anew == 0,,drop = FALSE])%*%(Xnew_i[Anew == 0,,drop = FALSE]) / nrow(Xnew_i[Anew == 0,,drop = FALSE])
    )
  
  # new gamma: gamma - xi/ja
  list(xi=xi, ha=ha)
  
}
