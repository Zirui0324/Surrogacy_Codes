# For Bootstrap
resample_data <- function(NuisanceFit, TargetFit) {
  
  n_src   <- length(NuisanceFit$m0_pred) 
  n_trg  <- length(TargetFit$m0_t) 
  
  # draw indices with replacement
  h <- n_src/2
  idx_src <- c(sample(h, h,TRUE), sample((h+1):n_src, n_src-h, TRUE) )
  
  idx_src  <- sample(seq_len(n_src),  n_src,  replace=TRUE)
  idx_trg <- sample(seq_len(n_trg), n_trg, replace=TRUE)
  
  # new bootstrap Out list
  NuisanceBoot <- list(
    Xnew = NuisanceFit$Xnew[idx_src, , drop=FALSE],
    Anew = NuisanceFit$Anew[idx_src],
    omega = NuisanceFit$omega[idx_src],
    Y0new = NuisanceFit$Y0new[idx_src], Y1new = NuisanceFit$Y1new[idx_src], 
    S0new = NuisanceFit$S0new[idx_src, , drop=FALSE], S1new = NuisanceFit$S1new[idx_src, , drop=FALSE], 
    m0_pred = NuisanceFit$m0_pred[idx_src], m1_pred = NuisanceFit$m1_pred[idx_src],
    r0_pred = NuisanceFit$r0_pred[idx_src, , drop=FALSE], r1_pred = NuisanceFit$r1_pred[idx_src, , drop=FALSE],
    b0_pred = NuisanceFit$b0_pred[idx_src], b1_pred = NuisanceFit$b1_pred[idx_src],
    c0_pred = NuisanceFit$c0_pred[idx_src, , drop=FALSE], c1_pred = NuisanceFit$c1_pred[idx_src, , drop=FALSE],
    d0_pred = NuisanceFit$d0_pred[idx_src, , drop=FALSE], d1_pred = NuisanceFit$d1_pred[idx_src, , drop=FALSE],
    # error terms:
    b0_pred_e = NuisanceFit$b0_pred_e,b1_pred_e = NuisanceFit$b1_pred_e[idx_src],
    c0_pred_e = NuisanceFit$c0_pred_e,c1_pred_e = NuisanceFit$c1_pred_e[idx_src, , drop=FALSE],
    d0_pred_e = NuisanceFit$d0_pred_e,d1_pred_e = NuisanceFit$d1_pred_e[idx_src, , drop=FALSE]
  )
  
  
  # new bootstrap Target list
  TargetBoot <- list(
    Xt  =  TargetFit$Xt[idx_src, , drop=FALSE],
    m0_t = TargetFit$m0_t[idx_src], m1_t = TargetFit$m1_t[idx_src], 
    r0_t = TargetFit$r0_t[idx_src, , drop=FALSE], r1_t = TargetFit$r1_t[idx_src, , drop=FALSE],
    b0_t = TargetFit$b0_t[idx_src], b1_t = TargetFit$b1_t[idx_src],
    c0_t = TargetFit$c0_t[idx_src, , drop=FALSE], c1_t = TargetFit$c1_t[idx_src, , drop=FALSE],
    d0_t = TargetFit$d0_t[idx_src, , drop=FALSE], d1_t = TargetFit$d1_t[idx_src, , drop=FALSE],
    # error terms:
    b0_t_e = TargetFit$b0_t_e[idx_src], b1_t_e = TargetFit$b1_t_e[idx_src],
    c0_t_e = TargetFit$c0_t_e[idx_src, , drop=FALSE], c1_t_e = TargetFit$c1_t_e[idx_src, , drop=FALSE],
    d0_t_e = TargetFit$d0_t_e[idx_src, , drop=FALSE], d1_t_e = TargetFit$d1_t_e[idx_src, , drop=FALSE]
  )
  
  # return
  list(NuisanceBoot = NuisanceBoot, TargetBoot = TargetBoot)
}




