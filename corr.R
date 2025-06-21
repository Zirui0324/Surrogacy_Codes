max_iter_V   <- 10
max_iter_gamma   <- 100
tol        <- 0.00001
step_size  <- 1

p <- 5
n <- 5000
rho_vals <- seq(-1, 1, by = 0.2)
grid <- expand.grid(rho_y = rho_vals, rho_s = rho_vals)

#  ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░ change corr(S,Y), store tau, beta estimate ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
results <- vector("list", nrow(grid))

for (g in seq_len(nrow(grid))) {
  
  set.seed(200427+g)
  
  rho_y  <- grid$rho_y[g]
  rho_s1 <- grid$rho_s[g] # rho_s1 = rho_s2
  rho_s2 <- grid$rho_s[g]
  
  tau_vec   <- numeric(p)
  beta1_vec <- numeric(p)
  beta2_vec <- numeric(p)
  
  for (k in seq_len(p)) {
    
    cat(sprintf("ρ_y = %4.1f   ρ_s = %4.1f   ρ_s = %4.1f   |  rep %2d / %2d\n",
                rho_y, rho_s1, rho_s2, k, p))        # ← updated message
    
    
    ## --- generate data with correlated errors
    DP     <- GD(n, a, b, rho_y, rho_s1, rho_s2)
    Index  <- Split2(K, n)
    
    ## --- K-fold fitting
    NuisanceFit_folds <- vector("list", K)
    TargetFit_folds   <- vector("list", K)
    
    for (j in seq_len(K)) {
      idNj <- Index != j
      idj  <- !idNj
      
      In  <- lapply(DP, function(z) if (is.matrix(z)) z[idNj, , drop = FALSE] else z[idNj])
      Out <- lapply(DP, function(z) if (is.matrix(z)) z[idj , , drop = FALSE] else z[idj ])
      Target <- DP$Xt
      
      fit_j <- tryCatch(
        estimate(In, Out, Target, givenestimator, NBoot),
        error = function(e) list(NuisanceFit = list(), TargetFit = list())
      )
      NuisanceFit_folds[[j]] <- fit_j$NuisanceFit
      TargetFit_folds[[j]]   <- fit_j$TargetFit
    }
    
    res <- tryCatch(
      iteration_V(NuisanceFit_folds, TargetFit_folds),
      error = function(e) NULL
    )
    
    if (is.null(res)) {
      tau_vec[k]   <- NA_real_
      beta1_vec[k] <- NA_real_
      beta2_vec[k] <- NA_real_
    } else {
      tau_vec[k]   <- res$tau
      beta1_vec[k] <- res$beta_opt[1]
      beta2_vec[k] <- res$beta_opt[2]
    }
    
  }
  
  ## save one row: ρ0, ρ1, τ₁ … τ₁₀
  results[[g]] <- c(rho_y = rho_y,
                    rho_s1 = rho_s1,
                    rho_s2 = rho_s2,
                    setNames(tau_vec,   paste0("tau_",   seq_len(p))),
                    setNames(beta1_vec, paste0("beta1_", seq_len(p))),
                    setNames(beta2_vec, paste0("beta2_", seq_len(p)))
                    )
  
  write.csv(do.call(rbind, lapply(results, as.data.frame.list)),
            "./output/corr/corr_results.csv",
            row.names = FALSE)
}


results_df <- do.call(rbind, lapply(results, as.data.frame.list))

# true lower bound 
corr = read_csv("./output/corr/corr_0620.csv") |> 
  as.data.frame() |> 
  mutate(rho_y = round(rho_y, 1),
         rho_s1 = round(rho_s1, 1),
         rho_s2 = round(rho_s2, 1))

summary_corr <- corr %>% 
  rowwise() %>%                           # process one row at a time
  mutate(
    ## τ statistics
    tau_mean  = mean(c_across(starts_with("tau_")),   na.rm = TRUE),
    tau_sd    =   sd(c_across(starts_with("tau_")),   na.rm = TRUE),
    ## β₁ statistics
    beta1_mean = mean(c_across(starts_with("beta1_")), na.rm = TRUE),
    beta1_sd   =   sd(c_across(starts_with("beta1_")), na.rm = TRUE),
    ## β₂ statistics
    beta2_mean = mean(c_across(starts_with("beta2_")), na.rm = TRUE),
    beta2_sd   =   sd(c_across(starts_with("beta2_")), na.rm = TRUE)
  ) %>% 
  ungroup() %>% 
  select(
    rho_y,
    rho_s1,
    rho_s2,
    tau_mean,  tau_sd,
    beta1_mean, beta1_sd,
    beta2_mean, beta2_sd
  )

#  ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░ change corr, truth ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
beta_s0 <- beta_s1 <- c(1, -1)
DP <- GD(n, a, b, rho_y, rho_s1, rho_s2) # doesn't matter, because only using X from Target below
Target <- rbind(DP$Xt)
grid_true <- expand.grid(rho_y = rho_vals, rho_s = rho_vals)

results_true <- data.frame() 

for (g in seq_len(nrow(grid_true))) {
  
  rho_y  <- grid_true$rho_y[g]
  rho_s1 <- grid_true$rho_s[g]
  rho_s2 <- grid_true$rho_s[g]

  
  ## ---- Y-errors: 2×2 block -------------------------------------------
  sigma_y <- matrix(c(
    b, rho_y*sqrt(b*b),
    rho_y*sqrt(b*b), b),
    nrow = 2, byrow = TRUE)
  
  errs_Y  <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = sigma_y)
  e_Y0    <- errs_Y[, 1]
  e_Y1    <- errs_Y[, 2]
  
  ## ---- S(·)_1 errors: 2×2 block --------------------------------------
  sigma_s1 <- matrix(c(
    a, rho_s1*sqrt(a*a),
    rho_s1*sqrt(a*a),  a),
    nrow = 2, byrow = TRUE)
  
  errs_S1 <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = sigma_s1)
  e_S0_1  <- errs_S1[, 1]
  e_S1_1  <- errs_S1[, 2]
  
  ## ---- S(·)_2 errors: 2×2 block --------------------------------------
  sigma_s2 <- matrix(c(
    a, rho_s2*sqrt(a*a),
    rho_s2*sqrt(a*a),  a),
    nrow = 2, byrow = TRUE)
  
  errs_S2 <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = sigma_s2)
  e_S0_2  <- errs_S2[, 1]
  e_S1_2  <- errs_S2[, 2]
  
  ## ---- assemble matrices used later ----------------------------------
  e_S0 <- cbind(e_S0_1, e_S0_2) 
  e_S1 <- cbind(e_S1_1, e_S1_2) 
  
  ## ---------- build T -----------------------------------------------
  T <- Target |>
    as.data.frame() |>
    mutate(
      S0_1 = -2*sin(X1) + 0.5*X3 + e_S0_1,
      S0_2 =  X2 + X1            + e_S0_2,
      
      S1_1 =  sin(X1) + X3       + e_S1_1,
      S1_2 = -2*X2   - X1        + e_S1_2,
      
      Y0 = S0_1*beta_s0[1] + S0_2*beta_s0[2] +
        (X1*X3 + 2*X1^2) + e_Y0,
      
      Y1 = S1_1*beta_s1[1] + S1_2*beta_s1[2] +
        (X2*X3 + X1*X2 + X1^2) + e_Y1) |> 
    mutate(
      delta_S1 = S1_1 - S0_1,
      delta_S2 = S1_2 - S0_2,
      delta_Y  = Y1   - Y0
    )
  
  ## ---------- θ_true  &  τ_true -------------------------------------
  fit          <- lm(delta_Y ~ 0 + delta_S1 + delta_S2, data = T)
  theta_true   <- coef(fit)
  
  var_n <- var(T$delta_Y -
                 (T$delta_S1*theta_true[1] + T$delta_S2*theta_true[2]))
  var_d <- var(T$delta_Y)
  tau_true <- 1 - var_n / var_d
  
  ## ---------- store & save ------------------------------------------
  new_row <- data.frame(
    rho_y = rho_y,
    rho_s1 = rho_s1,
    rho_s2 = rho_s2,
    tau_true = tau_true,
    beta1    = theta_true[1],
    beta2    = theta_true[2]
  )
  
  results_true <- rbind(results_true, new_row)
  write.csv(results_true, "./output/corr/corr_true.csv", row.names = FALSE)   # incremental save
  
  cat(sprintf("Done ρ0 = %4.1f, ρ1 = %4.1f  (row %d/%d)\n",
              rho0, rho1, g, nrow(grid_true)))
}

results_true = results_true |> 
  as.data.frame() |> 
  mutate(rho_y = round(rho_y, 1),
         rho_s1 = round(rho_s1, 1),
         rho_s2 = round(rho_s2, 1))

#  ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░


summary = summary_corr |> 
  left_join(results_true, by = c("rho_y", "rho_s1", "rho_s2")) |> 
  select(rho_y, rho_s1, rho_s2, tau_true, tau_mean, beta1, beta1_mean, beta2, beta2_mean) |> 
  mutate(tau_diff = tau_true-tau_mean)
