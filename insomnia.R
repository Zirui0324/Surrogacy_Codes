library(tidyverse)
library(rje)
library(caret)
library(ks)
library(kernlab)
library(splines) # for pte


N=100
n=1000
sample_sizes = c(n)
K=5
givenestimator = 'gbm' # c('glm', 'svmLinear2', 'rf', 'nnet)
NBoot=100
dimX = 3
dimS = 2
gamma0 = matrix(0, dimX+1, 1)

# data set-up:
a = 0.5           # var S
b = 0.5           # var Y
rho_y  <- 0       # Corr(Y0 , Y1)
rho_s1 <- 0       # Corr(S0_1 , S1_1)
rho_s2 <- 0       # Corr(S0_2 , S1_2)
z=1


Split2 <- function(K, input) {
  I1 <- rep(0, input)        # initiate the fold labels storage for each data point
  half <- input / 2
  
  # random permutation of [1~half]
  firstHalfIDs <- sample(1:half, half, replace = FALSE)
  
  # number of points per fold in the first half
  m <- half / K
  
  # assign folds to the first half
  for (i in 1:K) {
    # The slice for fold i
    theseIDs <- firstHalfIDs[((i - 1) * m + 1) : (i * m)]
    I1[theseIDs] <- i
  }
  
  # random permutation of [half+1~input]
  secondHalfIDs <- sample((half + 1) : input, half, replace = FALSE)
  
  # assign folds to the second half
  for (i in 1:K) {
    theseIDs <- secondHalfIDs[((i - 1) * m + 1) : (i * m)]
    I1[theseIDs] <- i
  }
  
  return(I1)
}



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



estimate <- function(In, Out, Target, givenestimator, NBoot) {
  
  # In and Out both from source data, defined through sample splitting
  
  # In: for training model  #(1-1/K)*n dp
  X <- In$X
  S0 <- In$S0
  S1 <- In$S1
  Y0 <- In$Y0
  Y1 <- In$Y1
  A <- In$A
  # new added for CS bound:
  b0 <- Y0^2 # keeping both A=0/1 now but separating when doing regression
  b1 <- Y1^2
  c0 <- S0*as.vector(Y0)
  c1 <- S1*as.vector(Y1)
  d0 <- t(apply(S0, 1, function(v) as.vector(v %o% v))) # 1*(p^2)
  d1 <- t(apply(S1, 1, function(v) as.vector(v %o% v))) 
  
  # Out: for estimation  #(1/K)*n dp
  Xnew <- Out$X
  S0new <- Out$S0
  S1new <- Out$S1
  Y0new <- Out$Y0
  Y1new <- Out$Y1
  Anew <- Out$A
  
  # Target: X from target population, for theta estimation   # n dp
  Xt <- Target
  
  ######################## omega estimate ########################
  Xc <- data.frame(rbind(cbind(X, source = 0), cbind(Xt, source = 1))) # pT = p1
  
  ### 1. glm:
  #Xc$source <- as.factor(Xc$source)
  #p <- glm(source ~ ., data = Xc, family = binomial) # logistic regression
  #p1 <- predict(p, newdata = data.frame(Xnew), type = "response")
  #omega <- (p1/(1-p1))*(1-1/K)
  ##omega <- rep(1, length(omega)) # for true tau calculation #*#
  #
  #fit_model <- function(outcome, predictor) {
  #  train(outcome ~ ., data = data.frame(outcome, predictor),
  #        method = 'glm',
  #        trControl = trainControl(method = "cv", number = 5, verboseIter = FALSE)
  #  )
  #}
  
  ### 2. machine learning:
  Xc$source <- factor(Xc$source, 
                      levels = c(0,1), 
                      labels = c("Class0", "Class1"))
  
  fitControl_p <- trainControl( # for p
    method = "cv",
    number = 5,
    verboseIter=FALSE,
    classProbs  = TRUE) # for classification
  
  fitControl <- trainControl( # for nuisance models
    method = "cv",
    number = 5,
    verboseIter=FALSE)
  
  # nuisance:
  fit_model <- function(outcome, predictor) {
    train(outcome ~ ., data = data.frame(outcome, predictor),
          method = givenestimator,
          trControl = fitControl,
          verbose = FALSE
    )}
  
  # p & omega, for different algrthms:
  if(givenestimator %in% c('svmLinear', 'svmLinear2')){
    
    p <- train(source ~., data = Xc, 
               method = givenestimator, 
               trControl =fitControl_p,
               verbose = FALSE,
               probability = TRUE) # unique to svmLinear2
    
  }else{
    
    p <- train(source ~., data = Xc, 
               method = givenestimator, 
               trControl =fitControl_p,
               verbose = FALSE) # c('gbm','rf','nnet')
    
  }
  
  p1 <- predict(p, newdata = data.frame(Xnew), type = "prob")$Class1
  omega <- (p1/(1-p1))*(1-1/K)
  
  #initialize model-storing list for high-d S (p*1), c (p*1), and d (p*p)
  dims <- ncol(S0)
  dims2 <- dims^2
  
  # Separate data by A=0/1 for training purpose
  X_A0 <- X[A == 0, , drop = FALSE]
  X_A1 <- X[A == 1, , drop = FALSE]
  S0_train <- S0[A == 0, , drop = FALSE]
  S1_train <- S1[A == 1, , drop = FALSE]
  Y0_train <- Y0[A == 0, , drop = FALSE]
  Y1_train <- Y1[A == 1, , drop = FALSE]
  # new added for CS bound:
  b0_train <- b0[A == 0, , drop = FALSE]
  b1_train <- b1[A == 1, , drop = FALSE]
  c0_train <- c0[A == 0, , drop = FALSE]
  c1_train <- c1[A == 1, , drop = FALSE]
  d0_train <- d0[A == 0, , drop = FALSE]
  d1_train <- d1[A == 1, , drop = FALSE]
  
  # creat list for n-d estimators
  modelS0 <- list()
  modelS1 <- list()
  modelc0 <- list()
  modelc1 <- list()
  modeld0 <- vector("list", dims2)
  modeld1 <- vector("list", dims2)
  
  # Regression:
  fitY0 <- fit_model(Y0_train, X_A0) # Y0 ~ X_{A=0}
  fitY1 <- fit_model(Y1_train, X_A1) # Y1 ~ X_{A=1}
  fitb0 <- fit_model(b0_train, X_A0)
  fitb1 <- fit_model(b1_train, X_A1)
  
  for (i in seq_len(dims)) { # p
    
    S0_dim <- S0_train[, i]
    S1_dim <- S1_train[, i]
    c0_dim <- c0_train[, i]
    c1_dim <- c1_train[, i]
    
    modelS0[[i]] <- fit_model(S0_dim, X_A0)
    modelS1[[i]] <- fit_model(S1_dim, X_A1) # store models for S_i
    modelc0[[i]] <- fit_model(c0_dim, X_A0)
    modelc1[[i]] <- fit_model(c1_dim, X_A1) 
    
  }
  
  for (j in seq_len(dims2)) { # p^2
    
    d0_dim <- d0_train[, j]
    d1_dim <- d1_train[, j]
    
    modeld0[[j]] <- fit_model(d0_dim, X_A0)
    modeld1[[j]] <- fit_model(d1_dim, X_A1)
    
  }
  # Regression ends #
  
  # New Regression:
  # (train nuisance_error)
  
  # format outcome (using source[-k] data), regardless of A = 0/1
  b0_e <- (Y0 - predict(fitY0, newdata = X))^2 # true Y(a) - m(a) # mse/var = 0.3
  b1_e <- (Y1 - predict(fitY1, newdata = X))^2
  c0_e <- as.vector(Y0 - predict(fitY0, newdata = X))*(S0 - sapply(modelS0, predict, newdata = X))
  c1_e <- as.vector(Y1 - predict(fitY1, newdata = X))*(S1 - sapply(modelS1, predict, newdata = X))
  d0_e <- t(apply((S0 - sapply(modelS0, predict, newdata = X)), 1, function(v) as.vector(v %o% v))) # mse/var = 0.3
  d1_e <- t(apply((S1 - sapply(modelS1, predict, newdata = X)), 1, function(v) as.vector(v %o% v)))
  
  b0_e_train <- b0_e[A == 0, , drop = FALSE]
  b1_e_train <- b1_e[A == 1, , drop = FALSE]
  c0_e_train <- c0_e[A == 0, , drop = FALSE]
  c1_e_train <- c1_e[A == 1, , drop = FALSE]
  d0_e_train <- d0_e[A == 0, , drop = FALSE]
  d1_e_train <- d1_e[A == 1, , drop = FALSE]
  
  # dim = 1
  fitb0_e <- fit_model(b0_e_train, X_A0)
  fitb1_e <- fit_model(b1_e_train, X_A1)
  
  # dim = p
  modelc0_e <- list()
  modelc1_e <- list()
  modeld0_e <- vector("list", dims2)
  modeld1_e <- vector("list", dims2)
  
  for (i in seq_len(dims)) {
    
    c0_dim <- c0_e_train[, i]
    c1_dim <- c1_e_train[, i]
    
    modelc0_e[[i]] <- fit_model(c0_dim, X_A0)
    modelc1_e[[i]] <- fit_model(c1_dim, X_A1)
    
  }
  
  for (j in seq_len(dims2)) { # p^2
    
    d0_dim <- d0_e_train[, j]
    d1_dim <- d1_e_train[, j]
    
    modeld0_e[[j]] <- fit_model(d0_dim, X_A0)
    modeld1_e[[j]] <- fit_model(d1_dim, X_A1)
    
  }
  # New Regression ends #
  
  # NUISANCE: Prediction for Source k-th X (regardless of A)
  m0_pred <- predict(fitY0, newdata = Xnew)
  m1_pred <- predict(fitY1, newdata = Xnew)
  b0_pred <- predict(fitb0, newdata = Xnew)
  b1_pred <- predict(fitb1, newdata = Xnew)
  r0_pred <- sapply(modelS0, predict, newdata = Xnew)
  r1_pred <- sapply(modelS1, predict, newdata = Xnew)
  c0_pred <- sapply(modelc0, predict, newdata = Xnew)
  c1_pred <- sapply(modelc1, predict, newdata = Xnew)
  d0_pred <- sapply(modeld0, predict, newdata = Xnew)
  d1_pred <- sapply(modeld1, predict, newdata = Xnew)
  
  # NUISANCE: Prediction for Target X (A unknown)
  # m,r,c,d
  m0_t <- predict(fitY0, newdata = Xt)
  m1_t <- predict(fitY1, newdata = Xt)
  b0_t <- predict(fitb0, newdata = Xt) # not used in numerator
  b1_t <- predict(fitb1, newdata = Xt) # not used in numerator
  r0_t <- sapply(modelS0, predict, newdata = Xt)
  r1_t <- sapply(modelS1, predict, newdata = Xt)
  c0_t <- sapply(modelc0, predict, newdata = Xt)
  c1_t <- sapply(modelc1, predict, newdata = Xt)
  d0_t <- sapply(modeld0, predict, newdata = Xt)
  d1_t <- sapply(modeld1, predict, newdata = Xt)
  
  #  NUISANCE_ERROR: Prediction for Source k-th X (regardless of A)
  b0_pred_e <- predict(fitb0_e, newdata = Xnew)
  b1_pred_e <- predict(fitb1_e, newdata = Xnew)
  c0_pred_e <- sapply(modelc0_e, predict, newdata = Xnew)
  c1_pred_e <- sapply(modelc1_e, predict, newdata = Xnew)
  d0_pred_e <- sapply(modeld0_e, predict, newdata = Xnew)
  d1_pred_e <- sapply(modeld1_e, predict, newdata = Xnew)
  
  #  NUISANCE_ERROR: Prediction for Target X
  b0_t_e <- predict(fitb0_e, newdata = Xt)
  b1_t_e <- predict(fitb1_e, newdata = Xt)
  c0_t_e <- sapply(modelc0_e, predict, newdata = Xt)
  c1_t_e <- sapply(modelc1_e, predict, newdata = Xt)
  d0_t_e <- sapply(modeld0_e, predict, newdata = Xt)
  d1_t_e <- sapply(modeld1_e, predict, newdata = Xt)
  
  # Prepare data list with nuisance estimators for beta & V calculation
  NuisanceFit <- list(
    Xnew = Xnew,
    Anew = Anew,
    omega = omega,
    Y0new = Y0new, Y1new = Y1new, S0new = S0new, S1new = S1new, 
    m0_pred = m0_pred, m1_pred = m1_pred,
    r0_pred = r0_pred, r1_pred = r1_pred,
    b0_pred = b0_pred, b1_pred = b1_pred,
    c0_pred = c0_pred, c1_pred = c1_pred,
    d0_pred = d0_pred, d1_pred = d1_pred,
    # error terms:
    b0_pred_e = b0_pred_e,b1_pred_e = b1_pred_e,
    c0_pred_e = c0_pred_e,c1_pred_e = c1_pred_e,
    d0_pred_e = d0_pred_e,d1_pred_e = d1_pred_e
  )
  
  TargetFit <- list(
    Xt = Xt,
    m0_t = m0_t, m1_t = m1_t, 
    r0_t = r0_t, r1_t = r1_t,
    b0_t = b0_t, b1_t = b1_t,
    c0_t = c0_t, c1_t = c1_t,
    d0_t = d0_t, d1_t = d1_t,
    # error terms:
    b0_t_e = b0_t_e, b1_t_e = b1_t_e,
    c0_t_e = c0_t_e, c1_t_e = c1_t_e,
    d0_t_e = d0_t_e, d1_t_e = d1_t_e
  )
  
  
  list(
    NuisanceFit = NuisanceFit,
    TargetFit = TargetFit
  )
  
}


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
  #b0_pred=NuisanceFit$b0_pred # these are unused
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
  
  # solve(P) %*% Q = beta
  
  # Return them
  list(P = P, Q = Q)
}




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

Vab = function(NuisanceFit, TargetFit, gamma, beta_opt) {
  
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
  alpha_k = exp(cbind(X0 = 1, Xnew) %*% gamma)
  
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
  alpha_t = exp(cbind(X0 = 1, Xt) %*% gamma)
  
  dims <- ncol(TargetFit$r0_t)
  
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
  
  # conditional debiased - error terms:
  ab0_e_target <- (1/alpha_t)*b0_t_e
  ab0_e_source <- (1/alpha_k)*omega*((Y0new-m0_pred)^2-b0_pred_e)
  ab1_e_target <- alpha_t*b1_t_e
  ab1_e_source <- alpha_k*omega*((Y1new-m1_pred)^2-b1_pred_e)
  
  ac0_e_target <- as.vector(1/alpha_t)*c0_t_e
  ac0_e_source <- as.vector(1/alpha_k)*omega*((S0new-r0_pred)*as.vector(Y0new-m0_pred)-c0_pred_e)
  ac1_e_target <- as.vector(alpha_t)*c1_t_e
  ac1_e_source <- as.vector(alpha_k)*omega*((S1new-r1_pred)*as.vector(Y1new-m1_pred)-c1_pred_e)
  
  ad0_e_target <- as.vector(1/alpha_t)*d0_t_e
  ad0_e_source <- as.vector(1/alpha_k)*omega*(t(apply(S0new-r0_pred, 1, function(v) as.vector(outer(v,v)))) - d0_pred_e)
  ad1_e_target <- as.vector(alpha_t)*d1_t_e
  ad1_e_source <- as.vector(alpha_k)*omega*(t(apply(S1new-r1_pred, 1, function(v) as.vector(outer(v,v)))) - d1_pred_e)
  
  # V:
  V <- # uncontional:
    B1 - 2*C1 %*% beta_opt +  t(D1 %*% beta_opt) %*% beta_opt - M1^2 + 2*M1*t(R1) %*% beta_opt - (t(R1) %*% beta_opt)^2 +
    B0 - 2*C0 %*% beta_opt +  t(D0 %*% beta_opt) %*% beta_opt - M0^2 + 2*M0*t(R0) %*% beta_opt - (t(R0) %*% beta_opt)^2 +
    2*M1*M0 - 2*M1*t(R0) %*% beta_opt - 2*M0*t(R1) %*% beta_opt + 2*(R1 %*% beta_opt)*(R0 %*% beta_opt) +
    #------------------------------------------------- conditional -------------------------------------------------#
    (-2)*(
      mean(m1m0_target) + m1m0_source - 
        (colMeans(m1r0_target) + m1r0_source) %*% beta_opt -
        (colMeans(m0r1_target) + m0r1_source) %*% beta_opt +
        t(matrix(colMeans(r1r0_target) + r1r0_source,
                 nrow = dims, ncol = dims
        ) %*% beta_opt) %*% beta_opt
    ) +
    
    (# A=1
      # ab1_e
      mean(ab1_e_target) + mean(ab1_e_source[Anew == 1,,drop = FALSE]) +
        # -2ac1_e
        (-2)*(colMeans(ac1_e_target) + colMeans(ac1_e_source[Anew == 1,,drop = FALSE])) %*% beta_opt +
        # ad1_e
        t(matrix(colMeans(ad1_e_target) + colMeans(ad1_e_source[Anew == 1,,drop = FALSE]),
                 nrow = dims, ncol = dims
        ) %*% beta_opt) %*% beta_opt
    ) +
    (# A=0
      # ab0_e
      mean(ab0_e_target) + mean(ab0_e_source[Anew == 0,,drop = FALSE]) +
        # -2ac0_e
        (-2)*(colMeans(ac0_e_target) + colMeans(ac0_e_source[Anew == 0,,drop = FALSE])) %*% beta_opt +
        # ad0_e
        t(matrix(colMeans(ad0_e_target) + colMeans(ad0_e_source[Anew == 0,,drop = FALSE]),
                 nrow = dims, ncol = dims
        ) %*% beta_opt) %*% beta_opt 
    )
  
  V
  
}



Tau <- function(NuisanceFit, TargetFit, beta_opt) {
  
  # mean only considering meaningful values:
  Mean <- function(x) {
    mean(x[is.finite(x)])
  }
  
  # exclude negative sqrt and non-positive denominators:
  Sqrt <- function(x, y) {
    out <- numeric(length(x))              # default 0
    ok  <- x >= 0 & y > 0                  # where the ratio is valid
    out[ok] <- sqrt(x[ok] / y[ok])         # compute only there
    out
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
  
  var_e_target <- Sqrt(var_1_t, 1)*Sqrt(var_0_t, 1)
  var_e_source_0 <- omega*Sqrt(var_1,var_0)*((Y0new - S0new %*% beta_opt - m0_pred + r0_pred %*% beta_opt)^2 - var_0)
  var_e_source_1 <- omega*Sqrt(var_0,var_1)*((Y1new - S1new %*% beta_opt - m1_pred + r1_pred %*% beta_opt)^2 - var_1)
  
  var_target <- Sqrt(b1_t_e, 1)*Sqrt(b0_t_e, 1)
  var_source_0 <- omega*Sqrt(b1_pred_e,b0_pred_e)*((Y0new-m0_pred)^2 - b0_pred_e)
  var_source_1 <- omega*Sqrt(b0_pred_e,b1_pred_e)*((Y1new-m1_pred)^2 - b1_pred_e)
  
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



iteration_gamma <- function(NuisanceFit_folds, TargetFit_folds, beta_opt, gamma){
  
  hit_max_gamma <- FALSE # indicator variable, indicates whether max_iter_gamma is reached before stopping
  
  ha_folds    <- array(NA, dim = c(dimX+1, dimX+1, K))
  xi_folds    <- matrix(NA, nrow=dimX+1, ncol=K)
  
  for (iter in 1:max_iter_gamma) {
    
    for (j in 1:K){
      ha_folds[,,j] <- new_gamma(NuisanceFit_folds[[j]], TargetFit_folds[[j]], gamma, beta_opt)$ha
      xi_folds[,j]  <- new_gamma(NuisanceFit_folds[[j]], TargetFit_folds[[j]], gamma, beta_opt)$xi
    }
    avg_ha    <- apply(ha_folds, c(1, 2), mean, na.rm = TRUE)
    avg_xi    <- rowMeans(xi_folds, na.rm = TRUE)
    gamma_new <- gamma - step_size*(solve(avg_ha) %*% avg_xi)
    
    delta_gamma <- max(abs(gamma_new - gamma), na.rm = TRUE)
    
    # ░░░░░░░░░░░░░░░░░░░░░░░░░ display block: γ  ░░░░░░░░░░░░░░░░░░░░░░░░░░
    #cat(sprintf(
    #  "\niter %3d | |Δγ| = %.5f\n", iter, delta_gamma))
    #
    #cat("γ")
    #print(round(gamma_new, 6))                              # keeps 4×1 layout
    #
    #cat("----------------------------------------------------------------\n")
    # ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
    
    if (delta_gamma < tol) {
      gamma <- gamma_new
      #message("Converged!")
      break
    }
    gamma <- gamma_new
  }
  
  if (iter == max_iter_gamma) 
    hit_max_gamma <- TRUE
  
  list(gamma = gamma, hit_max_gamma = hit_max_gamma)
  
}



iteration_V <- function(NuisanceFit_folds, TargetFit_folds) {
  
  gamma_reached_max <- FALSE # indicator variable for hit_max_gamma in iteration_gamma()
  no_drop <- 3 # plateau length: stops if does not decrease for 3 consecutive steps
  
  # storage:
  best_V      <- Inf # stores the smallest possible so far
  best_gamma  <- gamma0  
  best_beta   <- matrix(NA, nrow = dimS, ncol = 1)
  
  P_hat_folds <- array(NA, dim = c(dimS, dimS, K))
  Q_hat_folds <- matrix(NA, nrow=dimS, ncol=K)
  V           <- numeric(K)
  tau_numr    <- numeric(K)
  tau_denom   <- numeric(K)
  
  # initiate:
  gamma <- gamma0
  prev_avg_V  <- Inf
  plateau_cnt <- 0
  
  for (iter in 1:max_iter_V) {
    
    # P, Q:
    for (j in 1:K) {
      P_hat_folds[,,j]  <- calculation_pq(NuisanceFit_folds[[j]], TargetFit_folds[[j]], gamma)$P # can also move out
      Q_hat_folds[,j]   <- calculation_pq(NuisanceFit_folds[[j]], TargetFit_folds[[j]], gamma)$Q
    }
    
    # Beta: (average P, Q over folds)
    avg_P    <- apply(P_hat_folds, c(1, 2), mean, na.rm = TRUE)
    avg_Q    <- rowMeans(Q_hat_folds, na.rm = TRUE)
    beta_opt <- solve(avg_P) %*% avg_Q
    #print(beta_opt) # debug
    
    # Gamma iteration:
    iteration_gamma_result = iteration_gamma(NuisanceFit_folds, TargetFit_folds, beta_opt, gamma)
    gamma_new <- iteration_gamma_result$gamma
    gamma_reached_max <- gamma_reached_max || iteration_gamma_result$hit_max_gamma # once any gamma iteration hits max, stored as TRUE
    
    # V, ha, xi:
    for (j in 1:K){
      V[j] <- Vab(NuisanceFit_folds[[j]], TargetFit_folds[[j]], gamma_new, beta_opt)
    }
    
    avg_V <- mean(V) # then evaluate the convergence of avg_V
    
    # ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
    #cat(sprintf(
    #  "\niter %3d | V = %.6f  ΔV = %.6f  plateau_cnt = %d\n",
    #  iter, avg_V, avg_V - prev_avg_V, plateau_cnt))
    #
    #cat("gamma =")
    #print(round(gamma_new, 6), quote = FALSE)     # 4 × 1 layout
    #
    #cat("beta  =")
    #print(round(beta_opt, 6), quote = FALSE)      # 2 × 1 layout
    #cat("------------------------------------------------------------\n")
    # ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
    
    # store best positive V so far
    if (avg_V > 0 && avg_V < best_V) {
      best_V     <- avg_V             
      best_gamma <- gamma_new         
      best_beta  <- beta_opt          
    }
    
    # stop immediately if V ≤ 0
    if (avg_V <= 0) {
      message("V became non-positive; rolling back to the previous best V.")
      break # once break, output the best_XXXs
    }
    
    # stop if non-decrease for 3 consecutive rounds:
    
    if (avg_V >= prev_avg_V - .Machine$double.eps) { # if not decreasing
      plateau_cnt <- plateau_cnt + 1 # considering only consecutive
    } else {
      plateau_cnt <- 0
    }
    
    delta_V <- abs(avg_V - prev_avg_V)
    
    if (plateau_cnt >= no_drop) {
      break # do NOT touch best_XXX → roll back to global minimum
    }
    
    if (delta_V < tol) {
      best_V     <- avg_V
      best_gamma <- gamma_new
      best_beta  <- beta_opt
      break
    }
    
    prev_avg_V <- avg_V
    gamma      <- gamma_new
    
  }
  
  beta_opt <- best_beta # used in tau calculation later
  
  hit_V_max <- as.integer(iter == max_iter_V) # V only does one set of iteration, so compare with max_iter_V once (differentiate from gamma)
  hit_gamma_max <- as.integer(gamma_reached_max) 
  
  for (j in 1:K){
    tau_numr[j] <- Tau(NuisanceFit_folds[[j]], TargetFit_folds[[j]], beta_opt)$numr
    tau_denom[j] <- Tau(NuisanceFit_folds[[j]], TargetFit_folds[[j]], beta_opt)$denom
  }
  
  tau <- 1-mean(tau_numr)/mean(tau_denom)
  return(list(
    tau = tau,
    beta_opt = beta_opt,
    gamma = gamma,
    hit_gamma_max = hit_gamma_max,
    hit_V_max     = hit_V_max
  ))
  
}


# ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░ Data set-up ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
# S0
ES0 <- function(X) {
  X1 <- X[, 1]
  X2 <- X[, 2]
  X3 <- X[, 3]
  
  ES0_1 = -2*sin(X1) + 0.5*X3
  ES0_2 = X2 + X1
  
  ES0 <- cbind(ES0_1, ES0_2)
  ES0
}

# S1
ES1 <- function(X) {
  X1 <- X[, 1]
  X2 <- X[, 2]
  X3 <- X[, 3]
  
  ES1_1 = sin(X1) + X3
  ES1_2 = -2*X2 - X1
  
  ES1 <- cbind(ES1_1, ES1_2)
  ES1
}

# Y0
EY0 <- function(X, S0) {
  beta_s0 <- c(1, -1)
  EY0 <- S0 %*% beta_s0 + z*(X[,1]*X[, 3] + 2*X[,1]^2)
  EY0
}

# Y1
EY1 <- function(X, S1) {
  beta_s1 <- c(1, -1)
  EY1 <- S1 %*% beta_s1 + z*(X[,1]*X[, 2] + X[,2]*X[, 3] + X[,1]^2)
  EY1
}

## error generator
#err <- function(n, d, var) { # n dp, d dimension, variance
#  err_matrix <- matrix(rnorm(n * d, mean = 0, sd = sqrt(var)), nrow = n, ncol = d)
#  return(err_matrix)
#}

# ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░ Data generate function ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
# GD - generating data function
DP <- list()
GD <- function(n, a, b, rho_y, rho_s1, rho_s2){
  
  l = 2 #?# use Try() to ensure data generates enough
  Xa <- data.frame(X1 = rnorm(l*2*n, 0, 0.25),
                   X2 = rnorm(l*2*n, 0, 1),
                   X3 = runif(l*2*n, min = -1, max = 1)) |> 
    mutate(p1 = expit(X1 + X2 + X3)
    ) # almost symmetric around 0
  Xa$t <- rbinom(n = length(Xa$p1), size = 1, prob = Xa$p1)
  
  X <- Xa |> filter(t == 0) |> slice(1:n) |> select(X1, X2, X3) |> as.matrix()
  Xt <- Xa |> filter(t == 1) |> slice(1:n) |> select(X1, X2, X3) |> as.matrix()
  
  # ░░░░░░░░░░░░░░░░░░░░░░░ error correlation ░░░░░░░░░░░░░░░░░░░░░░░░░
  
  ## Y-errors: 2×2 block
  sigma_y <- matrix(c(
    b, rho_y*sqrt(b*b),
    rho_y*sqrt(b*b), b),
    nrow = 2, byrow = TRUE)
  
  errs_Y  <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = sigma_y)
  e_Y0    <- errs_Y[, 1]
  e_Y1    <- errs_Y[, 2]
  
  ## S(·)_1 errors: 2×2 block
  sigma_s1 <- matrix(c(
    a, rho_s1*sqrt(a*a),
    rho_s1*sqrt(a*a),  a),
    nrow = 2, byrow = TRUE)
  
  errs_S1 <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = sigma_s1)
  e_S0_1  <- errs_S1[, 1]
  e_S1_1  <- errs_S1[, 2]
  
  ## S(·)_2 errors: 2×2 block
  sigma_s2 <- matrix(c(
    a, rho_s2*sqrt(a*a),
    rho_s2*sqrt(a*a),  a),
    nrow = 2, byrow = TRUE)
  
  errs_S2 <- MASS::mvrnorm(n, mu = c(0, 0), Sigma = sigma_s2)
  e_S0_2  <- errs_S2[, 1]
  e_S1_2  <- errs_S2[, 2]
  
  ## assemble matrices used later
  e_S0 <- cbind(e_S0_1, e_S0_2) 
  e_S1 <- cbind(e_S1_1, e_S1_2) 
  
  # uncorrelated error:
  #e_S0 <- err(n, 2, 0.5) # d=2, var=0.5
  #e_S1 <- err(n, 2, 0.5)
  #e_Y0 <- err(n, 1, 0.5) # 0.5
  #e_Y1 <- err(n, 1, 0.5) # 0.5
  
  # ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░
  #S: n x 3 matrix
  S0 <- ES0(X) + e_S0
  S1 <- ES1(X) + e_S1
  #Y: n x 1 matrix
  Y0 <- EY0(X, S0) + e_Y0
  Y1 <- EY1(X, S1) + e_Y1
  
  #A: n x 1 matrix
  A <- c(rep(0, n/2), rep(1, n/2)) # first half: A=0, second half: A=1 # for better sample split
  
  DP <- list('A'=A, 'X'=X, 'Xt'=Xt, 'Y0'=Y0,'Y1'=Y1, 'S0'=S0, 'S1'=S1)
  return(DP)
}


sim_function <- function(n, K, givenestimator, NBoot) {
  
  DP <- GD(n, a, b, rho_y, rho_s1, rho_s2)
  Index <- Split2(K, n) # sample split
  
  NuisanceFit_folds <- vector("list", K)
  TargetFit_folds   <- vector("list", K)
  # Bootstrap storage:
  NuisanceBoot_folds <- lapply(seq_len(NBoot), function(...) vector("list", K))
  TargetBoot_folds <- lapply(seq_len(NBoot), function(...) vector("list", K))
  tau_boot = numeric(NBoot)
  
  #------------------------------------- Nuisance Estimate -------------------------------------#
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
  
  #----------------------------------------- Iteration -----------------------------------------#
  
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




output <- list()

for (n in sample_sizes) {
  
  output[[as.character(n)]] <- data.frame()
  
  for (i in 1:N) {
    
    set.seed(200427 + i)
    
    max_iter_V   <- 10
    max_iter_gamma   <- 100
    tol        <- 0.00001
    step_size  <- 1 # set at = 0.1/0.2 for gradient update
    result <- sim_function(n, K, givenestimator, NBoot)
    
    # result$theta[[1]] is a vector of length 2
    tau           <- result$tau
    tau_se        <- result$se_tau
    beta          <- as.vector(result$beta)
    gamma         <- as.vector(result$gamma)
    hit_gamma_max <- result$hit_gamma_max
    hit_V_max     <- result$hit_V_max
    
    df <- data.frame(
      replicate = i,
      
      tau = tau,
      tau_se = tau_se,
      beta1 = beta[1],
      beta2 = beta[2],
      gamma1 = gamma[1],
      gamma2 = gamma[2],
      gamma3 = gamma[3],
      gamma4 = gamma[4],
      hit_gamma_max = hit_gamma_max,
      hit_V_max = hit_V_max
    )
    
    
    # append the result for this replication
    output[[as.character(n)]] <- rbind(output[[as.character(n)]], df)
    
    csv_filename <- paste0(
      givenestimator, "_",        
      n, "_",                    
      Sys.Date(),                
      ".csv"
    )
    
    write.csv(output[[as.character(n)]],
              file = file.path("./output", csv_filename), 
              row.names = FALSE)
    
    message(sprintf("[%s] iteration %d done — beginning iteration %d",
                    format(Sys.time(), "%m‑%d %H:%M"),  # timestamp to the minute
                    i, i + 1))
  }
}











