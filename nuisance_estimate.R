
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
  
  #### 1. glm:
  Xc$source <- as.factor(Xc$source)
  p <- glm(source ~ ., data = Xc, family = binomial) # logistic regression
  p1 <- predict(p, newdata = data.frame(Xnew), type = "response")
  omega <- (p1/(1-p1))*(1-1/K)
  
  fit_model <- function(outcome, predictor) {
    train(outcome ~ ., data = data.frame(outcome, predictor),
          method = 'glm',
          trControl = trainControl(method = "cv", number = 5, verboseIter = FALSE)
    )
  }
  
  ##### 2. machine learning:
  #Xc$source <- factor(Xc$source, 
  #                  levels = c(0,1), 
  #                  labels = c("Class0", "Class1"))
  #
  #fitControl_p <- trainControl( # for p
  #    method = "cv",
  #    number = 5,
  #    verboseIter=FALSE,
  #    classProbs  = TRUE) # for classification
  #
  #fitControl <- trainControl( # for nuisance models
  #    method = "cv",
  #    number = 5,
  #    verboseIter=FALSE)
  #
  ## nuisance:
  #fit_model <- function(outcome, predictor) {
  #  train(outcome ~ ., data = data.frame(outcome, predictor),
  #        method = givenestimator,
  #        trControl = fitControl,
  #        verbose = FALSE
  #  )}
  #
  ## p & omega, for different algrthms:
  #if(givenestimator %in% c('svmLinear', 'svmLinear2')){
  #  
  #  p <- train(source ~., data = Xc, 
  #             method = givenestimator, 
  #             trControl =fitControl_p,
  #             verbose = FALSE,
  #             probability = TRUE) # unique to svmLinear2
  #  
  #}else{
  #  
  #  p <- train(source ~., data = Xc, 
  #             method = givenestimator, 
  #             trControl =fitControl_p,
  #             verbose = FALSE) # c('gbm','rf','nnet')
  #  
  #}
  #
  #p1 <- predict(p, newdata = data.frame(Xnew), type = "prob")$Class1
  #omega <- (p1/(1-p1))*(1-1/K)
  
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
    Xnew=Xnew,
    Anew=Anew,
    omega=omega,
    Y0new=Y0new, Y1new=Y1new, S0new=S0new, S1new=S1new, 
    m0_pred=m0_pred, m1_pred=m1_pred,
    r0_pred=r0_pred, r1_pred=r1_pred,
    b0_pred=b0_pred, b1_pred=b1_pred,
    c0_pred=c0_pred, c1_pred=c1_pred,
    d0_pred=d0_pred, d1_pred=d1_pred,
    # error terms:
    b0_pred_e=b0_pred_e,b1_pred_e=b1_pred_e,
    c0_pred_e=c0_pred_e,c1_pred_e=c1_pred_e,
    d0_pred_e=d0_pred_e,d1_pred_e=d1_pred_e
  )
  
  TargetFit <- list(
    Xt = Xt,
    m0_t=m0_t, m1_t=m1_t, 
    r0_t=r0_t, r1_t=r1_t,
    b0_t=b0_t, b1_t=b1_t,
    c0_t=c0_t, c1_t=c1_t,
    d0_t=d0_t, d1_t=d1_t,
    # error terms:
    b0_t_e=b0_t_e, b1_t_e=b1_t_e,
    c0_t_e=c0_t_e, c1_t_e=c1_t_e,
    d0_t_e=d0_t_e, d1_t_e=d1_t_e
  )

  
  list(
    NuisanceFit = NuisanceFit,
    TargetFit = TargetFit
  )
  
}
