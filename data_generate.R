# S0
ES0 <- function(X) {
  X1 <- X[, 1]
  X2 <- X[, 2]
  X3 <- X[, 3]
  
  # ES0_1 = -2*sin(X1) + 0.5*X3
  # ES0_2 = X2 + X1
  ES0_1 = -2*X1
  ES0_2 = X2
  
  
  ES0 <- cbind(ES0_1, ES0_2)
  ES0
}

# S1
ES1 <- function(X) {
  X1 <- X[, 1]
  X2 <- X[, 2]
  X3 <- X[, 3]
  
  #ES1_1 = sin(X1) + X3
  #ES1_2 = -2*X2 - X1
  ES1_1 = X1
  ES1_2 = -2*X2
  
  ES1 <- cbind(ES1_1, ES1_2)
  ES1
}

# Y0
EY0 <- function(X, S0) {
  #beta_s0 <- c(1, -2)
  beta_s0 <- c(1, -1)
  EY0 <- S0 %*% beta_s0 # + X[,1]*X[, 3] + 2*X[,1]^2
  EY0
}

# Y1
EY1 <- function(X, S1) {
  #beta_s1 <- c(0.5, -0.5)
  beta_s1 <- c(1, -1)
  EY1 <- S1 %*% beta_s1 #+ 0.2*X[, 3] # + X[,2]*X[, 3] + X[,1]*X[, 2] + X[,1]^2 + 5*X[, 3]
  EY1
}

# error generator
err <- function(n, d, var) { # n dp, d dimension, variance
  err_matrix <- matrix(rnorm(n * d, mean = 0, sd = sqrt(var)), nrow = n, ncol = d)
  return(err_matrix)
}

# GD - generating data function
DP <- list()
GD <- function(n){
  
  l = 2 #?# use Try() to ensure data generates enough
  Xa <- data.frame(X1 = rnorm(l*2*n, 0, 0.5),
                   X2 = rnorm(l*2*n, 0, 1),
                   X3 = runif(l*2*n, min = -1, max = 1)) |> 
    mutate(p1 = rbinom(l*2*n, 1, 0.5)
           #p1 = expit(X1 + X2 + X3)
           ) # almost symmetric around 0
  Xa$t <- rbinom(n = length(Xa$p1), size = 1, prob = Xa$p1)
  
  X <- Xa |> filter(t == 0) |> slice(1:n) |> select(X1, X2, X3) |> as.matrix()
  Xt <- Xa |> filter(t == 1) |> slice(1:n) |> select(X1, X2, X3) |> as.matrix()
  
  
  #S: n x 3 matrix
  e_S0 <- err(n, 2, 0.5) # d=2, var=0.5
  e_S1 <- err(n, 2, 0.5)
  S0 <- ES0(X) + e_S0
  S1 <- ES1(X) + e_S1
  
  #Y: n x 1 matrix
  e_Y0 <- err(n, 1, 0) # 0.5
  e_Y1 <- err(n, 1, 0) # 0.5
  Y0 <- EY0(X, S0) + e_Y0
  Y1 <- EY1(X, S1) + e_Y1
  
  #A: n x 1 matrix
  A <- c(rep(0, n/2), rep(1, n/2)) # first half: A=0, second half: A=1 # for better sample split
  
  
  # X <- as.data.frame(X)
  # Xt <- as.data.frame(Xt)
  DP <- list('A'=A, 'X'=X, 'Xt'=Xt, 'Y0'=Y0,'Y1'=Y1, 'S0'=S0, 'S1'=S1)
  return(DP)
}

