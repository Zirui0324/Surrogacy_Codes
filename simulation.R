library(tidyverse)
library(rje)
library(caret)
library(ks)
library(kernlab)

N=1000
n=1000
sample_sizes = c(n)
K=5
givenestimator = 'glm' # c('glm', 'svmLinear2', 'rf', 'nnet)
NBoot=100
dimX = 3
gamma0 = matrix(0, dimX+1, 1)

# -------------------------------------------- Simulation -------------------------------------------- #
output <- list()

for (n in sample_sizes) {
  
  output[[as.character(n)]] <- data.frame()
  
  for (i in 1:N) {
    
    set.seed(200427 + i)
    tol        <- 0.0001
    max_iter   <- 100
    step_size  <- 1 # set at = 0.1/0.2 for gradient update
    result <- sim_function(n, K, givenestimator, NBoot)
    
    # result$theta[[1]] is a vector of length 2
    tau       <- result$tau
    tau_se    <- result$se_tau
    beta      <- as.vector(result$beta)
    gamma     <- as.vector(result$gamma)
    
    df <- data.frame(
      replicate = i,
      
      tau = tau,
      tau_se = tau_se,
      beta1 = beta[1],
      beta2 = beta[2],
      gamma1 = gamma[1],
      gamma2 = gamma[2],
      gamma3 = gamma[3],
      gamma4 = gamma[4]
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

# -------------------------------------------- Output Organization -------------------------------------------- #
output_raw = read_csv("./output/glm_1000_2025-06-05.csv")







