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


result <- list()

for (n in sample_sizes) {
  
  result[[as.character(n)]] <- data.frame()
  
  for (i in 1:N) {
    set.seed(200427 + i)
    # N=1 inside sim_function so that each call returns one replication
    result <- sim_function(1, n, K, givenestimator, NBoot)
    
    # result$theta[[1]] is a vector of length 2
    tau       <- result$tau
    tau_se   <- result$se_tau
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
    result[[as.character(n)]] <- rbind(result[[as.character(n)]], df)
    
    csv_filename <- paste0(
      givenestimator, "_",        
      n, "_",                    
      Sys.Date(),                
      ".csv"
    )
    
    write.csv(result[[as.character(n)]], # save csv after each loop
              file = file.path("./output", csv_filename), 
              row.names = FALSE)
  }
}
