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
a = 0.5 # var S
b = 0.5 # var Y
rho_y  <- 0       # Corr(Y0 , Y1)
rho_s1 <- 0       # Corr(S0_1 , S1_1)
rho_s2 <- 0       # Corr(S0_2 , S1_2)

# -------------------------------------------- Simulation -------------------------------------------- #
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


# -------------------------------------------- Output Organization -------------------------------------------- #

output_wrap <- function(df, 
                        tau_conser, beta1_conser, beta2_conser,
                        tau_true,  beta1_true,  beta2_true
                        ) {
  
  #-------------------------------- FOR TAU: --------------------------------
  # bias:
  tau_bias_true   <- mean(df$tau, na.rm = TRUE) - tau_true
  tau_bias_conser <- mean(df$tau, na.rm = TRUE) - tau_conser
  # se:
  tau_ese  <- sd(df$tau, na.rm = TRUE)
  tau_ase  <- mean(df$tau_se, na.rm = TRUE)
  # coverage prob:
  lower_tau_ese <- df$tau - qnorm(0.95) * tau_ese
  lower_tau_ase <- df$ase_lower
  
  tau_cp_conser_ese <- mean(lower_tau_ese <= tau_conser, na.rm = TRUE)
  tau_cp_conser_ase <- mean(lower_tau_ase <= tau_conser, na.rm = TRUE)
  
  #-------------------------------- FOR BETA: --------------------------------
  # beta1:
  # bias:
  beta1_bias_true   <- mean(df$beta1, na.rm = TRUE) - beta1_true
  beta1_bias_conser <- mean(df$beta1, na.rm = TRUE) - beta1_conser
  # se:
  beta1_ese  <- sd(df$beta1, na.rm = TRUE)
  # coverage prob:
  lower_beta1_ese <- df$beta1 - qnorm(0.95) * beta1_ese
  beta1_cp_conser_ese <- mean(lower_beta1_ese <= beta1_conser, na.rm = TRUE)
  
  # beta2:
  # bias:
  beta2_bias_true   <- mean(df$beta2, na.rm = TRUE) - beta2_true
  beta2_bias_conser <- mean(df$beta2, na.rm = TRUE) - beta2_conser
  # se:
  beta2_ese  <- sd(df$beta2, na.rm = TRUE)
  # coverage prob:
  lower_beta2_ese <- df$beta2 - qnorm(0.95) * beta2_ese
  beta2_cp_conser_ese <- mean(lower_beta2_ese <= beta2_conser, na.rm = TRUE)

  
  # output:
  result <- rbind(
    data.frame(parameter    = "tau",
               truth_conser = tau_conser,
               bias_conser  = tau_bias_conser,
               se_type      = "ese",
               se           = tau_ese,
               cp_conser    = tau_cp_conser_ese,
               truth         = tau_true,
               bias_truth    = tau_bias_true,
               row.names = NULL),
    data.frame(parameter    = "tau",
               truth_conser = NA,
               bias_conser  = NA,
               se_type      = "ase",
               se           = tau_ase,
               cp_conser    = tau_cp_conser_ase,
               truth         = NA,
               bias_truth    = NA,
               row.names = NULL),
    data.frame(parameter    = "beta1",
               truth_conser = beta1_conser,
               bias_conser  = beta1_bias_conser,
               se_type      = "ese",
               se           = beta1_ese,
               cp_conser    = beta1_cp_conser_ese,
               truth         = beta1_true,
               bias_truth    = beta1_bias_true,
               row.names = NULL),
    data.frame(parameter    = "beta2",
               truth_conser = beta2_conser,
               bias_conser  = beta2_bias_conser,
               se_type      = "ese",
               se           = beta2_ese,
               cp_conser    = beta2_cp_conser_ese,
               truth         = beta2_true,
               bias_truth    = beta2_bias_true,
               row.names = NULL)
  ) |> 
    mutate(across(where(is.numeric), ~ round(.x, 4)))
  
  return(result)
}

# ---------------------------------------------------------------------------------------- #

output_lap_S5_2 = 
  as.data.frame(read_csv("./output/lap_S5_2.csv")) |> 
  mutate(ase_lower = tau - qnorm(0.95)*tau_se) |> 
  select(replicate, tau, tau_se, ase_lower, everything())

output_som_S5_2 = 
  as.data.frame(read_csv("./output/som_output/som_S5_2.csv")) |> 
  mutate(ase_lower = tau - qnorm(0.95)*tau_se) |> 
  select(replicate, tau, tau_se, ase_lower, everything())

output_som_S9_z_0.2 = 
  as.data.frame(read_csv("./output/som_output/som_S9_z_0.2_gbm.csv")) |> 
  mutate(ase_lower = tau - qnorm(0.95)*tau_se) |> 
  select(replicate, tau, tau_se, ase_lower, everything())

# S9:
tau_conser = 0.787788599
beta1_conser = 1.02566948
beta2_conser = -1.0185209
tau_true = 0.9007807
beta1_true = 1.029281
beta2_true = -1.022899

tau_conser = 
beta1_conser = 
beta2_conser = 
tau_true = 
beta1_true = 
beta2_true = 
  
# 6:
tau_conser = 0.82459251
beta1_conser = 1.01606806
beta2_conser = -1.0145708
tau_true = 0.9075922
beta1_true = 1.033031
beta2_true = -1.015984 

# 4:
tau_conser = 0.83173506
beta1_conser = 1.01072414
beta2_conser = -1.0096003
tau_true = 0.9140596
beta1_true = 1.024362 
beta2_true = -1.011264

# 2:

tau_conser = 0.83400753
beta1_conser = 1.00401278
beta2_conser = -0.9999239
tau_true = 0.9148054
beta1_true = 1.004012 
beta2_true = -1.004039

analysis_som_S6_V1_gbm <- output_wrap(output_som_S6_V1_gbm, tau_conser, tau_true)
analysis_som_S6_V2_gbm <- output_wrap(output_som_S6_V2_gbm |> filter(tau_se<=1), tau_conser, tau_true)

analysis2_som_S6_V1_gbm <- output_wrap(output_som_S6_V1_gbm
                                       |> filter(hit_gamma_max != 1, hit_V_max != 1), tau_conser, tau_true)

a6 = output_wrap(output_som_S9_z_0.6[-14, ], 
            tau_conser, beta1_conser, beta2_conser,
            tau_true,  beta1_true,  beta2_true)

