rm(list = ls())

library(matrixStats)
library(foreach)
library(doParallel)
registerDoParallel(cores = 8)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('fn/get_penetrance.R')
source('fn/get_true_penetrance.R')
source('get_lambda.R')

folders <- c("Fam-simulation-200-Miss-20", "Fam-simulation-200-Miss-50", 
             "Fam-simulation-300-Miss-20", "Fam-simulation-300-Miss-50",
             "Scenario-2", "Scenario-3", "Scenario-4", "Scenario-5", "Scenario-6",
             "Fam-simulation-500-Miss-20", "Fam-simulation-500-Miss-50")

true.beta <- c(4,2,3,3,3,2)
true.lambda <- c(0.001,0.001)
true.phi <- c(1,1)
true.param <- list(true.beta,true.lambda,true.phi)

true.penet.0 <- get_true_penetrance(0, true.param)
true.penet.1 <- get_true_penetrance(1, true.param)

bias <- array(0, dim = c(2,2,11,50))
CP <- array(0, dim = c(2,2,11,50))
MCIW <- array(0, dim = c(2,2,11,50))

for (k in 1:50) {
  for (i in 1:11) {
    print(c(k,i))
    
    beta <- read.csv(paste0(folders[i], '/post_beta_seed', k, '.csv'), header = TRUE)
    gamma <- read.csv(paste0(folders[i], '/post_gamma_seed', k, '.csv'), header = TRUE)
    phi <- read.csv(paste0(folders[i], '/post_phi_seed', k, '.csv'), header = TRUE)
    param <- list(beta, gamma, phi)
    
    penet.0 <- get_penetrance(0, c(0.025, 0.5, 0.975), param)
    penet.1 <- get_penetrance(1, c(0.025, 0.5, 0.975), param)
    
    bias[1,,i,k] <- apply(abs(penet.0[,,2] - true.penet.0), 2, mean)
    bias[2,,i,k] <- apply(abs(penet.1[,,2] - true.penet.1), 2, mean)
    CP[1,,i,k] <- apply((penet.0[,,1] < true.penet.0) & (penet.0[,,3] > true.penet.0), 2, mean)
    CP[2,,i,k] <- apply((penet.1[,,1] < true.penet.1) & (penet.1[,,3] > true.penet.1), 2, mean)
    MCIW[1,,i,k] <- apply(penet.0[,,3] - penet.0[,,1], 2, mean)
    MCIW[2,,i,k] <- apply(penet.1[,,3] - penet.1[,,1], 2, mean)
  }
}

bias.penet <- array(0, dim = c(2,2,11))
coverage.prob <- array(0, dim = c(2,2,11))
mean.CI.width <- array(0, dim = c(2,2,11))

for (i in 1:11) {
  bias.penet[1,,i] <- apply(bias[1,,i,], 1, mean)
  bias.penet[2,,i] <- apply(bias[2,,i,], 1, mean)
  
  coverage.prob[1,,i] <- apply(CP[1,,i,], 1, mean)
  coverage.prob[2,,i] <- apply(CP[2,,i,], 1, mean)
  
  mean.CI.width[1,,i] <- apply(MCIW[1,,i,], 1, mean)
  mean.CI.width[2,,i] <- apply(MCIW[2,,i,], 1, mean)
}

write.csv(cbind(bias.penet[1,,], bias.penet[2,,]), 
          file = 'summary/bias_first_penetrance.csv', row.names = FALSE)
write.csv(cbind(coverage.prob[1,,], coverage.prob[2,,]), 
          file = 'summary/CP_first_penetrance.csv', row.names = FALSE)
write.csv(cbind(mean.CI.width[1,,], mean.CI.width[2,,]), 
          file = 'summary/MCIW_first_penetrance.csv', row.names = FALSE)







