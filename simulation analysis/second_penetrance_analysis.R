rm(list = ls())

library(matrixStats)
library(foreach)
library(doParallel)
registerDoParallel(cores = 8)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('fn/get_second_penetrance.R')
source('fn/get_true_second_penetrance.R')
source('get_lambda.R')

folders <- c("Fam-simulation-200-Miss-20", "Fam-simulation-200-Miss-50", 
             "Fam-simulation-300-Miss-20", "Fam-simulation-300-Miss-50",
             "Scenario-2", "Scenario-3", "Scenario-4", "Scenario-5", "Scenario-6",
             "Fam-simulation-500-Miss-20", "Fam-simulation-500-Miss-50")

true.beta <- c(4,2,3,3,3,2)
true.lambda <- c(0.001,0.001)
true.phi <- c(1,1)
true.param <- list(true.beta,true.lambda,true.phi)

true.penet.01 <- get_true_second_penetrance(g = 0, t1 = 30, fp = 1, true.param)
true.penet.11 <- get_true_second_penetrance(g = 1, t1 = 30, fp = 1, true.param)
true.penet.02 <- get_true_second_penetrance(g = 0, t1 = 30, fp = 2, true.param)
true.penet.12 <- get_true_second_penetrance(g = 1, t1 = 30, fp = 2, true.param)

bias <- array(0, dim = c(4,2,11,50))
CP <- array(0, dim = c(4,2,11,50))
MCIW <- array(0, dim = c(4,2,11,50))

for (k in 1:50) {
  for (i in 1:11) {
    print(c(k,i))
    
    beta <- read.csv(paste0(folders[i], '/post_beta_seed', k, '.csv'), header = TRUE)
    gamma <- read.csv(paste0(folders[i], '/post_gamma_seed', k, '.csv'), header = TRUE)
    phi <- read.csv(paste0(folders[i], '/post_phi_seed', k, '.csv'), header = TRUE)
    param <- list(beta, gamma, phi)
    
    penet.01 <- get_second_penetrance(g = 0, t1 = 30, fp = 1, c(0.025, 0.5, 0.975), param)
    penet.11 <- get_second_penetrance(g = 1, t1 = 30, fp = 1, c(0.025, 0.5, 0.975), param)
    penet.02 <- get_second_penetrance(g = 0, t1 = 30, fp = 2, c(0.025, 0.5, 0.975), param)
    penet.12 <- get_second_penetrance(g = 1, t1 = 30, fp = 2, c(0.025, 0.5, 0.975), param)
    
    bias[1,,i,k] <- apply(abs(penet.01[,,2] - true.penet.01), 2, mean)
    bias[2,,i,k] <- apply(abs(penet.11[,,2] - true.penet.11), 2, mean)
    bias[3,,i,k] <- apply(abs(penet.02[,,2] - true.penet.02), 2, mean)
    bias[4,,i,k] <- apply(abs(penet.12[,,2] - true.penet.12), 2, mean)
    
    CP[1,,i,k] <- apply((penet.01[,,1] < true.penet.01) & (penet.01[,,3] > true.penet.01), 2, mean)
    CP[2,,i,k] <- apply((penet.11[,,1] < true.penet.11) & (penet.11[,,3] > true.penet.11), 2, mean)
    CP[3,,i,k] <- apply((penet.02[,,1] < true.penet.02) & (penet.02[,,3] > true.penet.02), 2, mean)
    CP[4,,i,k] <- apply((penet.12[,,1] < true.penet.12) & (penet.12[,,3] > true.penet.12), 2, mean)
    
    MCIW[1,,i,k] <- apply(penet.01[,,3] - penet.01[,,1], 2, mean)
    MCIW[2,,i,k] <- apply(penet.11[,,3] - penet.11[,,1], 2, mean)
    MCIW[3,,i,k] <- apply(penet.02[,,3] - penet.02[,,1], 2, mean)
    MCIW[4,,i,k] <- apply(penet.12[,,3] - penet.12[,,1], 2, mean)
  }
}

bias.penet <- array(0, dim = c(4,2,11))
coverage.prob <- array(0, dim = c(4,2,11))
mean.CI.width <- array(0, dim = c(4,2,11))

for (i in 1:11) {
  bias.penet[1,,i] <- apply(bias[1,,i,], 1, mean)
  bias.penet[2,,i] <- apply(bias[2,,i,], 1, mean)
  bias.penet[3,,i] <- apply(bias[3,,i,], 1, mean)
  bias.penet[4,,i] <- apply(bias[4,,i,], 1, mean)
  
  coverage.prob[1,,i] <- apply(CP[1,,i,], 1, mean)
  coverage.prob[2,,i] <- apply(CP[2,,i,], 1, mean)
  coverage.prob[3,,i] <- apply(CP[3,,i,], 1, mean)
  coverage.prob[4,,i] <- apply(CP[4,,i,], 1, mean)
  
  mean.CI.width[1,,i] <- apply(MCIW[1,,i,], 1, mean)
  mean.CI.width[2,,i] <- apply(MCIW[2,,i,], 1, mean)
  mean.CI.width[3,,i] <- apply(MCIW[3,,i,], 1, mean)
  mean.CI.width[4,,i] <- apply(MCIW[4,,i,], 1, mean)
}

write.csv(rbind(bias.penet[1,,], bias.penet[2,,], bias.penet[3,,], bias.penet[4,,]), 
          file = 'summary/bias_second_penetrance.csv', row.names = FALSE)
write.csv(rbind(coverage.prob[1,,], coverage.prob[2,,], coverage.prob[3,,], coverage.prob[4,,]), 
          file = 'summary/CP_second_penetrance.csv', row.names = FALSE)
write.csv(rbind(mean.CI.width[1,,], mean.CI.width[2,,], mean.CI.width[3,,], mean.CI.width[4,,]), 
          file = 'summary/MCIW_second_penetrance.csv', row.names = FALSE)







