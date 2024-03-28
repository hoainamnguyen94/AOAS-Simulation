rm(list = ls())

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source('get_lambda.R')

seed.list <- 1:50
nn <- length(seed.list)

folders <- c("Fam-simulation-200-Miss-20", "Fam-simulation-200-Miss-50", 
             "Fam-simulation-300-Miss-20", "Fam-simulation-300-Miss-50",
             "Scenario-2", "Scenario-3", "Scenario-4", "Scenario-5", "Scenario-6",
             "Fam-simulation-500-Miss-20", "Fam-simulation-500-Miss-50")

#######beta#######
est.beta <- array(0, dim = c(11,6,nn))

for (k in 1:nn){
  seed <- seed.list[k]
  print(seed)
  
  for (i in 1:11) {
    beta.file <- paste0(folders[i], '/post_beta_seed', seed, '.csv')
    beta <- read.csv(beta.file, header = TRUE)
    est.beta[i,,k]  <- apply(beta[-(1:15000),], 2, median)
  }
}

bias.beta <- matrix(0, nrow = 11, ncol = 6)
rmse.beta <- matrix(0, nrow = 11, ncol = 6)

true.beta <- c(4, 2, 3, 3, 3, 2)

for (i in 1:6) {
  bias.beta[,i] <- apply(abs(est.beta[,i,] - true.beta[i]), 1, mean)
  rmse.beta[,i] <- sqrt(apply((est.beta[,i,] - true.beta[i])^2, 1, mean))
}

write.csv(bias.beta, file = 'summary/bias_beta.csv', row.names = FALSE)
write.csv(rmse.beta, file = 'summary/rmse_beta.csv', row.names = FALSE)

#######phi#######
est.phi <- array(0, dim = c(11,2,nn))

for (k in 1:nn){
  seed <- seed.list[k]
  print(seed)
  
  for (i in 1:11) {
    phi.file <- paste0(folders[i], '/post_phi_seed', seed, '.csv')
    phi <- read.csv(phi.file, header = TRUE)
    est.phi[i,,k]  <- apply(phi[-(1:15000),], 2, median)
  }
}

bias.phi <- matrix(0, nrow = 11, ncol = 2)
rmse.phi <- matrix(0, nrow = 11, ncol = 2)

true.phi <- c(1,1)

for (i in 1:2) {
  bias.phi[,i] <- apply(abs(est.phi[,i,] - true.phi[i]), 1, mean)
  rmse.phi[,i] <- sqrt(apply((est.phi[,i,] - true.phi[i])^2, 1, mean))
}

write.csv(bias.phi, file = 'summary/bias_phi.csv', row.names = FALSE)
write.csv(rmse.phi, file = 'summary/rmse_phi.csv', row.names = FALSE)

#######lambda#######
est.lambda1 <- est.lambda2 <- array(0, dim = c(11, 105, nn))

for (k in 1:nn){
  seed <- seed.list[k]
  print(seed)
  
  for (i in 1:11) {
    gamma.file <- paste0(folders[i], '/post_gamma_seed', seed, '.csv')
    gamma <- read.csv(gamma.file, header = TRUE)
    lambda1 <- get.lambda(gamma[,1:2], 105, 2)
    lambda2 <- get.lambda(gamma[,3:4], 105, 2)
    
    est.lambda1[i,,k] <- lambda1[-1,2]
    est.lambda2[i,,k] <- lambda2[-1,2]
  }
}

bias.lambda1 <- bias.lambda2 <- NULL
rmse.lambda1 <- rmse.lambda2 <-  NULL

for (j in 1:11){
  bias.lambda1[j] <- mean(abs(est.lambda1[j,,] - 0.001))
  bias.lambda2[j] <- mean(abs(est.lambda2[j,,] - 0.001))
  rmse.lambda1[j] <- sqrt(mean((est.lambda1[j,,] - 0.001)^2))
  rmse.lambda2[j] <- sqrt(mean((est.lambda2[j,,] - 0.001)^2))
}

write.csv(rbind(bias.lambda1, rmse.lambda1), file = 'summary/bias_rmse_lambda1.csv', row.names = FALSE)
write.csv(rbind(bias.lambda2, rmse.lambda2), file = 'summary/bias_rmse_lambda2.csv', row.names = FALSE)
