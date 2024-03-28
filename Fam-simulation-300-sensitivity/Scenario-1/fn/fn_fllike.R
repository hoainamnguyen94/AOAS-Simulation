fllike <- function(beta, gamma, xii, n1, n2, data1, data2, cancer, cancer_type,
                   Ft, ft, G, nG, M, K, range.t, allef, mRate)
{
  id <- data2$id
  time <- data2$time
  time[time == 0] <- 1.0e-12
  tilde.t <- data2$tilde.t
  gender <- data2$gender
  test <- data2$test
  dp.1 <- data2$dp.1
  dp.2 <- data2$dp.2
  
  ##### Key - likelihood part! ########
  #####################################
  test0 <- rep(0, n2)
  test1 <- rep(1, n2)
  
  xp.test0 <- cbind(test0, dp.1, dp.2)
  xp.test1 <- cbind(test1, dp.1, dp.2)
  m <- nrow(xp.test0)
  
  xpbeta.test0 <- xpbeta.test1 <- matrix(0, nrow = K, ncol = m)
  Lambda <- lambda <- matrix(0, nrow = K, ncol = m)
  unique.id <- unique(id)
  diff.Lambda <- matrix(0, nrow = K, ncol = m)
  log.l0 <- log.l1 <- matrix(0, nrow = K, ncol = m)
  log.L0 <- log.L1 <- matrix(0, nrow = K, ncol = m)
  
  for (k in 1:K){
    xpbeta.test0[k,] <- xp.test0 %*% beta[, k]
    xpbeta.test1[k,] <- xp.test1 %*% beta[, k]
    
    Lambda[k,] <- (Ft %*% gamma[, k])         # cumulative baseline
    lambda[k,] <- (ft %*% gamma[, k]) / range.t # baseline
    
    temp.diff.Lambda <- NULL
    for (ii in unique.id) {
      temp.diff.Lambda <- c(temp.diff.Lambda, diff(c(0, Lambda[k, id == ii])))
    }
    diff.Lambda[k,] <- temp.diff.Lambda
    
    log.l0[k,] <- log(lambda[k, ]) + log(xii[k]) + xpbeta.test0[k, ]
    log.l1[k,] <- log(lambda[k, ]) + log(xii[k]) + xpbeta.test1[k, ]
    
    log.L0[k,] <- diff.Lambda[k, ] * xii[k] * exp(xpbeta.test0[k, ])
    log.L1[k,] <- diff.Lambda[k, ] * xii[k] * exp(xpbeta.test1[k, ])
  }
  
  llike0 <- llike1 <- NULL
  for (ii in unique.id) {
    temp.id <- which(id == ii)
    L <- length(temp.id)
    
    if (L == 1){
      temp.llike0 <- -sum(log.L0[, temp.id])
      temp.llike1 <- -sum(log.L1[, temp.id])
    } else {
      temp.llike0 <- 0
      temp.llike1 <- 0
      for (l in 1:(L - 1)){
        ctype <- cancer[cancer$id == ii,]$cancer.type[l]
        temp.llike0 <- temp.llike0 + log.l0[ctype, temp.id[l]] - sum(log.L0[, temp.id[l]])
        temp.llike1 <- temp.llike1 + log.l1[ctype, temp.id[l]] - sum(log.L1[, temp.id[l]])
      }
      temp.llike0 <- temp.llike0 - sum(log.L0[, temp.id[L]])
      temp.llike1 <- temp.llike1 - sum(log.L1[, temp.id[L]])
    }
    llike0 <- c(llike0, temp.llike0)
    llike1 <- c(llike1, temp.llike1)
  }
  
  lik <- cbind(exp(llike0), exp(llike1), exp(llike1))
  
  i.test <- NULL
  for (ii in unique.id) {
    temp.id <- which(id == ii)
    i.test <- c(i.test, as.integer(mean(test[temp.id])))
  }
  id0 <- which(i.test == 0)
  id1 <- which(i.test == 1)
  
  lik[id0, 2:3] <- 0
  lik[id1, 1] <- 0
  
  obj <- CalPeelingProbLIK(allef, LIK = lik, ped = data1, counselee.id = 1, 1, mRate)
  obj
}
