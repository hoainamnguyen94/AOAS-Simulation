get_true_second_penetrance <- function(g,t1,fp,param) {
  M <- 2
  K <- 2
  range.t <- 105
  
  if (fp == 1) {
    covariate <- c(g, 1, 0)
  } else {
    covariate <- c(g, 0, 1)
  }
  
  beta <- matrix(param[[1]], ncol = K)
  lambda0 <- matrix(param[[2]], ncol = K)
  phi <- param[[3]]
  
  lambda <- rep(1, range.t*10+1) %*% (lambda0/10)
  Lambda <- apply(lambda, 2, cumsum)
  
  ratio <- covariate %*% beta
  init <- Lambda[t1*10,]
  lambda <- lambda[(t1*10+1):(range.t*10+1),]
  Lambda <- Lambda[(t1*10+1):(range.t*10+1),]
  
  S <- t(apply(Lambda, 1, function(x){phi / (phi + (x - init) * exp(ratio))}))
  penet <- apply(S * t(apply(lambda, 1, function(x){x * exp(ratio)})) * 
                   apply(S, 1, function(x){prod(x ^ phi)}), 2, cumsum)
  
  return(penet)
}














