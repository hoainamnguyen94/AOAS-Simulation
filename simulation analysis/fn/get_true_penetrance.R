get_true_penetrance <- function(g,param) {
  M <- 2
  K <- 2
  range.t <- 105
  
  covariate <- c(g,0,0)
  
  beta <- matrix(param[[1]], ncol = K)
  lambda0 <- matrix(param[[2]], ncol = K)
  phi <- param[[3]]
  
  lambda <- rep(1, range.t*10+1) %*% (lambda0/10)
  Lambda <- apply(lambda, 2, cumsum) 
  
  ratio <- covariate %*% beta
  
  S <- t(apply(Lambda, 1, function(x){phi / (phi + x * exp(ratio))}))
  penet <- apply(S * t(apply(lambda, 1, function(x){x * exp(ratio)})) * 
                   apply(S, 1, function(x){prod(x ^ phi)}), 2, cumsum)
  
  return(penet)
}
























