get_penetrance <- function(g,q,param) {
  M <- 2
  K <- 2
  range.t <- 105
  grid <- seq(0, range.t, length = range.t*10+1)
  
  W <- unlist(lapply(1:M, function(k) pbeta(grid / range.t, k, M - k + 1)))
  w <- unlist(lapply(1:M, function(k) dbeta(grid / range.t, k, M - k + 1)))
  
  Ft <- matrix(W, ncol = M)
  ft <- matrix(w, ncol = M)
  
  covariate <- c(g,0,0)
  
  obj <- foreach(i = 1:15000) %dopar% {
    beta <- matrix(as.numeric(param[[1]][i + 15000,]), ncol = K)
    gamma <- matrix(as.numeric(param[[2]][i + 15000,]), ncol = K)
    phi <- as.numeric(param[[3]][i + 15000,])
    
    Lambda <-  Ft %*% gamma
    lambda <- (ft %*% gamma) / (range.t*10+1)
    
    ratio <- covariate %*% beta
    
    S <- t(apply(Lambda, 1, function(x){phi / (phi + x * exp(ratio))}))
    penet <- apply(S * t(apply(lambda, 1, function(x){x * exp(ratio)})) * 
                     apply(S, 1, function(x){prod(x ^ phi)}), 2, cumsum)
    
    penet
  }
  
  post.cpe <- array(unlist(obj), dim = c(range.t*10+1, K, 15000))
  cpe <- array(0, dim = c(range.t*10+1, K, length(q)))
  
  for (j in 1:K) {
    cpe[,j,] <- rowQuantiles(post.cpe[,j,], probs = q)
  }
  
  return(cpe)
}
























