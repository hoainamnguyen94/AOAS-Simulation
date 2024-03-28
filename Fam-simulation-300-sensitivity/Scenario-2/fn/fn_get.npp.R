get.npp <- function(pi, beta, lambda, phi, theta.c, xi = NULL, G = NULL, S = NULL) {
  # censored time
  cc <- runif(1) * theta.c
  
  # frailty
  if (is.null(xi)) {
    xi <- c(rgamma(1, phi[1], phi[1]), rgamma(1, phi[2], phi[2]))
  } 
  
  # genotype
  if (is.null(G)) {
    G <- rbinom(1, 1, 0.005)   
  } 
  
  # gender
  if (is.null(S)) {
    S <- 2 - rbinom(1, 1, 1/2)
  } else if (S == 0) {
    S <- 2 
  } else if (S == 1) {
    S <- 1
  } else if (S == 2) {
    S <- 2
  } else stop("gender should be coded as 1 for male and 2 for female")
  
  # number obs to be sampled
  nn <- 100
  
  x1 <- cbind(G, 0, 0)
  x2 <- cbind(G, 1, 0)
  x3 <- cbind(G, 0, 1)
  
  h1.1 <- lambda[1] * xi[1] * exp(x1 %*% beta[, 1])
  h2.1 <- lambda[2] * xi[2] * exp(x1 %*% beta[, 2])
  
  w1.1 <- rexp(1, h1.1)
  w1.2 <- rexp(1, h2.1)
  w1 <- min(w1.1, w1.2)
  c.type1 <- ifelse(w1.2 > w1.1, 1, 2)
  
  if (w1.2 > w1.1) {
    h1.2 <- lambda[1] * xi[1] * exp(x2 %*% beta[, 1])
    h2.2 <- lambda[2] * xi[2] * exp(x2 %*% beta[, 2])
    
    w2.1 <- rexp(nn - 1, h1.2)
    w2.2 <- rexp(nn - 1, h2.2)
    w2 <- apply(cbind(w2.1, w2.2), MARGIN = 1, FUN = min)
    c.type2 <- apply(cbind(w2.1, w2.2), MARGIN = 1, 
                     FUN = function(x){ifelse(x[2] > x[1], 1, 2)})
  } else {
    h1.3 <- lambda[1] * xi[1] * exp(x3 %*% beta[, 1])
    h2.3 <- lambda[2] * xi[2] * exp(x3 %*% beta[, 2])
    
    w2.1 <- rexp(nn - 1, h1.3)
    w2.2 <- rexp(nn - 1, h2.3)
    w2 <- apply(cbind(w2.1, w2.2), MARGIN = 1, FUN = min)
    c.type2 <- apply(cbind(w2.1, w2.2), MARGIN = 1, 
                     FUN = function(x){ifelse(x[2] > x[1], 1, 2)})
    
  }
  
  w <- c(w1, w2)
  tt <- cumsum(w)
  c.type <- c(c.type1, c.type2)
  
  y <- c(tt[tt < cc], cc)
  n <- length(y)
  
  if (n == 1) {
    x <- x1 
  } else {
    if (w1.2 > w1.1) {
      x <- rbind(x1, matrix(rep(x2, n - 1), nrow = n - 1, byrow = T))
    } else {
      x <- rbind(x1, matrix(rep(x3, n - 1), nrow = n - 1, byrow = T))
    }
  }
  
  obj <- cbind(y, x[1:n,,drop = F], rep(xi[1], n), rep(xi[2], n), c.type[1:n])
  colnames(obj) <- c("time", "g", "dp.1", "dp.2", "xi.1", "xi.2", "ct")
  
  if (nrow(obj) > 10) {
    obj <- obj[1:11,]
    obj[11,1] <- min(obj[10,1] + 0.5, theta.c)
  }
  
  return(obj)
}