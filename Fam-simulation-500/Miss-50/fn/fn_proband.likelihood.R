proband.likelihood <- function(beta, gamma, xi, prelim.proband, M, K, range.t, allef) {
  n <- nrow(prelim.proband$obj)
  
  test0 <- rep(0, n)
  test1 <- rep(1, n) 
  d.1 <- prelim.proband$obj$D.1
  d.2 <- prelim.proband$obj$D.2
  tmp.vec0 <- cbind(test0, d.1, d.2)
  tmp.vec1 <- cbind(test1, d.1, d.2)
  
  pG <- allef[[1]][2]
  p0 <- 1- 0.005
  
  lambda <- (prelim.proband$ft %*% gamma) / range.t
  Lambda <- prelim.proband$Ft %*% gamma
  
  ratio0 <- tmp.vec0 %*% beta
  ratio1 <- tmp.vec1 %*% beta
  
  like0.1 <- (lambda[, 1] * xi[, 1] * exp(ratio0[, 1])) * exp(- apply(Lambda * xi * exp(ratio0), 1, sum))
  like1.1 <- (lambda[, 1] * xi[, 1] * exp(ratio1[, 1])) * exp(- apply(Lambda * xi * exp(ratio1), 1, sum))
  
  like0.2 <- (lambda[, 2] * xi[, 2] * exp(ratio0[, 2])) * exp(- apply(Lambda * xi * exp(ratio0), 1, sum))
  like1.2 <- (lambda[, 2] * xi[, 2] * exp(ratio1[, 2])) * exp(- apply(Lambda * xi * exp(ratio1), 1, sum))
  
  like1 <- p0 * like0.1 + (1 - p0) * like1.1
  like2 <- p0 * like0.2 + (1 - p0) * like1.2
  
  obj <- log(like1 + like2)
  
  obj
}


