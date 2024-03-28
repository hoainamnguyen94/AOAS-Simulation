get.lambda <- function(post.gamma, range.t)
{
  time <- seq(0, range.t - range.t / 100, length = 101)
  tilde.time <- time / range.t
  # Bernstein Basis
  W <- unlist(lapply(1:M, function(k) pbeta(tilde.time, k, M - k + 1)))
  w <- unlist(lapply(1:M, function(k) dbeta(tilde.time, k, M - k + 1)))    
  
  Ft <- matrix(W, ncol = M)
  ft <- matrix(w, ncol = M)
  
  
  lambda <- apply(ft %*% t(post.gamma) / range.t, 1, median)
  Lambda <- apply(Ft %*% t(post.gamma), 1, median)
  data.frame(time = time, hazard = lambda, Hazard = Lambda)
}
