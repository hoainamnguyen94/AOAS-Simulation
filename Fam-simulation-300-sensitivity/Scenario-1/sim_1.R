library(foreach)
library(doParallel)
registerDoParallel(cores = 25)
library(LFSPRO)

sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

sourceDir("fn")

batch <- 1
n.sim <- 25
start <- n.sim * (batch - 1) + 1
end <- n.sim * batch

K <- 2
N <- 300
phi <- rep(1, K)

# estimation setting
n.sample <- 30000  # of posterior samples
M <- 2             # degrees of Bernstein polynomials

# Set-up
true.beta <- beta <- matrix(c(4, 2, 3, 3, 3, 2), nrow = 3) # regression coefficients (test, D_1, D_2)
lambda <- c(0.001, 0.001)     # baseline intensity
theta.c <- 80                 # hazard for censoring time
pG <- pi <- 0.001             # maf prevalence

cancer_type <- c(1, 2)
  
OBJ <- foreach(iter = start:end) %dopar% {
  dyn.load("fn/cpp/source/CalPeelingProbLIK.so")
  
  temp <- data.gen(seed = iter, N, pi, beta, lambda, theta.c, phi, miss = T)
    
  data.obj1 <- temp$data.obj1 # pedigree
  data.obj2 <- temp$data.obj2 # disease history
  cancer <- temp$cancer.obj # cancer
    
  # rescaling constant for Bernstein Polynomials
  range.t <- max(unlist(lapply(data.obj2, function(x) max(x[,2])))) + 1.0e-3
    
  # Initial Values
  init  <- list(beta = true.beta - 5, gamma = matrix(1, nrow = M, ncol = K), 
                xi = matrix(0.1, nrow = N, ncol = K), phi = phi)
  delta <- list(beta = matrix(c(0.7, 1, 0.8, 0.7, 1, 0.8), ncol = K), gamma = matrix(1, nrow = M, ncol = K), 
                xi = matrix(0.1, nrow = N, ncol = K), phi = rep(1, K))
    
  # MCMC sampling
  obj <- MCMC(data.obj1, data.obj2, cancer, cancer_type, 
              n.sample, M, pG, mRate = 0, K, init, delta, frailty = T, abc = T, range.t, check = F)
    
  # posterior
  post.beta  <- t(obj$posterior$beta)
  post.gamma <- t(obj$posterior$gamma)
  write.csv(post.beta, paste0("result/post_beta_seed", iter, ".csv"), row.names = FALSE)
  write.csv(post.gamma, paste0("result/post_gamma_seed", iter, ".csv"), row.names = FALSE)
  post.phi <- t(obj$posterior$phi)
  write.csv(post.phi, paste0("result/post_phi_seed", iter, ".csv"), row.names = FALSE)
  
  # acceptance ratio
  id <- 10001:n.sample
  
  beta.est <- apply(post.beta[id, ], 2, mean)
  gamma.est <- apply(post.gamma[id, ], 2, mean)
    
  temp <- list(beta = beta.est, gamma = gamma.est)
  
  temp
}  

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

# a <- 0
# for (i in 1:200) {
#   fam <- cancer[[i]]
#   id <- unique(fam$id)
#   for (j in 1:length(id)) {
#     if (nrow(fam[fam$id == id[j],]) > 2) {
#       a <- a + 1
#     }
#   }
# }














