sourceDir <- function(path, trace = TRUE, ...) {
  for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
    if(trace) cat(nm,":")
    source(file.path(path, nm), ...)
    if(trace) cat("\n")
  }
}

sourceDir("fn")

n.sim <- 16

K <- 2
N <- 300
phi <- rep(1, K)

# estimation setting
n.sample <- 36000 # of posterior samples
M <- 2 # degrees of Bernstein polynomials

# Set-up
true.beta <- beta <- matrix(c(5, 1, 2, 4, 2, 1), nrow = 3) # regression coefficients (test, D_1, D_2)
lambda <- c(0.001, 0.001) # baseline intensity
theta.c <- 10  # hazard for censoring time
pG <- pi <- 0.001 # maf prevalence

cancer_type <- c(1, 2)

temp <- data.gen(seed = 1, N, pi, beta, lambda, theta.c, phi, miss = T)
  
data.obj1 <- temp$data.obj1 # disease history
data.obj2 <- temp$data.obj2 # pedigree
  
range.t <- max(unlist(lapply(data.obj2, function(x) max(x[, 2])))) + 1.0e-3

write(range.t, "result/range_t.txt")
