MCMC <- function(data.obj1, data.obj2, cancer, cancer_type,
                 n.sample, M, pG, mRate, K, init, delta, frailty, abc, range.t, check = T){
   ########################
   # mutation information #
   ########################
   maf.rate <- pG
   allef <- list(c(1 - maf.rate, maf.rate))
   nloci <- 1
   G <- c(0, 1, 1)
   nG <- length(G)
   
   ###################
   # hyper parameter #
   ###################
   sigma.beta <- sigma.alpha <- 20  # for beta prior
   a.phi <- b.phi <- 1              # for phi prior

   n.family <- length(data.obj1) # of families: 189
  
   # prepare_data 
   prelim.data <- pre.data(data.obj1, data.obj2, range.t, M) 
   prelim.proband <- pre.proband(data.obj2, range.t, M)
  
   # initial value
   beta  <- init$beta 
   gamma <- init$gamma
   phi   <- init$phi
   xi    <- init$xi
  
   if (!frailty) {
    xi <- matrix(rep(1, n.family * K), ncol = K)
    phi <- rep(1, K)
   }
   
  llike.pro <- rep(0, n.family)
  if (abc) {
    llike.pro <- proband.likelihood(beta, gamma, xi, prelim.proband, M, K, range.t, allef)
  }
  llike <- FamilyLikelihood(beta, gamma, xi, n.family, data.obj1, data.obj2, cancer, cancer_type, prelim.data, 
                            llike.pro, M, K, range.t, allef, mRate, sum = T)
  
  # step   
  delta.beta <- delta$beta
  delta.gamma <- delta$gamma
  delta.phi <- delta$phi
  delta.xi <- delta$xi
  
  # storage
  p <- nrow(beta)
  post.beta  <- acpt.beta  <- matrix(0, p * K, n.sample)
  post.gamma <- acpt.gamma <- matrix(0, M * K, n.sample)
  
  if (frailty)
  {
    post.phi <- acpt.phi <- matrix(0, K, n.sample)
  }
  
  likelihood <- NULL
  
  # Start: MCMC
  set.seed(1)
  
  for (iter in 1:n.sample) {
    print(iter)
    ttt <- Sys.time()
    updated <- posterior.update(beta, gamma, phi, xi, llike, n.family, data.obj1, data.obj2, cancer, cancer_type,
                                sigma.beta, a.phi, b.phi, delta.beta, delta.gamma, delta.phi, delta.xi, 
                                prelim.data, prelim.proband, frailty, abc, G, nG, M, K, range.t, allef, mRate)
    
    beta  <- updated$beta 
    post.beta[, iter]  <- c(updated$beta)
    gamma  <- updated$gamma
    post.gamma[, iter] <- c(updated$gamma)
    print(beta)
    
    acpt.beta[, iter]  <- c(updated$ac.beta)
    acpt.gamma[, iter] <- c(updated$ac.gamma)
    
    if (frailty)   {
      post.phi[, iter]  <- phi <- updated$phi
      acpt.phi[, iter]  <- updated$ac.phi
      xi <- updated$xi
    }
    
    likelihood[iter]  <- llike <- updated$llike
    if (check) 
    {
      cat(iter, " takes ", round(Sys.time() - ttt, 3), ", like: ", llike, ".\n", sep = "")
      print(round(beta, 3))
    }
  } # End: MCMC
  
  if (frailty) 
  {
    posterior <- list(beta = post.beta, 
                      gamma = post.gamma, 
                      phi = post.phi)
    
    acpt.ratio <- list(beta  = apply(acpt.beta , 1, mean), 
                       gamma = apply(acpt.gamma, 1, mean),
                       phi   = mean(acpt.phi))
  } else {
    posterior <- list(beta = post.beta, 
                      gamma = post.gamma)
    acpt.ratio <- list(beta  = apply(acpt.beta , 1, mean), 
                       gamma = apply(acpt.gamma, 1, mean))
  } 
    
  obj <- list(posterior = posterior, acpt.ratio = acpt.ratio, likelihood = likelihood, range.t = range.t)
}
