posterior.update <- function(beta, gamma, phi, xi, llike, n.family, data.obj1, data.obj2, cancer, cancer_type,
                             sigma.beta, a.phi, b.phi, delta.beta, delta.gamma, delta.phi, delta.xi, 
                             prelim.data, prelim.proband, frailty, abc, G, nG, M, K, 
                             range.t, allef, mRate) 
{
  p <- nrow(beta)
  ac.beta <- matrix(0, nrow = p, ncol = K)
  ac.gamma <- matrix(0, nrow = M, ncol = K)
  ac.xi <- matrix(0, nrow = n.family, ncol = K)
  ac.phi <- rep(0, K)
  
  for (j in 1:K){
    # update beta
    for (i in 1:p){
      new.beta <- beta
      # draw new value
      new.beta[i, j] <- rnorm(1, mean = beta[i, j], sd = delta.beta[i, j])    
      
      # prior likelihood
      lprior <- dnorm(beta[i, j], mean = 0, sd = sigma.beta, log = T)
      new.lprior <- dnorm(new.beta[i, j], mean = 0, sd = sigma.beta, log = T)
      
      # likelihood
      llike.pro <- rep(0, n.family) 
      if (abc) {
        llike.pro <- proband.likelihood(new.beta, gamma, xi, prelim.proband, M, K, range.t, allef)
      }
      new.llike <- FamilyLikelihood(new.beta, gamma, xi, n.family, data.obj1, data.obj2, cancer, cancer_type,
                                    prelim.data, llike.pro, M, K, range.t, allef, mRate)
      
      # update sample
      r <- exp((new.llike + new.lprior) - (llike + lprior))
      if (r > runif(1)) 
      {
        beta[i, j] <- new.beta[i, j]
        ac.beta[i, j] <- 1
        llike <- new.llike
      }
    }
    
    # update gamma
    for (i in 1:M){
      new.gamma <- gamma
      # draw new value
      tmp <- gamma[i, j] + rnorm(50, 0, min(c(2 * gamma[i, j], delta.gamma[i, j])))
      tmp.id <- min(which(tmp > 0))
      new.gamma[i, j] <- tmp[tmp.id]
      
      # likelihood
      llike.pro <- rep(0, n.family)
      if (abc) {
        llike.pro <- proband.likelihood(beta, new.gamma, xi, prelim.proband, M, K, range.t, allef)
      }
      new.llike <- FamilyLikelihood(beta, new.gamma, xi, n.family, data.obj1, data.obj2, cancer, cancer_type,
                                    prelim.data, llike.pro, M, K, range.t, allef, mRate)
      
      # update sample
      r <- exp(new.llike - llike)
      if (r > runif(1)) {
        gamma[i, j] <- new.gamma[i, j]
        ac.gamma[i, j] <- 1
        llike <- new.llike
      }
    }
    
    if (frailty){
      #############
      # update xi #
      #############
      
      tmp.likelihood <- FamilyLikelihood(beta, gamma, xi, n.family, data.obj1, data.obj2, cancer, cancer_type, 
                                         prelim.data, llike.pro, M, K, range.t, allef, mRate, sum = F)
      for (f in 1:n.family){
        n1 <- prelim.data$n1[[f]]
        n2 <- prelim.data$n2[[f]]
        
        data1 <- data.obj1[[f]]
        data2 <- prelim.data$obj[[f]]
        
        cancer_f <- cancer[[f]]
        
        Ft <- prelim.data$Ft[[f]]
        ft <- prelim.data$ft[[f]]
        
        # draw xi
        new.xi <- xi
        tmp <- xi[f, j] + rnorm(50, 0, delta.xi[f, j])
        tmp.id <- min(which(tmp > 0))
        new.xi[f, j] <- tmp[tmp.id]
        
        f.new.xi <- new.xi[f, ]
        
        # compute prior
        lprior <- dgamma(xi[f, j], shape = phi[j], rate = phi[j], log = T)
        new.lprior <- dgamma(new.xi[f, j], shape = phi[j], rate = phi[j], log = T)
        
        # compute likelihood
        f.llike.pro <- rep(0, n.family)
        if (abc) {
          f.llike.pro <- proband.likelihood(beta, gamma, new.xi, prelim.proband, M, K, range.t, allef)
        }
        f.new.llike.unadjusted <- fllike(beta, gamma, f.new.xi, n1, n2, data1, data2, cancer_f, cancer_type,
                                         Ft, ft, G, nG, M, K, range.t, allef, mRate)
        
        f.llike <- tmp.likelihood[f]
        f.new.llike <- f.new.llike.unadjusted - f.llike.pro[f]
        
        # decision 
        r <- exp((f.new.llike + new.lprior) - (f.llike + lprior))
        if (r > runif(1)) {
          xi[f, j] <- new.xi[f, j]
          ac.xi[f, j] <- 1
          tmp.likelihood[f] <- f.new.llike
        }
      }
      llike <- sum(tmp.likelihood)
      
      ##############
      # update phi #
      ##############
      
      new.phi <- phi
      
      # draw sample
      new.phi[j] <- exp(log(phi[j]) + rnorm(1, 0, delta.phi[j]))
      
      # compute posterior
      post <- (n.family * phi[j] + a.phi - 1) * log(phi[j]) - n.family * lgamma(phi[j]) - 
        b.phi * phi[j] + phi[j] * (sum(log(xi[, j])) - sum(xi[, j]))
      new.post <- (n.family * new.phi[j] + a.phi - 1) * log(new.phi[j]) - n.family * lgamma(new.phi[j]) - 
        b.phi * new.phi[j] + new.phi[j] * (sum(log(xi[, j])) - sum(xi[, j]))
      
      ## decision
      r <- exp(new.post - post)
      if (r > runif(1)) {
        phi[j] <- new.phi[j]
        ac.phi[j] <- 1
      }
    } # end frailty
  }
  
  obj <- list(llike = llike, 
              beta = beta, gamma = gamma, phi = phi, xi = xi,
              ac.beta = ac.beta, ac.gamma = ac.gamma, ac.phi = ac.phi, ac.xi = ac.xi)
  
  return(obj)
}
