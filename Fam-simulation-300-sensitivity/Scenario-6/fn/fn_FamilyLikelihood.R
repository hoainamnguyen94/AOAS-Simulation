FamilyLikelihood <- function(beta, gamma, xi, n.family, data.obj1, data.obj2, cancer, cancer_type,
                             prelim.data, llike.pro, M, K, range.t, allef, mRate, sum = T)
{
  unadjusted <- NULL
  
  for (f in 1:n.family) {
    n1 <- prelim.data$n1[[f]]
    n2 <- prelim.data$n2[[f]]
    
    data1 <- data.obj1[[f]]
    data2 <- prelim.data$obj[[f]]
    
    cancer_f <- cancer[[f]]
    
    Ft <- prelim.data$Ft[[f]]
    ft <- prelim.data$ft[[f]]
    
    xii <- xi[f, ]
    
    unadjusted[f] <- fllike(beta, gamma, xii, n1, n2, data1, data2, cancer_f, cancer_type,
                            Ft, ft, G, nG, M, K, range.t, allef, mRate)
    
  }
  obj <- if (sum == T) sum(unadjusted - llike.pro)
  else unadjusted - llike.pro
  return(obj)
}
  