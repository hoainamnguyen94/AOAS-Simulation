data.gen <- function(seed, N, pi, beta, lambda, theta.c, phi, miss = T) {
  set.seed(seed)
  
  count <- count.mut <- 0
  data.obj1 <- data.obj2 <- cancer.obj <- as.list(1:N)
  
  while (count < N) {
    proband <- get.npp(pi, beta, lambda, phi, theta.c)
    
    if (nrow(proband) > 2) {
      data1 <- data.frame(ID = 1:30, 
                          FatherID = c(3,5,0,0,0,0,
                                       1,1,1,1,
                                       3,3,5,5,
                                       0,0,0,0,
                                       11,11,11,
                                       16,16,16,
                                       13,13,13,
                                       18,18,18), 
                          MotherID = c(4,6,0,0,0,0,
                                       2,2,2,2,
                                       4,4,6,6,
                                       0,0,0,0,
                                       15,15,15,
                                       12,12,12,
                                       17,17,17,
                                       14,14,14), 
                          Gender = c(1,2,1,2,1,2,
                                     1,2,1,2,
                                     1,2,1,2,
                                     2,1,2,1,
                                     1,2,1,2,1,2,2,1,2,1,2,1))
      n1 <- nrow(data1)
      # ID
      ID <- data1$ID
      # gender
      S <- data1$Gender
      # genotype
      G <- rep(0, n1)
      G[1] <- unique(proband[, "g"])
      if (G[1] == 1) {
        # spouse family
        spouse.index <- c(2,5,6,12,14,15,16,17,18,25,26,27,28,29,30)
        G[spouse.index] <- 0
        # parents
        G[3] <- rbinom(1, 1, 1/2)
        G[4] <- 1 - G[3]
        # offspring
        G[7:10] <- rbinom(4, 1, 1/2)
        # sibling
        G[11:12] <- rbinom(2, 1, 1/2)
        # sib-offspring1
        if (G[11] == 1) G[19:21] <- rbinom(3,1,1/2)
        # sib-offspring2
        if (G[12] == 1) G[22:24] <- rbinom(3,1,1/2)
      }
      
      # generate data
      d.1 <- dp.1 <- proband[, "dp.1"]
      d.2 <- dp.2 <- proband[, "dp.2"]
      if (length(d.1) > 1) d.1[2] <- 0
      if (length(d.2) > 1) d.2[2] <- 0
      
      data2 <- cbind(ID = rep(ID[1], nrow(proband)), 
                     time = proband[, "time"],
                     test = proband[, "g"],
                     D.1 = d.1, D.2 = d.2, 
                     Dp.1 = dp.1, Dp.2 = dp.2)
      
      cancer <- cbind(id = rep(ID[1], nrow(proband) - 1),
                      cancer.type = proband[1:(nrow(proband) - 1), "ct"],
                      diag.age = proband[1:(nrow(proband) - 1), "time"])
      
      for (k in 2:n1) {
        temp <- get.npp(pi, beta, lambda, phi, theta.c, as.numeric(proband[1, c("xi.1", "xi.2")]), 
                        G[k], S[k])
        
        d.1 <- dp.1 <- temp[, "dp.1"]
        d.2 <- dp.2 <- temp[, "dp.2"]
        if (length(d.1) > 1) d.1[2] <- 0
        if (length(d.2) > 1) d.2[2] <- 0
        
        Temp <- cbind(ID = rep(ID[k], nrow(temp)), 
                      time = temp[, "time"],
                      test = temp[, "g"], 
                      D.1 = d.1, D.2 = d.2, 
                      Dp.1 = dp.1, Dp.2 = dp.2)
        
        data2 <- rbind(data2, Temp)
        
        if (nrow(temp) > 1){
          cancer.Temp <- cbind(id = rep(ID[k], nrow(temp) - 1),
                               cancer.type = temp[1:(nrow(temp) - 1), "ct"],
                               diag.age = temp[1:(nrow(temp) - 1), "time"])
          
          cancer <- rbind(cancer, cancer.Temp)
        }
      }
      
      rownames(data2) <- NULL
      rownames(cancer) <- NULL
      
      # generate missing (20%)
      if (miss)
      {
        m.index <- sample(ID[-1], 6)
        na.index <- !is.na(match(data2[, 1], m.index))
        data2[na.index, "test"] <- NA
      }
      
      if (proband[,"g"][1] == 1) {
        count.mut <- count.mut+1
      }
      
      count <- count + 1
      print(count)
      
      data.obj1[[count]] <- data1
      data.obj2[[count]] <- data2
      cancer.obj[[count]] <- data.frame(cancer)
    }
  }
  
  print(count.mut)
  list(data.obj1 = data.obj1, data.obj2 = data.obj2, cancer.obj = cancer.obj)
}






