Pest <- function(mxy) {
  return(colSums(mxy) / (ny + colSums(mxy)))
}

#Log likelihood based on Pest
ll <- function(mxy) {
  r <- (mxy) / (ny + (mxy))
  sum(log(r)) * 10000
}

#Estimation
Eest <- function(mxy) {
  res <- t(Pest(mxy) * t(as.matrix(nxy + mxy)))
  return(res)
  
}
#Maximization
mxy <- function(mxy, Est) {
  res <- (mx * Eest(mxy)) / (rowSums(Est))
  return(res)
  
}



#maxiter <- 100000

EM <- function(mxy_init,
               maxit = maxiter,
               tol = 0.002,
               iter = F,
               ad = F) {
  flag <- 0
  
  mxy_cur <- mxy_init
  
  for (i in 1:maxit) {
    ll1 <- ll(mxy_cur)
    
    new <- mxy(mxy_cur, Eest(mxy_cur))
    
    ll2 <- ll(new)
    
    
    if (abs(ll2 - ll1) < tol) {
      flag <- 1
      
      
      if (ad) {
        return(c(abs(ll1), abs(ll2)))
      }
      
      print(paste("iteration #", i))
      
      
      if (!iter) {
        print(paste("Estimated mxy"))
        print(round(mxy_cur))
        print(paste("Actual mxy"))
        print(realmxy)
        
        print(paste("Col sums"))
        print(round(colSums(mxy_cur)))
        print(paste("Actual col sums"))
        print(colSums(realmxy))
      }
      
      break
    }
    
    
    mxy_cur <- new
    
  }
  if (ad) {
    return(c(abs(ll1), abs(ll2)))
  }
  
  if (!iter) {
    if (!flag) {
      warning("Didn't converge\n")
      print(paste("iteration #", i))
      print(paste("Estimated mxy"))
      print(round(mxy_cur))
      print(paste("Actual mxy"))
      print(realmxy)
      
      print(paste("Col sums"))
      print(round(colSums(mxy_cur)))
      print(paste("Actual col sums"))
      print(colSums(realmxy))
      
    }
    print(paste("ll1,ll2"))
    print(c(ll1, ll2))
    return(mxy_cur)
  }
  if (iter) {
    return(i)
  }
  
}

#Cascade iteration to find the lowest tol level, for which EM converges
tolE <- function(maxiter = maxiter) {
  tol1 <- 50
  dlt <- c(0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 1)
  
  ite <- c(100, 20)
  
  
  t1 <- Sys.time()
  
  
  for (i in 1:1000) {
    it <- EM(mxy_init,
             maxit = maxiter,
             tol = tol1,
             iter = T)
    
    if (it != maxiter) {
      tol1 = tol1 -  dlt[8]
    }
    else{
      tol1 = tol1 +  dlt[8]
      print(paste("iteration0 #", i))
      break
    }
  }
  
  for (i in 1:1000) {
    it <- EM(mxy_init,
             maxit = maxiter,
             tol = tol1,
             iter = T)
    
    if (it != maxiter) {
      tol1 = tol1 -  dlt[1]
    }
    else{
      print(paste("iteration1 #", i))
      break
    }
  }
  
  
  for (i in 1:20) {
    it <- EM(mxy_init,
             maxit = maxiter,
             tol = tol1,
             iter = T)
    
    if (it == maxiter) {
      tol1 = tol1 +  dlt[2]
    }
    else{
      tol1 = tol1 -  dlt[2]
      print(paste("iteration2 #", i))
      break
    }
  }
  
  
  for (i in 1:20) {
    it <- EM(mxy_init,
             maxit = maxiter,
             tol = tol1,
             iter = T)
    
    if (it == maxiter) {
      tol1 = tol1 +  dlt[3]
    }
    else{
      tol1 = tol1 -  dlt[3]
      print(paste("iteration3 #", i))
      break
    }
  }
  
  for (i in 1:20) {
    it <- EM(mxy_init,
             maxit = maxiter,
             tol = tol1,
             iter = T)
    
    if (it == maxiter) {
      tol1 = tol1 +  dlt[4]
    }
    else{
      tol1 = tol1 -  dlt[4]
      print(paste("iteration4 #", i))
      break
    }
  }
  
  for (i in 1:20) {
    it <- EM(mxy_init,
             maxit = maxiter,
             tol = tol1,
             iter = T)
    
    if (it == maxiter) {
      tol1 = tol1 +  dlt[5]
    }
    else{
      tol1 = tol1 -  dlt[5]
      print(paste("iteration4 #", i))
      break
    }
  }
  for (i in 1:20) {
    it <- EM(mxy_init,
             maxit = maxiter,
             tol = tol1,
             iter = T)
    
    if (it == maxiter) {
      tol1 = tol1 +  dlt[6]
    }
    else{
      tol1 = tol1 -  dlt[6]
      print(paste("iteration4 #", i))
      break
    }
  }
  for (i in 1:20) {
    it <- EM(mxy_init,
             maxit = maxiter,
             tol = tol1,
             iter = T)
    
    if (it == maxiter) {
      tol1 = tol1 +  dlt[7]
    }
    else{
      print(paste("iteration4 #", i))
      break
    }
  }
  
  
  
  Sys.time() - t1
  return(tol1)
}

ad <- function(tol,maxiter){
  
  
  v1 <- EM(mxy_init,maxit = maxiter,tol = tol,ad = T)
  
  res <- sum(all(trunc(diff(v1)*10)/10 == 0))
  
  return(res)
}

