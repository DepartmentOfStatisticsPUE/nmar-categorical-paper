## y -- target variable(s)
## x -- propensity score variables
## z -- variables for totals (calibrarion variables)
## d -- weights
## link -- link function
## controls --- currently only eps and maxiter

nmar_chang_kott <- function(Y, X, Z, d, totals, link = "logit", control) {
  
  
  if (link != "logit") {
    stop("Currently, only logit link function is implemented")
  }
  
  phi <- rep(0, ncol(X))
  k<-0
  n <- NROW(X)
  
  repeat {
    
    phi0 <- phi
    k  <- k+1
    
    rho <- as.numeric(exp(X %*% phi) / (1 + exp(X %*% phi))) 
    U <- colSums( d*(1/ rho - 1)*Z) - totals
    S <- matrix(0, ncol(Z), ncol(X))
    
    for (i in 1:ncol(Z)) {
      for (j in 1:ncol(X)) {
        S[i, j] = -sum(d *(1 - rho) * Z[,i]*X[,j] / rho)
      }
    }
    
    phi <- phi0 - MASS::ginv(t(S) %*% S) %*% (t(S) %*% U)
    
    dif <- sum((phi - phi0) ^ 2)
    
    if (dif < control$eps | k >= control$maxiter) {
      break
    }
    
  }
  rownames(phi) <- colnames(X)
  colnames(phi) <- "phi"
  
  rho <- as.numeric(exp(X %*% phi) / (1 + exp(X %*% phi))) 
  esty <- colSums(d * Y / rho) / sum(d / rho)
  
  SM <- d*(1/ rho - 1)*Z
  
  S <- matrix(0, ncol(Z), ncol(X))
  for (i in 1:nrow(S)) {
    for (j in 1:ncol(S)) {
      S[i, j] = -sum(d *(1 - rho) * Z[,i]*X[,j] / rho)
    }
  }
  
  BM1 <- matrix(0, ncol(Y), ncol(X))
  for (i in 1:nrow(BM1)) {
    for (j in 1:ncol(BM1)) {
      BM1[i,j] <-  -sum(d * as.numeric(Y[,i] - esty[i]) * X[,j] * (1 - rho) / rho)
    }
  }
  
  K <- BM1 %*% MASS::ginv(t(S) %*% S)
  
  eta <- d * (Y - esty) / rho - matrix(data = SM %*% S %*% t(K), ncol = ncol(Y), nrow = n, byrow = T)
  tau <- -sum(d / rho)
  
  varesty <- n * var(eta) / tau ^ 2
  varphi <- n * diag(MASS::ginv(t(S) %*% S) %*% t(S) %*% var(SM) %*% S %*% MASS::ginv(t(S) %*% S))
  
  return(list(phi=phi,
              phi_se = sqrt(varphi),
              y_hat = esty,
              y_se = sqrt(diag(varesty)),
              iter = k, 
              totals_met = sum(abs(U)) < 1, 
              totals_diff_sum = sum(abs(U)),
              totals_diff = U,
              rho = rho))
}

