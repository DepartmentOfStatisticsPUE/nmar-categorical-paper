##  helper

nonignorable <- function(target = y ~ x1 + x2, 
                         missing = resp ~ 1, 
                         data = df, 
                         weights = NULL,
                         totals = NULL, 
                         family = "multinomial", 
                         model = "lee-marsh-2",
                         sparse = FALSE) {
  
  mf <- as.Formula(target, missing)
  mm <- model.frame(formula = mf, data = data)
  X <- model.matrix(mf, data = data, rhs = 1)
  Z <- model.matrix(mf, data = data, rhs = 2)
  R <- as.vector(as.matrix(model.part(mf, data = mm, lhs = 2)))
  Y <- as.vector(as.matrix(model.part(mf, data = mm, lhs = 1)))
  Y[R == 0] <- -999
  
  if (family == "multinomial") {
    Y <- factor(Y)
    Y <- as.matrix(model.matrix(~ -1 + Y))
    Y <- Y[,-1]
  }
  
  w <- weights
  
  if (is.null(weights)) {
    w <- rep(1, nrow(data))
  } 
  
  return(list(Y = Y, R = R,  X = X, Z = Z, weights = w))
}

result <- nonignorable(target = y ~ x1 + x2, missing = resp ~ 1, data = df)


## multinomial model from Lee and Marsh

model3_ll <- function(par, data) {
  
  y_obs <- data$R == 1
  y_miss <- data$R == 0
  X <- data$X
  Y <- data$Y*y_obs
  w <- data$weights
  par_end <- NROW(par)
  y_k <- ncol(Y)
  x_k <- ncol(X)
  beta_pars <- x_k*(y_k - 1)
  start_alpha <- par[(beta_pars+1):par_end]^2 
  
  eta_1 <- 0
  eta_k <- list()
  k <- 1
  for (i in 1:(y_k-1)) {
    eta_k[[i]] <- X %*% par[k:(k + x_k-1)]
    k <- k + y_k
  }
  
  etas <- cbind(eta_1, do.call("cbind", eta_k))
  sumex <- rowSums(exp(etas))
  
  
  alphas <- matrix(data = start_alpha, ncol = NROW(start_alpha), nrow = nrow(etas), byrow = T)
  sumcexc <- rowSums(alphas*exp(etas)/(1+alphas))
  loglike <- y_miss*log(sumcexc) +  rowSums(Y*(etas - log(1+alphas))) - log(sumex) 
  ll <- sum(loglike)
  
  return(ll)
  
}


model3_ll_grad <- function(par, data) {
  
  
  y_obs <- data$R == 1
  y_miss <- data$R == 0
  X <- data$X
  Y <- data$Y*y_obs
  w <- data$weights
  par_end <- NROW(par)
  y_k <- ncol(Y)
  x_k <- ncol(X)
  beta_pars <- x_k*(y_k - 1)
  start_alpha <- par[(beta_pars+1):par_end]^2 
  
  eta_1 <- 0
  eta_k <- list()
  k <- 1
  for (i in 1:(y_k-1)) {
    eta_k[[i]] <- X %*% par[k:(k + x_k-1)]
    k <- k + y_k
  }
  
  
  etas <- cbind(eta_1, do.call("cbind", eta_k))
  
  alphas <- matrix(data = start_alpha, ncol = NROW(start_alpha), nrow = nrow(X), byrow = T)
  alphas_frac <- alphas/(1+alphas)
  alphas_pow2 <- 1/(1+alphas)^2
  
  ## betas deriv
  probs <- exp(etas) / rowSums(exp(etas))
  probs_new <- ( exp(etas) * alphas_frac ) / rowSums( exp(etas) * alphas_frac )
  
  
  probs_obs <- as.matrix(Y[,-1] - probs[,-1])*y_obs
  beta_grad_obs <- t(probs_obs) %*% (X*w)
  beta_grad_obs_t <- as.vector(t(beta_grad_obs))
  
  probs_miss <- as.matrix(probs_new[,-1] - probs[,-1])*y_miss
  beta_grad_miss <- t(probs_miss) %*% (X*w)
  beta_grad_miss_t <- as.vector(t(beta_grad_miss))
  
  beta_grad <- beta_grad_obs_t+beta_grad_miss_t
  
  ## alpha
  probs_alpha <- ( exp(etas) * alphas_pow2 ) / rowSums( exp(etas) * alphas_frac )
  alphas_grad <- y_obs*as.matrix(Y*(-1/(1 + alphas))) + y_miss*probs_alpha
  
  grads <- c(beta_grad, colSums(alphas_grad))
  return(grads)
  
}


model3_result <- maxLik(logLik = model3_ll, 
                        grad = model3_ll_grad,
                        ## good starting points is the key
                        start = c(coef(model2_result)[1:6], rep(coef(model2_result)[7],3)),
                        data = result,
                        #df = df,
                        #resp = df$resp,
                        #weight = rep(, nrow(df)),
                        control = list(printLevel = 2),
                        method = "NR")

summary(model3_result)

