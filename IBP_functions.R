
# Sample from Dirichlet distribution
rdirichlet <- function (n, a){
  # parameters in columns:
  a <- rbind(a)
  M <- ncol(a)
  # for a single set of parameters, extend to matrix
  if (n > nrow(a))
    a <- matrix(a, n, M, byrow = TRUE)
  x <- matrix(rgamma(M*n, a), nrow=n, ncol=M)
  # normalize to sum up to one:
  x/rowSums(x)
}

## function to compute the logarithm of the EFPF
log.EFPF <- function(W, eta, delta, gamma){
  # Number of active trees K_n
  k=dim(W)[2]
  # Sample size n
  n=dim(W)[1]
  # Number of observations belonging to each tree (vector of size K_n, with integer values between 0 and n)
  m=colSums(W)
  # First, compute the sum that appears in the exp in the EFPF
  aaa = rep(0,n)
  for (i in 1:n){
    aaa[i]=lgamma(eta + delta + i - 1) -lgamma(eta + delta) +lgamma( delta + 1 )-lgamma( delta + i )
  }
  aaa=exp(aaa)
  aaa1=sum(aaa)
  # Compute the log of each term in the product l=1 to k of the factorials  (the 10^-100 is just to avoid numerical problems)
  bbb = lgamma( m - eta ) - lgamma(1- eta) + lgamma(delta + eta + n -m  + 10^-100) - lgamma( delta + eta )
  # Compute the log of the EFPF
  # log.EFPF = - lgamma(k+1) + k * ( log(gamma) - lgamma( delta + n) + lgamma(delta +1)) - gamma * aaa1 + sum(bbb)
  log.EFPF = k * ( log(gamma) - lgamma( delta + n) + lgamma(delta +1)) - gamma * aaa1 + sum(bbb)
  return(log.EFPF)
}


### Update gamma - Close full conditional. Prior gamma \sim Ga(gamma.alpha,gamma.beta)
update_gamma <- function(W, eta,delta,gamma, gamma.alpha,gamma.beta){
  n = nrow(W)
  aaa = rep(0,n)
  for (i in 1:n){
    aaa[i]=lgamma(eta + delta + i - 1) -lgamma(eta + delta) +lgamma( delta + 1 )-lgamma( delta + i )
  }
  aaa=exp(aaa)
  aaa1=sum(aaa) 
  gamma = rgamma(1, dim(W)[2] + gamma.alpha, aaa1  + gamma.beta)
  return(gamma)
}



# ### Update delta - Multiplicative Random Walk Metropolis for eta+delta
# update_delta <- function(W,accept.delta, eta,delta,gamma, sigma.eta,sigma.delta){
# 
#   # sample  log( eta + delta^prop ) \sim Normal ( log( eta + delta^curr ) , sigma.delta )
#   prop = rnorm(1, log(delta+eta), sigma.delta)
#   delta.prop = exp(prop) - eta
# 
#   if (eta + delta.prop >0){
#   # compute acceptance probability
#   accept =  log.EFPF(W, eta, delta.prop, gamma)  -  log.EFPF(W, eta, delta, gamma) +                       # 'Likelihood' (EFPF)
#     log(dgamma(delta.prop + eta ,delta.alpha ,delta.beta ))  - log(dgamma(delta + eta,delta.alpha, delta.beta))  +   # prior
#     log( delta.prop + eta)  - log( delta + eta )                                             # Proposal
# 
#   accepttt = exp(accept)
#   acceptprob=min(c(accepttt,1))
#   bern_accept=rbinom(1,1,acceptprob)
#   if (bern_accept==1) {
#      delta = delta.prop
#      accept.delta = accept.delta+1
#      }
#   }
#   return(c(delta,accept.delta))
# }
# 
# 
# 
# ### update eta  - Multiplicative Random Walk Metropolis
# 
# update_eta <- function(W, accept.eta, eta,delta,gamma, sigma.eta,sigma.delta){
# 
#   # new parametrization
#   eta.par = 1 - eta
#   delta.par = delta + eta
# 
#   # Sample log(1-eta^prop)
# 
#   prop = rnorm(1, log( eta.par ), sigma.eta)
# 
#   eta.par.prop = exp(prop)
# 
#     # compute the log of the acceptance probability
#     accept = log.EFPF(W, 1- eta.par.prop, delta.par + eta.par.prop - 1, gamma) - log.EFPF(W, 1- eta.par, delta.par + eta.par -1 , gamma) +          # 'Likelihood' (EFPF)
#       log( dgamma( eta.par.prop,eta.alpha,eta.beta)) - log( dgamma( eta.par ,eta.alpha,eta.beta)) +   # prior,  1-eta \sim Ga(eta.alpha,eta.beta)
#       log( eta.par.prop) - log( eta.par )                                         # Additional part for Multiplicative RW
# 
#     accepttt=exp(accept)
#     acceptprob=min(c(accepttt,1))
#     bern_accept=rbinom(1,1,acceptprob)
#     if (bern_accept==1){
#       eta= 1- eta.par.prop
#       accept.eta = accept.eta +1
#     }
#   return(c(eta,accept.eta))
# }




# ### update eta  - Independent Metropolis 
# 
# # update_eta <- function(W,eta,delta,gamma, sigma.eta,sigma.delta){
# #   
# #   eta.prop = runif(1,-delta,1)
# #   
# #    # compute the log of the acceptance probability
# #     accept = log.EFPF(W, eta.prop, delta, gamma) - log.EFPF(W, eta, delta, gamma) +          # 'Likelihood' (EFPF)
# #       log( dgamma(1-eta.prop,eta.alpha,eta.beta) ) - log( dgamma(1-eta,eta.alpha,eta.beta) ) +   # prior,  1-eta \sim Ga(eta.alpha,eta.beta)
# #       log( dgamma(delta+eta.prop,delta.alpha,delta.beta) )  - log(dgamma(delta+eta,delta.alpha,delta.beta))  + # prior
# #       log( dunif(eta,-delta,1)) - log(dunif(eta.prop,-delta,1))                                         # proposal
# #     
# #     accepttt=exp(accept)
# #     acceptprob=min(c(accepttt,1))
# #     bern_accept=rbinom(1,1,acceptprob)
# #     if (bern_accept==1) eta=eta.prop 
# #     
# #   return(eta)
# # }


# Update eta - Random Walk Metropolis
update_eta <- function(W,accept.eta,eta,delta,gamma, sigma.eta,sigma.delta, eta.alpha,eta.beta,delta.alpha,delta.beta){
  eta.prop = rnorm(1,eta,sigma.eta)

  if (eta.prop < 1  & -delta < eta.prop){
  # compute the log of the acceptance probability
  accept = log.EFPF(W, eta.prop, delta, gamma) - log.EFPF(W, eta, delta, gamma) +          # 'Likelihood' (EFPF)
    log( dgamma(1-eta.prop,eta.alpha,eta.beta) ) - log( dgamma(1-eta,eta.alpha,eta.beta) ) +   # prior,  1-eta \sim Ga(eta.alpha,eta.beta)
    log( dgamma(delta+eta.prop,delta.alpha,delta.beta) )  - log(dgamma(delta+eta,delta.alpha,delta.beta))   # prior
                                       # proposal
 # accepttt=exp(accept)
#  acceptprob=min(c(accepttt,1))
#  bern_accept=rbinom(1,1,acceptprob)
  
  if (log(runif(1)) < accept) {
    eta=eta.prop
    accept.eta = accept.eta + 1
  }
  }
  return(c(eta,accept.eta))
}



### Update delta - Random Walk Metropolis
update_delta <- function(W,accept.delta,eta,delta,gamma, sigma.eta,sigma.delta, delta.alpha ,delta.beta){

  # sample  log( eta + delta^prop ) \sim Normal ( log( eta + delta^curr ) , sigma.delta )
  delta.prop  = rnorm(1, delta, sigma.delta)

  if (eta + delta.prop >0){
    # compute acceptance probability
    accept =  log.EFPF(W, eta, delta.prop, gamma)  -  log.EFPF(W, eta, delta, gamma) +                       # 'Likelihood' (EFPF)
      log(dgamma(delta.prop + eta ,delta.alpha ,delta.beta ))  - log(dgamma(delta + eta,delta.alpha, delta.beta))    # prior
                                               # Proposal
    #accepttt = exp(accept)
    #acceptprob=min(c(accepttt,1))
    #bern_accept=rbinom(1,1,acceptprob)
    if (log(runif(1)) < accept) {
      delta = delta.prop
      accept.delta = accept.delta + 1
    }
  }
  return(c(delta,accept.delta))
}


#### update delta and eta together

# Generic 1D slice sampler
slice_step <- function(x0, logpost_x, w = 1.0, m = 50) {
  fx0 <- logpost_x(x0)
  if (!is.finite(fx0)) return(x0)  # safeguard
  
  # vertical level
  logy <- fx0 - rexp(1)
  
  # Step out
  u <- runif(1, 0, w)
  L <- x0 - u
  R <- L + w
  J <- floor(runif(1, 0, m))
  K <- (m - 1) - J
  
  while (J > 0 && is.finite(logpost_x(L)) && logpost_x(L) > logy) { L <- L - w; J <- J - 1 }
  while (K > 0 && is.finite(logpost_x(R)) && logpost_x(R) > logy) { R <- R + w; K <- K - 1 }
  
  # Shrinkage
  repeat {
    x1 <- runif(1, L, R)
    fx1 <- logpost_x(x1)
    if (is.finite(fx1) && fx1 >= logy) return(x1)
    if (x1 < x0) L <- x1 else R <- x1
    if (abs(R - L) < .Machine$double.eps * max(1, abs(x0))) return(x0)
  }
}


update_eta_delta_block <- function(W,accept.eta,eta,delta,gamma, sigma.eta,sigma.delta, eta.alpha,eta.beta,delta.alpha,delta.beta,
                                   w1=1, w2=1, m=100)
{
    # --- propose independently in original parametrisation ---
  u1 <- eta
  u2 <- log(eta + delta)  # requires eta+delta > 0
  if (!is.finite(u2) || !(eta < 1)) stop("Invalid starting state")
  
  # log posterior in (u1,u2)
  logpost <- function(u1_val, u2_val){
    eta_val <- u1_val
    delta_val <- exp(u2_val) - eta_val
    if (!(eta_val < 1 && (eta_val + delta_val) > 0)) return(-Inf)
    
    ll <- log.EFPF(W, eta_val, delta_val, gamma)
    if (!is.finite(ll)) return(-Inf)
    
    one_minus_eta <- 1 - eta_val
    if (!(one_minus_eta > 0)) return(-Inf)
    pr_eta <- dgamma(one_minus_eta, shape = eta.alpha, rate = eta.beta, log = TRUE)
    
    ll + pr_eta + u2_val   # Jacobian term
  }
  
  # 1. update u1 | u2
  logpost_u1 <- function(u1_val) logpost(u1_val, u2)
  u1 <- slice_step(u1, logpost_u1, w = w1, m = m)
  
  # 2. update u2 | u1
  logpost_u2 <- function(u2_val) logpost(u1, u2_val)
  u2 <- slice_step(u2, logpost_u2, w = w2, m = m)
  
  # map back
  eta_prop   <- u1
  delta_prop <- exp(u2) - u1
  # Map back
  
   # if (log(runif(1)) < log_acc) {
      eta <- eta_prop
      delta <- delta_prop
      accept.eta <- accept.eta + 1
  #  }
    
    return(c(eta = eta, delta = delta, accept.eta= accept.eta))
  }
