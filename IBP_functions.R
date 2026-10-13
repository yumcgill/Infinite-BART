
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

## Sum_{i=1}^n Gamma(eta+delta+i-1) Gamma(delta+1) / (Gamma(eta+delta) Gamma(delta+i))
## OPTIMISED: vectorised (the original used an R loop over i = 1..n)
ibp_harmonic = function(n, eta, delta) {
  i = seq_len(n)
  sum(exp(lgamma(eta + delta + i - 1) - lgamma(eta + delta) + lgamma(delta + 1) - lgamma(delta + i)))
}

## log EFPF from the column counts m = colSums(W) (so that the counts are not
## recomputed at every evaluation inside the slice sampler)
log.EFPF.counts <- function(m, n, eta, delta, gamma){
  k = length(m)
  aaa1 = ibp_harmonic(n, eta, delta)
  bbb = lgamma( m - eta ) - lgamma(1- eta) + lgamma(delta + eta + n -m  + 10^-100) - lgamma( delta + eta )
  k * ( log(gamma) - lgamma( delta + n) + lgamma(delta +1)) - gamma * aaa1 + sum(bbb)
}

## function to compute the logarithm of the EFPF (same interface as before)
log.EFPF <- function(W, eta, delta, gamma){
  W = as.matrix(W)
  log.EFPF.counts(colSums(W), nrow(W), eta, delta, gamma)
}


### Update gamma - Close full conditional. Prior gamma \sim Ga(gamma.alpha,gamma.beta)
update_gamma <- function(W, eta,delta,gamma, gamma.alpha,gamma.beta){
  W = as.matrix(W)
  aaa1 = ibp_harmonic(nrow(W), eta, delta)
  rgamma(1, ncol(W) + gamma.alpha, aaa1  + gamma.beta)
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
# OPTIMISED: the original evaluated logpost_x(L) / logpost_x(R) twice per
# step-out iteration; each value is now computed once. Same algorithm.
slice_step <- function(x0, logpost_x, w = 1.0, m = 50) {
  fx0 <- logpost_x(x0)
  if (!is.finite(fx0)) return(x0)  # safeguard

  logy <- fx0 - rexp(1)

  u <- runif(1, 0, w)
  L <- x0 - u
  R <- L + w
  J <- floor(runif(1, 0, m))
  K <- (m - 1) - J

  if (J > 0) { fL <- logpost_x(L)
    while (J > 0 && is.finite(fL) && fL > logy) { L <- L - w; J <- J - 1; if (J > 0) fL <- logpost_x(L) } }
  if (K > 0) { fR <- logpost_x(R)
    while (K > 0 && is.finite(fR) && fR > logy) { R <- R + w; K <- K - 1; if (K > 0) fR <- logpost_x(R) } }

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
  W = as.matrix(W)
  m_counts = colSums(W)      # computed once, not at every log-posterior evaluation
  n = nrow(W)

  u1 <- eta
  u2 <- log(eta + delta)  # requires eta+delta > 0
  if (!is.finite(u2) || !(eta < 1)) stop("Invalid starting state")

  logpost <- function(u1_val, u2_val){
    eta_val <- u1_val
    delta_val <- exp(u2_val) - eta_val
    if (!(eta_val < 1 && (eta_val + delta_val) > 0)) return(-Inf)
    ll <- log.EFPF.counts(m_counts, n, eta_val, delta_val, gamma)
    if (!is.finite(ll)) return(-Inf)
    pr_eta <- dgamma(1 - eta_val, shape = eta.alpha, rate = eta.beta, log = TRUE)
    ll + pr_eta + u2_val   # Jacobian term
  }

  u1 <- slice_step(u1, function(v) logpost(v, u2), w = w1, m = m)
  u2 <- slice_step(u2, function(v) logpost(u1, v), w = w2, m = m)

  eta   <- u1
  delta <- exp(u2) - u1
  accept.eta <- accept.eta + 1

  return(c(eta = eta, delta = delta, accept.eta= accept.eta))
}
