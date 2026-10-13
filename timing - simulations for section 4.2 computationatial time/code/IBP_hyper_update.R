# -------------------------------------------------------------------------#
# Improved updates and elicitation for the IBP hyperparameters             #
# (eta, delta, gamma).                                                      #
#                                                                          #
# Model (as in log.EFPF):  obs i joins existing tree k w.p.                 #
#   (m_k - eta)/(i - 1 + delta),  and opens Poisson(gamma * r_i) new trees, #
#   r_i = Gamma(1+delta)Gamma(i-1+delta+eta) / (Gamma(i+delta)Gamma(delta+eta))
# Constraints: eta < 1, eta + delta > 0.                                    #
# Priors (same hyperparameters as in bart.ibp):                             #
#   1 - eta       ~ Gamma(eta.alpha,   eta.beta)                            #
#   eta + delta   ~ Gamma(delta.alpha, delta.beta)                          #
#   gamma         ~ Gamma(gamma.alpha, gamma.beta)                          #
# -------------------------------------------------------------------------#

## H_n(eta, delta) = sum_{i=1}^n r_i, in CLOSED FORM (O(1) instead of O(n)):
##   H_n = [ (eta+delta) R - delta ] / eta,
##   R   = Gamma(eta+delta+n) Gamma(delta+1) / ( Gamma(eta+delta+1) Gamma(delta+n) )
## with the eta -> 0 limit 1 + delta * (digamma(delta+n) - digamma(delta+1)).
## Written so that every lgamma argument is positive (delta may be negative).
ibp_H = function(n, eta, delta) {
  if (eta == 0) return(1 + delta * (digamma(delta + n) - digamma(delta + 1)))
  logR = if (abs(eta) < 1e-2) {
    sum(log1p(eta / (delta + seq_len(n - 1))))      # accurate when eta is near 0 (no cancellation)
  } else lgamma(eta + delta + n) - lgamma(eta + delta + 1) - lgamma(delta + n) + lgamma(delta + 1)
  exp(logR) + delta * expm1(logR) / eta
}

## log p(W | eta, delta) with gamma INTEGRATED OUT (up to constants in W).
## m = colSums(W), K = length(m).
ibp_log_marg = function(m, n, eta, delta, gamma.alpha, gamma.beta) {
  K = length(m)
  H = ibp_H(n, eta, delta)
  K * (lgamma(delta + 1) - lgamma(delta + n)) +
    sum(lgamma(m - eta) - lgamma(1 - eta) + lgamma(delta + eta + n - m) - lgamma(delta + eta)) +
    lgamma(K + gamma.alpha) - (K + gamma.alpha) * log(H + gamma.beta)
}

## Unconstrained coordinates: u1 = log(1 - eta), u2 = log(eta + delta)
ibp_to_u   = function(eta, delta) unname(c(log(1 - eta), log(eta + delta)))
ibp_from_u = function(u) { u = unname(u); eta = 1 - exp(u[1]); c(eta = eta, delta = exp(u[2]) - eta) }

## log posterior of u (collapsed over gamma), including both priors and Jacobians
ibp_log_post_u = function(u, m, n, pr) {
  if (any(!is.finite(u))) return(-Inf)
  ed = ibp_from_u(u)
  lp = ibp_log_marg(m, n, ed[["eta"]], ed[["delta"]], pr$gamma.alpha, pr$gamma.beta) +
    dgamma(exp(u[1]), pr$eta.alpha,   pr$eta.beta,   log = TRUE) + u[1] +
    dgamma(exp(u[2]), pr$delta.alpha, pr$delta.beta, log = TRUE) + u[2]
  if (is.finite(lp)) lp else -Inf
}

## Slice step with stepping out (Neal 2003), one coordinate
.slice1 = function(x0, f, w, m = 50, f0 = f(x0)) {
  logy = f0 - rexp(1)
  L = x0 - runif(1, 0, w); R = L + w
  J = floor(runif(1, 0, m)); Kk = m - 1 - J
  while (J > 0 && f(L) > logy) { L = L - w; J = J - 1 }
  while (Kk > 0 && f(R) > logy) { R = R + w; Kk = Kk - 1 }
  repeat {
    x1 = runif(1, L, R); f1 = f(x1)
    if (f1 >= logy) return(c(x1, f1))
    if (x1 < x0) L = x1 else R = x1
  }
}

## MAIN UPDATE.
##  1) (eta, delta) | W with gamma INTEGRATED OUT, by slice sampling on
##     u = (log(1-eta), log(eta+delta)) along the principal axes of the posterior
##     of u. The axes are estimated during burn-in (adapt = TRUE) and then frozen,
##     so the chain after burn-in is an exact MCMC. No tuning parameters.
##  2) gamma | eta, delta, W ~ Gamma(K + gamma.alpha, H_n + gamma.beta) (exact).
##  Steps 1+2 form an exact blocked Gibbs step for p(eta, delta | W) p(gamma | eta, delta, W).
##
##  `state` carries the adaptation; start with state = NULL and pass back what is returned.
update_ibp_hyper = function(W, eta, delta, gamma, pr, state = NULL, adapt = TRUE, n_rep = 2) {
  W = as.matrix(W); m = colSums(W); n = nrow(W)
  if (is.null(state)) state = list(t = 0, mean = c(0, 0), M2 = matrix(0, 2, 2),
                                   D = list(c(1, 0), c(0, 1), c(1, -1) / sqrt(2), c(1, 1) / sqrt(2)))
  u = ibp_to_u(eta, delta)
  f = function(uu) ibp_log_post_u(uu, m, n, pr)
  fu = f(u)
  for (r in seq_len(n_rep)) for (d in state$D) {
    s = .slice1(0, function(t) f(u + t * d), 1, f0 = fu)     # width 1 = one posterior sd once adapted
    u = u + s[1] * d; fu = s[2]
  }
  if (adapt) {                                   # Welford running covariance of u
    state$t = state$t + 1
    dlt = u - state$mean; state$mean = state$mean + dlt / state$t
    state$M2 = state$M2 + tcrossprod(dlt, u - state$mean)
    if (state$t >= 50 && state$t %% 25 == 0) {
      S = state$M2 / (state$t - 1) + diag(1e-8, 2)
      e = eigen(S, symmetric = TRUE)
      state$D = lapply(1:2, function(j) e$vectors[, j] * sqrt(e$values[j]))
    }
  }
  ed = ibp_from_u(u)
  gamma = rgamma(1, ncol(W) + pr$gamma.alpha, ibp_H(n, ed[["eta"]], ed[["delta"]]) + pr$gamma.beta)
  list(eta = ed[["eta"]], delta = ed[["delta"]], gamma = gamma, state = state)
}


# -------------------------------------------------------------------------#
# Elicitation / initialisation                                             #
# -------------------------------------------------------------------------#
## Interpretable prior quantities (exact under the model):
##   trees per observation      E[sum_k W_ik]            = gamma
##   pairwise sharing           P(obs j uses a tree of obs i) = (1 - eta)/(1 + delta) = rho
##   total number of trees      E[K_n]                   = gamma * H_n(eta, delta)
## rho = 1 is standard BART (everybody shares every tree); small rho = little sharing.
## Given (trees_per_obs, rho, total_trees) this solves for (gamma, eta, delta).
## Feasibility: total_trees must exceed about trees_per_obs / rho.
ibp_elicit = function(n, trees_per_obs, rho, total_trees) {
  stopifnot(rho > 0, rho < 1, total_trees > trees_per_obs)
  gamma = trees_per_obs
  target = total_trees / gamma
  delta_of = function(eta) (1 - eta) / rho - 1
  g = function(t) { eta = 1 - exp(t); log(ibp_H(n, eta, delta_of(eta))) - log(target) }  # t = log(1-eta)
  lo = log(1e-6); hi = log(1e6)
  if (g(hi) > 0) stop(sprintf("total_trees too small: with rho = %.3f the minimum is about %.1f",
                              rho, gamma * ibp_H(n, 1 - 1e6, delta_of(1 - 1e6))))
  if (g(lo) < 0) stop("total_trees too large for this rho (would need eta >= 1).")
  eta = 1 - exp(uniroot(g, c(lo, hi), tol = 1e-10)$root)
  c(gamma = gamma, eta = eta, delta = delta_of(eta))
}

## Gamma priors centred (in mean) on elicited values, with coefficient of variation cv.
## Returns the six prior hyperparameters in the naming used by bart.ibp.
ibp_prior_from_elicit = function(el, cv = 0.5) {
  shape = 1 / cv^2
  one_m_eta = 1 - el[["eta"]]; ed = el[["eta"]] + el[["delta"]]
  list(eta.alpha = shape,   eta.beta = shape / one_m_eta,
       delta.alpha = shape, delta.beta = shape / ed,
       gamma.alpha = shape, gamma.beta = shape / el[["gamma"]])
}

## Summary of what given hyperparameters imply (useful for checking a setting)
ibp_implied = function(n, eta, delta, gamma) {
  eta = unname(eta); delta = unname(delta); gamma = unname(gamma)
  c(trees_per_obs = gamma, rho = (1 - eta) / (1 + delta), total_trees = gamma * ibp_H(n, eta, delta))
}
