
# library(Matrix)   # not used by any function here; loading it slowed down base functions (S4 dispatch)


### UPDATE sigma2_mu

# prior half cauchy
log_prior_halfcauchy <- function(s, scale) -log1p((s / scale)^2)

## Random-walk MH on log(sigma_mu); target includes Jacobian sigma_mu
update_sigma2_mu <- function(sig_mu, mu_list, step = 0.1, log_prior = NULL) {
  K <- length(mu_list)
  if (is.null(log_prior))
    log_prior <- function(s) log_prior_halfcauchy(s, scale = 0.25 / sqrt(K))
  mu <- unlist(mu_list, use.names = FALSE); L <- length(mu); ssq <- sum(mu^2)
  lt <- function(s) -L * log(s) - ssq / (2 * s^2) + log_prior(s) + log(s)
  
  sig_prop <- exp(log(sig_mu) + step * rnorm(1))
  accepted <- log(runif(1)) < lt(sig_prop) - lt(sig_mu)
  list(sig_mu = if (accepted) sig_prop else sig_mu, accepted = accepted)
}






## =====================================================================
##  MH update of (alpha, beta) for trees encoded as
##    tree$tree_matrix  : one row per node, columns terminal, child_left,
##                        child_right, parent, split_variable, split_value,
##                        mu, node_size
##    tree$node_indices : for each observation, the row of its leaf
##    tree$var          : not used here
##
##  trees : list of K such trees (the current state of the sampler)
##  Prior on tree shape: node at depth d is split with prob alpha (1+d)^(-beta)
## =====================================================================

## ---------------------------------------------------------------------
##  Depth of every node (root = 0), computed from the parent column.
##  Does not assume that parents appear before children in the matrix,
##  so it stays correct after prune moves reorder or delete rows.
## ---------------------------------------------------------------------
node_depths <- function(tm) {
  parent <- tm[, "parent"]
  depth  <- rep(NA_integer_, nrow(tm))
  depth[is.na(parent)] <- 0L
  while (anyNA(depth)) {
    todo <- which(is.na(depth) & !is.na(depth[parent]))
    if (!length(todo)) stop("tree_matrix has nodes not connected to the root")
    depth[todo] <- depth[parent[todo]] + 1L
  }
  depth
}

## ---------------------------------------------------------------------
##  Is each leaf splittable (has at least one valid splitting rule)?
##  A leaf with no valid rule has p_split = 0 under the CGM prior and must
##  NOT contribute a (1 - p) factor.
##
##  With X: a leaf is splittable if some variable admits a cut leaving at
##  least `min_size` observations on each side, i.e. if the sorted values
##  x_(min_size) < x_(n - min_size + 1).
##  Without X: falls back to node_size >= 2 * min_size (ignores ties).
##
##  Use the SAME rule your grow/prune moves use; if your tree MH ignores
##  splittability altogether, set use_splittable = FALSE in depth_counts().
## ---------------------------------------------------------------------
leaf_splittable <- function(tree, X = NULL, min_size = 1L) {
  tm     <- tree$tree_matrix
  leaves <- which(tm[, "terminal"] == 1)
  if (is.null(X))
    return(setNames(tm[leaves, "node_size"] >= 2 * min_size, leaves))
  
  obs_by_leaf <- split(seq_along(tree$node_indices),
                       factor(tree$node_indices, levels = leaves))
  vapply(obs_by_leaf, function(idx) {
    nl <- length(idx)
    if (nl < 2 * min_size) return(FALSE)
    for (v in seq_len(ncol(X))) {
      xs <- sort(X[idx, v])
      if (xs[min_size] < xs[nl - min_size + 1]) return(TRUE)
    }
    FALSE
  }, logical(1))
}

## ---------------------------------------------------------------------
##  Pooled counts over all K trees, by depth:
##    n_int[d+1]  = number of internal (split) nodes at depth d
##    n_leaf[d+1] = number of splittable leaves at depth d
## ---------------------------------------------------------------------
depth_counts <- function(trees, X = NULL, min_size = 1L, use_splittable = TRUE, W = NULL) {
  per_tree <- lapply(seq_along(trees), function(k) {
    tr    <- trees[[k]]
    tm    <- tr$tree_matrix
    depth <- node_depths(tm)
    term  <- tm[, "terminal"] == 1
    keep  <- term
    if (use_splittable) {
      # BUG FIX: node_indices refer to the observations that USE tree k, so the
      # splittability check must look at those rows of X (not the first rows of X)
      Xk <- if (!is.null(X) && !is.null(W)) X[which(W[, k] == 1), , drop = FALSE] else X
      if (!is.null(Xk) && nrow(Xk) != length(tr$node_indices)) Xk <- NULL   # fall back to node sizes
      spl <- leaf_splittable(tr, Xk, min_size)          # named by leaf row
      keep[as.integer(names(spl))] <- spl
    }
    list(d_int = depth[!term], d_leaf = depth[term & keep])
  })
  d_int  <- unlist(lapply(per_tree, `[[`, "d_int"))
  d_leaf <- unlist(lapply(per_tree, `[[`, "d_leaf"))
  if (length(d_int) + length(d_leaf) == 0)            # nothing informative (e.g. all trees
    return(list(depth = 0L, n_int = 0L, n_leaf = 0L))   # are unsplittable stumps): avoids max() = -Inf
  D <- max(c(d_int, d_leaf))
  list(depth  = 0:D,
       n_int  = tabulate(d_int  + 1L, nbins = D + 1L),
       n_leaf = tabulate(d_leaf + 1L, nbins = D + 1L))
}

## sum_k log p(T_k | alpha, beta), dropping split-variable/cutpoint terms
log_tree_prior_ab <- function(alpha, beta, cnt) {
  p <- alpha * (1 + cnt$depth)^(-beta)
  sum(cnt$n_int * log(p) + cnt$n_leaf * log1p(-p))
}

## alpha ~ Beta(9.5, 0.5) (mean 0.95), beta ~ Gamma(4, 2) (mean 2)
log_prior_ab_default <- function(alpha, beta)
  dbeta(alpha, 9.5, 0.5, log = TRUE) + dgamma(beta, 4, 2, log = TRUE)

## ---------------------------------------------------------------------
##  Joint random-walk MH on (logit alpha, log beta).
##  Depends on the data only through the tree shapes (not on mu, W, y).
## ---------------------------------------------------------------------
update_alpha_beta <- function(alpha, beta, trees, X = NULL, min_size = 1L,
                              use_splittable = TRUE, W = NULL,
                              Sigma_prop = diag(c(0.3, 0.3)^2),
                              log_prior = log_prior_ab_default) {
  cnt <- depth_counts(trees, X, min_size, use_splittable, W)
  lt <- function(a, b) log_tree_prior_ab(a, b, cnt) + log_prior(a, b) +
    log(a) + log1p(-a) + log(b)      # Jacobian a(1-a)b
  
  uv_p <- c(qlogis(alpha), log(beta)) +
    as.vector(t(chol(Sigma_prop)) %*% rnorm(2))
  a_p <- plogis(uv_p[1]); b_p <- exp(uv_p[2])
  
  accepted <- log(runif(1)) < lt(a_p, b_p) - lt(alpha, beta)
  if (accepted) { alpha <- a_p; beta <- b_p }
  list(alpha = alpha, beta = beta, accepted = accepted, counts = cnt)
}

## Burn-in tuning of Sigma_prop from draws of (logit alpha, log beta)
adapt_cov <- function(uv_draws, eps = 1e-6)
  (2.38^2 / 2) * cov(uv_draws) + eps * diag(2)

## ---------------------------------------------------------------------
##  Usage (current state of the K trees, e.g. inside your sampler loop):
##    ab <- update_alpha_beta(alpha, beta, trees_current, X = X, min_size = 1)
##    alpha <- ab$alpha; beta <- ab$beta
##  For a stored fit, the K trees of iteration it are bart.fit$trees[[it]].
## ---------------------------------------------------------------------
