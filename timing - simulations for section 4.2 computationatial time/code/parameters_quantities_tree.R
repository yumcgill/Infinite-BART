# -------------------------------------------------------------------------#
# Description: this script contains 2 functions that are used to generate  #
#              the predictions, update variance and compute the tree prior #
#              and the marginalised likelihood                             #
# -------------------------------------------------------------------------#

# 1. simulate_mu: generate the predicted values (mu's)
# 2. update_sigma2: updates the parameters sigma2
# 3. get_tree_prior: returns the tree log prior score
# 4. tree_full_conditional: computes the marginalised likelihood for all nodes for a given tree
# 5. get_number_distinct_cov: counts the number of distinct covariates that are used in a tree to create the splitting rules
# Compute the full conditionals -------------------------------------------------

# Helper: sum of R within each terminal node, in the order of `leaves`.
# rowsum() is compiled code; the original aggregate() goes through
# data.frames and was ~50x slower. Also robust to empty leaves (the
# original silently mis-aligned nj and S_j if a leaf was empty).
sum_by_leaf = function(R, node_indices, leaves) {
  s = rowsum(R, node_indices, reorder = FALSE)
  out = numeric(length(leaves))
  out[match(as.integer(rownames(s)), leaves)] = s[, 1]
  out
}

tree_full_conditional = function(tree, R, sigma2, sigma2_mu) {

  # log marginal likelihood of the partial residuals R for an individual tree
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)
  nj = tree$tree_matrix[which_terminal,'node_size']
  S_j = sum_by_leaf(R, tree$node_indices, which_terminal)
  # (the original also computed sum(R^2) per leaf, but it cancels in the MH ratio and was never used)

  log_post = 0.5 * ( sum(log( sigma2 / (nj*sigma2_mu + sigma2))) +
                       sum( (sigma2_mu* S_j^2) / (sigma2 * (nj*sigma2_mu + sigma2))))
  return(log_post)
}


# Simulate_par -------------------------------------------------------------

simulate_mu = function(tree, R, sigma2, sigma2_mu) {

  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)
  nj = tree$tree_matrix[which_terminal,'node_size']
  sumR = sum_by_leaf(R, tree$node_indices, which_terminal)

  mu = rnorm(length(nj) ,
             mean = (sumR / sigma2) / (nj/sigma2 + 1/sigma2_mu),
             sd = sqrt(1/(nj/sigma2 + 1/sigma2_mu)))

  tree$tree_matrix[,'mu'] = NA
  tree$tree_matrix[which_terminal,'mu'] = mu
  return(tree)
}

# Update sigma2 -------------------------------------------------------------

update_sigma2 <- function(S, n, nu, lambda){
  u = 1/rgamma(1, shape = (n + nu)/2, rate = (S + nu*lambda)/2)
  return(u)
}


# Update the latent variable z --------------- (THIS FUNCTION IS NEEDED ONLY FOR CLASSIFICATION - ALBERT CHIB AUGMENTATION)

update_z = function(y, prediction){
  # z_i ~ N(prediction_i, 1) truncated to (0, Inf) if y_i = 1 and to (-Inf, 0) if y_i = 0.
  # Inverse-cdf sampling on the log scale: exact, vectorised, stable even when the
  # truncation region is far in the tail, and no dependency on the truncnorm package.
  n = length(y)
  lu = log(runif(n))
  z = numeric(n)
  i1 = y == 1; i0 = !i1
  m1 = prediction[i1]; m0 = prediction[i0]
  # y = 1: x = z - m > -m;  P(X > x) = u P(X > -m)  =>  x = -qnorm(log u + log pnorm(m))
  z[i1] = m1 - qnorm(lu[i1] + pnorm(m1, log.p = TRUE), log.p = TRUE)
  # y = 0: x = z - m < -m;  P(X < x) = u P(X < -m)  =>  x =  qnorm(log u + log pnorm(-m))
  z[i0] = m0 + qnorm(lu[i0] + pnorm(-m0, log.p = TRUE), log.p = TRUE)
  return(z)
}

# Get tree priors ---------------------------------------------------------

get_tree_prior = function(tree, alpha, beta) {

  tm = tree$tree_matrix
  nr = nrow(tm)
  if (nr == 1) return(log(1 - alpha))   # stump

  # depth of each node (parents precede children in the matrix)
  parent = tm[, 'parent']
  level = numeric(nr)
  for (i in 2:nr) level[i] = level[parent[i]] + 1

  term = tm[, 'terminal'] == 1
  sum(log(alpha) - beta * log(1 + level[!term])) +
    sum(log(1 - alpha * (1 + level[term])^(-beta)))
}

# get_num_cov_prior <- function(tree, lambda_cov, nu_cov, penalise_num_cov){
#
#   # Select the rows that correspond to internal nodes
#   which_terminal = which(tree$tree_matrix[,'terminal'] == 0)
#   # Get the covariates that are used to define the splitting rules
#   num_distinct_cov = length(unique(tree$tree_matrix[which_terminal,'split_variable']))
#
#   if (penalise_num_cov == TRUE){
#     if (num_distinct_cov > 0){
#       # A zero-truncated double Poisson
#       log_prior_num_cov = ddoublepois(num_distinct_cov, lambda_cov, nu_cov, log = TRUE) -
#         log(1-ddoublepois(0,lambda_cov, nu_cov))
#     } else {
#       log_prior_num_cov = 0
#     }
#   } else {
#     log_prior_num_cov = 0 # no penalisation
#   }
#
#   return(log_prior_num_cov)
#
# }