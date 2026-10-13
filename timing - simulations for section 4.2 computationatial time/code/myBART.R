#' @export
#' @importFrom stats 'rgamma' 'runif' 'dnorm' 'sd' 'rnorm' 'pnorm' 'qnorm'

### Load Auxiliary functions
### (no external packages needed any more: rdirichlet is in IBP_functions.R and
###  update_z no longer uses truncnorm)
source("aux_tree.R")
source("modify_tree.R")
source("predict_tree.R")
source("parameters_quantities_tree.R")
source("IBP_functions.R")           # rdirichlet (used by update_s)
source("update_hyperpar_bart.R")
source("parallel_utils.R")


# -------------------------------------------------------------------------#
# bart: same algorithm as before, with the changes made for bart.ibp:       #
#  - fast tree routines (fill_tree_details, likelihood, simulate_mu, prior) #
#    from aux_tree.R / parameters_quantities_tree.R                         #
#  - the fit of tree j is read off its node indices (all observations use   #
#    every tree), instead of re-routing x through the tree                  #
#  - input checks, optional progress / printing                             #
#  - predictions computed after the MCMC from the stored trees (parallel    #
#    over draws), or inline as before                                       #
#  - optional sigma2_mu and (alpha, beta) updates as in bart.ibp            #
#  - several chains in parallel: bart.chains(); pooled with combine_chains()#
#  - family = "probit" gives the classification model (used by cl_bart)     #
# -------------------------------------------------------------------------#

bart = function(x,y,x.pred = NULL,sparse = TRUE,ntrees = 200,node_min_size = 5,alpha = 0.95,beta = 2,
                nu = 3,lambda = NULL,mu_mu = 0,sigma2 = 1,sigma2_mu = NULL,nburn = 20000,npost = 1000,nthin = 10,
                penalise_num_cov = TRUE,lambda_cov = 0.4,nu_cov = 2,
                # ---- new arguments (as in bart.ibp) ----
                sigma2_mu_update = FALSE, alpha_beta_update = FALSE,
                pred_mode = c("posthoc", "inline"),
                n_cores = 1, seed = NULL, verbose = FALSE, progress = TRUE,
                family = c("gaussian", "probit"),
                store_node_indices = FALSE) {   # TRUE: keep each stored tree's node_indices (length n: large memory)

  family = match.arg(family)
  pred_mode = match.arg(pred_mode)
  if (!is.null(seed)) set.seed(seed)

  # ---- input checks
  x = .check_x(x, "x")
  if (!is.null(x.pred)) { x.pred = .check_x(x.pred, "x.pred")
    if (ncol(x.pred) != ncol(x)) stop("x.pred must have the same columns as x.") }
  y = .check_y(y, x, family)

  # Extract MCMC details
  TotIter = nburn + npost*nthin # Total of iterations

  # Storage containers
  store_size = npost
  tree_store = vector('list', store_size)
  sigma2_store = rep(NA, store_size)
  sigma2_mu_store = rep(NA, store_size)
  alpha_beta_store = matrix(NA, store_size, 2)
  y_hat_store = matrix(NA, ncol = length(y), nrow = store_size)
  var_count = rep(0, ncol(x))
  var_count_store = matrix(0, ncol = ncol(x), nrow = store_size)
  s_prob_store = matrix(0, ncol = ncol(x), nrow = store_size)
  tree_fits_store = matrix(0, ncol = ntrees, nrow = length(y))
  if (!is.null(x.pred)) y_pred_store = matrix(NA, ncol = nrow(x.pred), nrow = store_size)

  if (family == "gaussian") {
    # Scale the response target variable
    y_mean = mean(y)
    y_sd = sd(y)
    y_scale = (y - y_mean)/y_sd
  } else {
    # Probit (Albert & Chib): work on the latent z, sigma2 fixed at 1
    y_mean = 0; y_sd = 1; sigma2 = 1
    y_scale = ifelse(y == 0, -3, 3)     # initial latent values (as in the original cl_bart)
  }
  n = length(y_scale)
  p = ncol(x)
  s = rep(1/p, p)

  if (is.null(lambda)){
    q = 0.9
    lambda = qgamma(1-q,nu, rate =1/2)*(y_sd^2)/nu
  }
  if (is.null(sigma2_mu)){
    sigma2_mu = if (family == "gaussian") (0.5/(sqrt(ntrees)*2))^2 else (3/(2*sqrt(ntrees)))^2
  }

  # Create a list of trees for the initial stump
  curr_trees = create_stump(num_trees = ntrees, y = y_scale, X = x)

  # Initialise the predicted values to zero
  y_hat = rowSums(tree_fits_store)

  if (progress) pb = utils::txtProgressBar(min = 1, max = TotIter, style = 3, width = 60, title = 'Running rBART...')

  # Start the MCMC iterations loop
  for (i in 1:TotIter) {

    if (progress) utils::setTxtProgressBar(pb, i)

    # If at the right place, store everything
    if((i > nburn) & ((i - nburn) %% nthin) == 0) {
      curr = (i - nburn)/nthin
      tree_store[[curr]] = if (store_node_indices) curr_trees else .strip_trees(curr_trees)
      sigma2_store[curr] = sigma2
      sigma2_mu_store[curr] = sigma2_mu
      alpha_beta_store[curr, ] = c(alpha, beta)
      y_hat_store[curr,] = if (family == "probit") pnorm(y_hat) else y_hat
      var_count_store[curr,] = var_count
      s_prob_store[curr,] = s
      if (!is.null(x.pred) && pred_mode == "inline"){
        f_pred = get_predictions(curr_trees, x.pred, single_tree = FALSE)
        y_pred_store[curr,] = if (family == "probit") pnorm(f_pred) else f_pred*y_sd + y_mean
      }
    }

    # Start looping through trees
    for (j in 1:ntrees) {

      current_partial_residuals = y_scale - y_hat + tree_fits_store[,j]

      # Propose a new tree via grow/change/prune/swap
      type = sample_move(curr_trees[[j]], i, nburn)

      # Generate a new tree based on the current
      new_tree = update_tree(y = y_scale, X = x, type = type, curr_tree = curr_trees[[j]],
                             node_min_size = node_min_size, s = s)

      # log marginalised likelihood + log tree prior, current and proposed tree
      l_old = tree_full_conditional(curr_trees[[j]], current_partial_residuals, sigma2, sigma2_mu) +
        get_tree_prior(curr_trees[[j]], alpha, beta)
      l_new = tree_full_conditional(new_tree, current_partial_residuals, sigma2, sigma2_mu) +
        get_tree_prior(new_tree, alpha, beta)

      if(exp(l_new - l_old) > runif(1)) {
        curr_trees[[j]] = new_tree
      }

      # Update mu whether tree accepted or not
      curr_trees[[j]] = simulate_mu(curr_trees[[j]], current_partial_residuals, sigma2, sigma2_mu)

      # Updating BART predictions. Every observation uses every tree, so the
      # tree's node_indices already give the leaf of each observation: no need
      # to route x through the tree again.
      current_fit = as.numeric(curr_trees[[j]]$tree_matrix[curr_trees[[j]]$node_indices, 'mu'])
      y_hat = y_hat - tree_fits_store[,j] + current_fit
      tree_fits_store[,j] = current_fit

    } # End loop through trees

    if (family == "gaussian") {
      sum_of_squares = sum((y_scale - y_hat)^2)
      # Update sigma2 (variance of the residuals)
      sigma2 = update_sigma2(sum_of_squares, n = length(y_scale), nu, lambda)
    } else {
      # Update z (latent variable); sigma2 stays 1
      y_scale = update_z(y, y_hat)
    }

    if (sigma2_mu_update){
      mu_list = unlist(lapply(curr_trees, function(t) t$tree_matrix[, 'mu']))
      mu_list = mu_list[!is.na(mu_list)]
      sigma2_mu = update_sigma2_mu(sqrt(sigma2_mu), mu_list)[[1]]^2   # works on the standard deviation
    }
    if (alpha_beta_update){
      ab = update_alpha_beta(alpha, beta, curr_trees, X = x)          # all observations use every tree
      alpha = ab$alpha; beta = ab$beta
    }

    #  TO DO !!!
    # var_count =

    # Update s = (s_1, ..., s_p), where s_p is the probability that predictor p is used to create new terminal nodes
    if (sparse == 'TRUE' & i > floor(TotIter*0.1)){
      s = update_s(var_count, p, 1)
    }

    if (verbose) print(i)
  } # End iterations loop

  if (progress) { close(pb); cat('\n') }

  out = list(trees = tree_store,
             sigma2 = sigma2_store*y_sd^2,
             y_hat = y_hat_store*y_sd + y_mean,
             npost = npost,
             nburn = nburn,
             nthin = nthin,
             x = x,
             ntrees = ntrees,
             y_mean = y_mean,
             y_sd = y_sd,
             var_count_store = var_count_store,
             s = s_prob_store,
             sigma2_mu = sigma2_mu_store,
             alpha_beta = alpha_beta_store,
             y_predicted = if (is.null(x.pred)) 0 else y_pred_store,
             family = family)

  if (!is.null(x.pred) && pred_mode == "posthoc")
    out$y_predicted = predict_bart(out, x.pred, n_cores = n_cores)
  out

} # End main function


# -------------------------------------------------------------------------#
# Predictions from the stored trees (deterministic given the draws, so the  #
# parallel version gives exactly the same numbers as the serial one)        #
# -------------------------------------------------------------------------#
.predict_bart_chunk = function(trees_chunk, x.pred) {
  t(vapply(trees_chunk, function(tr) get_predictions(tr, x.pred, single_tree = FALSE), numeric(nrow(x.pred))))
}

predict_bart = function(object, x.pred, n_cores = 1) {
  x.pred = .check_x(x.pred, "x.pred")
  idx = seq_along(object$trees)
  k = max(1, min(n_cores, length(idx)))
  chunks = split(idx, sort(rep_len(seq_len(k), length(idx))))
  res = ibp_lapply(lapply(chunks, function(ii) object$trees[ii]), .predict_bart_chunk,
                   n_cores = n_cores, x.pred = x.pred)
  f = do.call(rbind, res)
  if (identical(object$family, "probit")) pnorm(f) else f * object$y_sd + object$y_mean
}


# Several independent chains in parallel (pool them with combine_chains)
bart.chains = function(x, y, ..., n_chains = 4, n_cores = n_chains, seed = NULL) {
  args = list(...)
  args$progress = FALSE
  run_chain = function(ch, x, y, args) do.call(bart, c(list(x = x, y = y), args))
  ibp_lapply(seq_len(n_chains), run_chain, n_cores = n_cores, seed = seed, x = x, y = y, args = args)
}


# -------------------------------------------------------------------------#
# cl_bart: probit BART (Albert & Chib). Same defaults as before.            #
# CHANGES w.r.t. the original cl_bart:                                      #
#  - sigma2 is FIXED at 1. The original updated sigma2 from                 #
#    sum((y - pnorm(y_hat))^2), which is not part of the probit model (the  #
#    latent z has unit variance by construction) and distorts the fit.     #
#  - sigma2_mu = NULL (default) uses Chipman et al.'s probit choice          #
#    (3 / (2 sqrt(ntrees)))^2; pass sigma2_mu = 1 for the old default.      #
#  - x.pred is supported (the original returned an undefined y_pred_store). #
#  - y_hat / y_predicted are P(y = 1 | x) for each stored draw.             #
#  - lambda, nu, sigma2 are accepted for compatibility but not used.       #
# -------------------------------------------------------------------------#
cl_bart = function(x, y, x.pred = NULL, sparse = TRUE, ntrees = 10, node_min_size = 5, alpha = 0.95, beta = 2,
                   nu = 3, lambda = 0.1, mu_mu = 0, sigma2 = 1, sigma2_mu = NULL,
                   nburn = 1000, npost = 1000, nthin = 1,
                   penalise_num_cov = TRUE, lambda_cov = 0.4, nu_cov = 2, ...) {
  bart(x, y, x.pred = x.pred, sparse = sparse, ntrees = ntrees, node_min_size = node_min_size,
       alpha = alpha, beta = beta, mu_mu = mu_mu, sigma2_mu = sigma2_mu,
       nburn = nburn, npost = npost, nthin = nthin, family = "probit", ...)
}

cl_bart.chains = function(x, y, ..., n_chains = 4, n_cores = n_chains, seed = NULL) {
  args = utils::modifyList(list(ntrees = 10, nburn = 1000, npost = 1000, nthin = 1), list(...))
  do.call(bart.chains, c(list(x = x, y = y), args,
                         list(family = "probit", n_chains = n_chains, n_cores = n_cores, seed = seed)))
}
