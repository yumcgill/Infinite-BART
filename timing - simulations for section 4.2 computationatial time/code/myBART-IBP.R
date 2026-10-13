

### Load Auxiliary functions and packages
source("aux_tree.R")
source("modify_tree.R")
source("predict_tree.R")
source("parameters_quantities_tree.R")
source("IBP_functions.R")
source("simulate_data.R")
source("predict_y.R")
source("update_hyperpar_bart.R")
source("IBP_hyper_update.R")
source("parallel_utils.R")

bart.ibp = function(x,y,x.pred = NULL,sparse = TRUE,node_min_size = 5,alpha = 0.95,beta = 2,
          nu = 3,lambda = NULL,mu_mu = 0, sigma2 = 1,sigma2_mu = NULL,lambda_cov = 0.4,nu_cov = 2,
          nburn = 20000,npost = 1000,nthin = 10,penalise_num_cov = TRUE,
          eta.alpha = 0.05,eta.beta = 0.01,delta.alpha = 0.1,delta.beta = 0.01,gamma.alpha = 0.05,gamma.beta = 0.01, # prior hyperparmaters
          sigma.eta= 1,sigma.delta= 1,K.new.max = 5,M.pred = 10,eta = -100,delta = 101,gamma = 100,# MH variances, MCMC par and initial values
          sigma2_mu_update = FALSE, alpha_beta_update = FALSE,
          # ---- new arguments ----
          pred_mode = c("posthoc", "inline"), # "posthoc": predictions computed after the MCMC from the stored draws
                                              #   (can run in parallel, same distribution); "inline": as in the original code
          n_cores = 1,                        # cores for the post-hoc prediction step
          seed = NULL,                        # optional seed (the post-hoc parallel step uses its own RNG streams)
          reset_singleton_mu = TRUE,          # TRUE = original behaviour (see notes)
          use_cpp = TRUE,                     # use the compiled W-sweep kernel if it could be built
          verbose = FALSE, progress = TRUE,
          ibp_update = c("collapsed", "original"), # "collapsed": new sampler for (eta, delta, gamma), see IBP_hyper_update.R
          ibp_n_rep = 2,                           # slice sweeps for (eta, delta) per MCMC iteration (cheap)
          ibp_init = NULL,      # e.g. list(trees_per_obs = 50, rho = 0.5, total_trees = 200): sets eta, delta, gamma
                                #   via ibp_elicit(); explicitly supplied eta / delta / gamma take precedence
          ibp_prior_cv = 0.5,   # with ibp_init: Gamma priors centred on the elicited values with this CV
                                #   (NULL = keep the priors as given); explicitly supplied prior values take precedence
          family = c("gaussian", "probit"),  # "probit": binary y via Albert-Chib augmentation (use cl_bart.ibp)
          store_node_indices = FALSE         # TRUE: keep each stored tree's node_indices (not needed for prediction
          ) {                                #   or the summaries; for large n they dominate the memory use)
  ibp_update = match.arg(ibp_update)
  family = match.arg(family)

  ## ---- elicited initialisation / priors for the IBP (see IBP_hyper_update.R)
  ibp_elicited = NULL
  if (!is.null(ibp_init)) {
    need = c("trees_per_obs", "rho", "total_trees")
    if (!is.list(ibp_init) || !all(need %in% names(ibp_init)))
      stop("ibp_init must be a list with elements trees_per_obs, rho and total_trees.")
    ibp_elicited = ibp_elicit(nrow(as.matrix(x)), ibp_init$trees_per_obs, ibp_init$rho, ibp_init$total_trees)
    if (missing(eta))   eta   = ibp_elicited[["eta"]]
    if (missing(delta)) delta = ibp_elicited[["delta"]]
    if (missing(gamma)) gamma = ibp_elicited[["gamma"]]
    if (!is.null(ibp_prior_cv)) {
      pri = ibp_prior_from_elicit(ibp_elicited, cv = ibp_prior_cv)
      if (missing(eta.alpha))   eta.alpha   = pri$eta.alpha
      if (missing(eta.beta))    eta.beta    = pri$eta.beta
      if (missing(delta.alpha)) delta.alpha = pri$delta.alpha
      if (missing(delta.beta))  delta.beta  = pri$delta.beta
      if (missing(gamma.alpha)) gamma.alpha = pri$gamma.alpha
      if (missing(gamma.beta))  gamma.beta  = pri$gamma.beta
    }
  }
  if (!(eta < 1 && eta + delta > 0 && gamma > 0))
    stop("Invalid IBP starting values: need eta < 1, eta + delta > 0 and gamma > 0.")
  ibp_prior = list(eta.alpha = eta.alpha, eta.beta = eta.beta, delta.alpha = delta.alpha,
                   delta.beta = delta.beta, gamma.alpha = gamma.alpha, gamma.beta = gamma.beta)
  ibp_state = NULL
  eta0 = eta; delta0 = delta; gamma0 = gamma

  pred_mode = match.arg(pred_mode)
  if (!is.null(seed)) set.seed(seed)
  row_kernel = if (use_cpp && exists(".ibp_row_kernel_cpp", mode = "function")) .ibp_row_kernel_cpp else .ibp_row_kernel_R

  # ---- input checks
  x = .check_x(x, "x")
  if (!is.null(x.pred)) { x.pred = .check_x(x.pred, "x.pred")
    if (ncol(x.pred) != ncol(x)) stop("x.pred must have the same columns as x.") }
  y = .check_y(y, x, family)

  if (family == "gaussian") {
    # Scale the response target variable
    y_mean = mean(y)
    y_sd = sd(y)
    y_scale = (y - y_mean)/y_sd
  } else {
    # Probit (Albert & Chib): y = 1(z > 0), z = sum_k W_k f_k(x) + N(0, 1).
    # The sampler works on the latent z (no scaling) and sigma2 is fixed at 1.
    y_mean = 0; y_sd = 1
    sigma2 = 1
    y_scale = ifelse(y == 0, -3, 3)       # initial latent values, as in cl_bart
  }

  n = length(y_scale)
  p = ncol(x)
  s = rep(1/p, p)

  if (is.null(lambda)){
    q = 0.9
    lambda = qgamma(1-q,nu, rate =1/2)*(y_sd^2)/nu
  }

  # Initialize W (training set n, only) from the IBP prior
  W = as.matrix(prior.IBP(n,0,eta,delta,gamma))
  ntrees = ncol(W)

  if (is.null(sigma2_mu)){
    sigma2_mu = if (family == "gaussian") (0.5/(sqrt(ntrees)*2))^2 else
      (3/(2*sqrt(ntrees)))^2      # Chipman et al. (2010) probit default, k = 2: f(x) in (-3, 3) w.p. 0.95
  }

  # Initialize Trees (stumps with mu = 0)
  curr_trees = create_stump_ibp(num_trees = ntrees, y = y_scale, X = x, W)

  # Fx[i,k] = f_k(x_i) (unmasked), TF = Fx * W (= tree_fits_store of the original code)
  Fx = matrix(0, n, ntrees)
  TF = matrix(0, n, ntrees)
  y_hat = rowSums(TF)

  TotIter = nburn + npost*nthin

  store_size = npost
  tree_store = vector('list', store_size)
  sigma2_store = rep(NA, store_size)
  sigma2_mu_store = rep(NA, store_size)
  y_hat_store = matrix(NA, ncol = length(y_scale), nrow = store_size)
  var_count = rep(0, ncol(x))
  var_count_store = matrix(0, ncol = ncol(x), nrow = store_size)
  s_prob_store = matrix(0, ncol = ncol(x), nrow = store_size)
  IBP.hyper.store= matrix(rep(0,store_size*3),ncol=3)
  IBP.accept.store= matrix(rep(0,TotIter*2),ncol=2)
  ntrees.store = rep(NA, store_size)
  W_store = vector('list',store_size)
  alpha_beta_store = matrix(rep(0,store_size*2),ncol=2)

  accept.eta=0
  accept.delta=0
  if (!is.null(x.pred)){
    n.pred = nrow(x.pred)
    y_pred_store = matrix(0, ncol = n.pred, nrow = store_size)
    y_hat_pred_store = matrix(NA, ncol = n.pred, nrow = store_size)
  } else{
    y_pred_store = 0}

  if (progress) pb = utils::txtProgressBar(min = 1, max = TotIter,style = 3, width = 60, title = 'Running IBP-BART...')

  ####################
  ####### MCMC Sampler
  ####################

  for (i in 1:TotIter) {

    if (progress) utils::setTxtProgressBar(pb, i)

    ### Store (same position in the iteration as the original code)
    if((i > nburn) & ((i - nburn) %% nthin) == 0) {
      curr = (i - nburn)/nthin
      tree_store[[curr]] = if (store_node_indices) curr_trees else .strip_trees(curr_trees)
      W_store[[curr]] = W
      sigma2_store[curr] = sigma2
      y_hat_store[curr,] = if (family == "probit") pnorm(y_hat) else y_hat   # probit: fitted P(y = 1)
      var_count_store[curr,] = var_count
      s_prob_store[curr,] = s
      IBP.hyper.store[curr,]=c(eta,delta,gamma)
      ntrees.store[curr]=ntrees
      sigma2_mu_store[curr]= sigma2_mu
      alpha_beta_store[curr,]=c(alpha,beta)

      if (!is.null(x.pred) && pred_mode == "inline"){
        y_hat_pred = predict_draw(curr_trees, W, x, x.pred, eta, delta, gamma, sigma2, sigma2_mu,
                                  mu_mu, K.new.max, M.pred, row_kernel = row_kernel)
        y_hat_pred_store[curr,] = if (family == "probit") pnorm(y_hat_pred) else y_hat_pred*y_sd + y_mean
      }
    }

    ################
    ### UPDATE THE TREES  (sequential backfitting: see notes on parallelisation)
    ################
    for (j in seq_len(ntrees)) {
      if (length(curr_trees[[j]]$node_indices) > node_min_size) {
        idx = which(W[, j] == 1)
        R = y_scale[idx] - y_hat[idx] + TF[idx, j]      # partial residuals of the obs using tree j
        type = sample_move(curr_trees[[j]], i, nburn)
        new_tree = update_tree(y = y_scale[idx], X = x[idx, , drop = FALSE], type = type,
                               curr_tree = curr_trees[[j]], node_min_size = node_min_size, s = s)
        l_old = tree_full_conditional(curr_trees[[j]], R, sigma2, sigma2_mu) + get_tree_prior(curr_trees[[j]], alpha, beta)
        l_new = tree_full_conditional(new_tree, R, sigma2, sigma2_mu) + get_tree_prior(new_tree, alpha, beta)
        a = exp(l_new - l_old)
        if(a > runif(1)) {
          curr_trees[[j]] = new_tree
        }
        curr_trees[[j]] = simulate_mu(curr_trees[[j]], R, sigma2, sigma2_mu)

        f = predict_tree_full(curr_trees[[j]], x)
        new_fit = f * W[, j]
        y_hat = y_hat - TF[, j] + new_fit
        TF[, j] = new_fit
        Fx[, j] = f
      }
    }

    ##################
    ###  Update the Matrix W
    ##################
    K.new.par = exp( log(gamma) + lgamma(1+delta) + lgamma(n-1+delta+eta ) - lgamma(n+delta) - lgamma(delta+eta) )

    sw = ibp_sweep_rows(W, Fx, TF, curr_trees, seq_len(n), y_scale, n, eta, delta, sigma2,
                        K.new.max, K.new.par, mu_mu, sigma2_mu, guard_all_empty = TRUE,
                        row_kernel = row_kernel)
    rf = ibp_refresh_trees(sw$W, sw$Fr, sw$TFr, sw$trees, x, y, reset_singleton_mu)
    W = rf$W; Fx = rf$Fx; TF = rf$TF; curr_trees = rf$trees
    ntrees = ncol(W)

    ################
    ### Update other BART parameters
    ################
    y_hat = rowSums(TF)
    if (family == "gaussian") {
      sum_of_squares = sum((y_scale - y_hat)^2)
      sigma2 = update_sigma2(sum_of_squares, n = length(y_scale), nu, lambda)
    } else {
      y_scale = update_z(y, y_hat)        # latent z | trees, W, y  (sigma2 stays 1)
    }

    if (sigma2_mu_update){
      mu_list <- unlist(lapply(curr_trees, function(x) x$tree_matrix[, 7]))
      mu_list <- mu_list[!is.na(mu_list)]
      # BUG FIX: update_sigma2_mu samples the standard deviation sigma_mu, not the variance
      sigma2_mu = update_sigma2_mu(sqrt(sigma2_mu), mu_list)[[1]]^2
    }
    if (alpha_beta_update){
      alpha_beta = update_alpha_beta(alpha, beta, curr_trees, X = x, W = W)
      alpha = alpha_beta$alpha
      beta = alpha_beta$beta
    }

    if (sparse == 'TRUE' & i > floor(TotIter*0.1)){
      s = update_s(var_count, p, 1)
    }

    ##############
    ### Update IBP hyperparameters (eta,delta,gamma)
    ##############
    if (ibp_update == "collapsed") {
      # (eta, delta) | W with gamma integrated out, then gamma | eta, delta, W.
      # Slice directions adapt during burn-in only, then are frozen.
      up = update_ibp_hyper(W, eta, delta, gamma, ibp_prior, ibp_state, adapt = (i <= nburn), n_rep = ibp_n_rep)
      eta = up$eta; delta = up$delta; gamma = up$gamma; ibp_state = up$state
      accept.eta = accept.eta + 1; accept.delta = accept.delta + 1   # slice sampling: every step "accepted"
    } else {
      up.etadelta = update_eta_delta_block(W,accept.eta,eta,delta,gamma, sigma.eta,sigma.delta, eta.alpha,eta.beta,delta.alpha,delta.beta)
      eta = up.etadelta[1]
      accept.eta = up.etadelta[3]
      delta = up.etadelta[2]
      accept.delta = up.etadelta[3]
      gamma = update_gamma(W,eta,delta,gamma, gamma.alpha,gamma.beta)
    }

    IBP.accept.store[i,]= c(accept.eta,accept.delta)

    if ( i%%100 == 0 & i < nburn){
      accept_rate.eta = (IBP.accept.store[i,1] - IBP.accept.store[i-99,1])/100
      accept_rate.delta  = (IBP.accept.store[i,2] - IBP.accept.store[i-99,2])/100
      if (accept_rate.eta >0.5){sigma.eta=sigma.eta*1.1}
      if (accept_rate.eta <0.1){sigma.eta=sigma.eta*0.9}
      if (accept_rate.delta >0.5){sigma.delta=sigma.delta*1.1}
      if (accept_rate.delta <0.1){sigma.delta=sigma.delta*0.9}
    }

    if (verbose) print(c(i,ntrees))

  } # End MCMC iterations loop

  if (progress) { close(pb); cat('\n') }

  out = list(trees = tree_store,
             sigma2 = sigma2_store*y_sd^2,
             y_hat = y_hat_store*y_sd + y_mean,
             npost = npost,
             nburn = nburn,
             nthin = nthin,
             ntrees = ntrees.store,
             x = x,
             y_mean = y_mean,
             y_sd = y_sd,
             var_count_store = var_count_store,
             s = s_prob_store,
             IBP.hyper = IBP.hyper.store,
             y_fitted = if (is.null(x.pred)) NULL else y_pred_store,
             y_predicted = if (is.null(x.pred)) NULL else y_hat_pred_store,
             IBP.accept.store = IBP.accept.store,
             W = W_store,
             sigma2_mu = sigma2_mu_store,
             alpha_beta = alpha_beta_store,
             mu_mu = mu_mu, K.new.max = K.new.max, M.pred = M.pred,
             ibp_start = c(eta = eta0, delta = delta0, gamma = gamma0),   # starting values actually used
             ibp_prior = ibp_prior,                                       # priors actually used
             ibp_elicited = ibp_elicited,
             family = family)

  if (!is.null(x.pred) && pred_mode == "posthoc") {
    out$y_predicted = predict_bart_ibp(out, x.pred, n_cores = n_cores, seed = seed, use_cpp = use_cpp)
  }
  out

} # End main function


##########################################################################
#### Post-hoc prediction from the stored draws (parallel over draws)
##########################################################################
# The prediction step never feeds back into the chain (the original code
# discarded the updated W / trees after predicting), so it can be done after
# the MCMC, independently for every stored draw -> embarrassingly parallel.

.predict_one_stored = function(d, x, x.pred, y_sd, y_mean, mu_mu, K.new.max, M.pred, use_cpp, family = "gaussian") {
  rk = if (use_cpp && exists(".ibp_row_kernel_cpp", mode = "function")) .ibp_row_kernel_cpp else .ibp_row_kernel_R
  yp = predict_draw(d$trees, d$W, x, x.pred, d$eta, d$delta, d$gamma, d$sigma2, d$sigma2_mu,
                    mu_mu, K.new.max, M.pred, row_kernel = rk)
  if (family == "probit") pnorm(yp) else yp * y_sd + y_mean
}

predict_bart_ibp = function(object, x.pred, M.pred = object$M.pred, K.new.max = object$K.new.max,
                            n_cores = 1, seed = NULL, use_cpp = TRUE) {
  x.pred = as.matrix(x.pred)
  draws = lapply(seq_len(length(object$trees)), function(t) list(
    trees = object$trees[[t]], W = object$W[[t]],
    eta = object$IBP.hyper[t, 1], delta = object$IBP.hyper[t, 2], gamma = object$IBP.hyper[t, 3],
    sigma2 = object$sigma2[t] / object$y_sd^2,          # back to the scale used inside the sampler
    sigma2_mu = object$sigma2_mu[t]))
  res = ibp_lapply(draws, .predict_one_stored, n_cores = n_cores, seed = seed,
                   x = object$x, x.pred = x.pred, y_sd = object$y_sd, y_mean = object$y_mean,
                   mu_mu = object$mu_mu, K.new.max = K.new.max, M.pred = M.pred, use_cpp = use_cpp,
                   family = if (is.null(object$family)) "gaussian" else object$family)
  do.call(rbind, res)
}


##########################################################################
#### Several independent chains in parallel
##########################################################################
bart.ibp.chains = function(x, y, ..., n_chains = 4, n_cores = n_chains, seed = NULL) {
  args = list(...)
  args$progress = FALSE
  run_chain = function(ch, x, y, args) do.call(bart.ibp, c(list(x = x, y = y), args))
  ibp_lapply(seq_len(n_chains), run_chain, n_cores = n_cores, seed = seed, x = x, y = y, args = args)
}

##########################################################################
#### Classification: probit BART-IBP (Albert & Chib augmentation, as in cl_bart)
##########################################################################
# y must be 0/1 (or a two-level factor). Same arguments as bart.ibp; the
# defaults for the MCMC length follow cl_bart. Output:
#   y_hat       : fitted P(y = 1 | x) for every stored draw (npost x n)
#   y_predicted : predicted P(y = 1 | x.pred) for every stored draw
#   sigma2      : always 1 (fixed by the probit model)
cl_bart.ibp = function(x, y, x.pred = NULL, nburn = 1000, npost = 1000, nthin = 1, ...) {
  bart.ibp(x, y, x.pred = x.pred, nburn = nburn, npost = npost, nthin = nthin, family = "probit", ...)
}

cl_bart.ibp.chains = function(x, y, ..., n_chains = 4, n_cores = n_chains, seed = NULL) {
  bart.ibp.chains(x, y, ..., family = "probit", n_chains = n_chains, n_cores = n_cores, seed = seed)
}


##########################################################################
#### Optional compiled kernel (Rcpp) for the inner loop of the W sweep.
#### Falls back silently to the pure-R version if Rcpp / a compiler is not
#### available (on Windows you need Rtools).
##########################################################################
.ibp_cpp_source = '
#include <Rcpp.h>
// [[Rcpp::export]]
Rcpp::NumericVector ibp_row_kernel_cpp(Rcpp::NumericVector Ww, Rcpp::NumericVector Fw,
    Rcpp::NumericVector TFw, Rcpp::NumericVector a, Rcpp::NumericVector lu,
    double yw, double inv2s) {
  Rcpp::NumericVector W = Rcpp::clone(Ww);
  int K = W.size();
  double fw = 0.0;
  for (int k = 0; k < K; ++k) fw += TFw[k];
  std::vector<double> tf(TFw.begin(), TFw.end());
  for (int k = 0; k < K; ++k) {
    double r0 = yw - (fw - tf[k]);
    double r1 = r0 - Fw[k];
    double z  = (lu[k] < a[k] - (r1 * r1 - r0 * r0) * inv2s) ? 1.0 : 0.0;
    if (z != W[k]) {
      W[k] = z;
      double nf = z * Fw[k];
      fw += nf - tf[k];
      tf[k] = nf;
    }
  }
  return W;
}'

.ibp_try_compile = function() {
  if (!requireNamespace("Rcpp", quietly = TRUE)) return(invisible(FALSE))
  ok = tryCatch({
    e = new.env()
    Rcpp::sourceCpp(code = .ibp_cpp_source, env = e)
    assign(".ibp_row_kernel_cpp", e$ibp_row_kernel_cpp, envir = globalenv())
    TRUE
  }, error = function(err) FALSE)
  invisible(ok)
}
if (!exists(".ibp_row_kernel_cpp", mode = "function")) {
  if (!isTRUE(.ibp_try_compile()))
    message("bart.ibp: Rcpp kernel not available, using the pure-R W sweep (slower but identical).")
}
