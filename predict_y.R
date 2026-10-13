
# -------------------------------------------------------------------------#
# Gibbs sweep over rows of W, and the out-of-sample prediction step.       #
#                                                                          #
# KEY OPTIMISATION. The original code called                               #
#     get_predictions_ibp(curr_trees[[k]], x, 1, single_tree = TRUE)[w]    #
# inside the double loop over rows w and trees k, i.e. it routed ALL n     #
# observations down tree k just to read off ONE value: O(n^2 K) work per   #
# MCMC iteration. The trees do not change during the sweep, so here the    #
# unmasked predictions F[i,k] = f_k(x_i) are computed once (one tree       #
# traversal per tree) and then only looked up. The fit of row w without    #
# tree k is also kept as a running sum instead of a K-length dot product.  #
# -------------------------------------------------------------------------#

# A stump with a given leaf value (same structure as create_stump_ibp)
make_stump = function(mu, size = 1) {
  tm = matrix(c(1, NA, NA, NA, NA, NA, mu, size), nrow = 1,
              dimnames = list(NULL, c('terminal','child_left','child_right','parent',
                                      'split_variable','split_value','mu','node_size')))
  list(tree_matrix = tm, node_indices = rep(1, size))
}

# Optional compiled kernel for the inner loop over trees (set in myBART-IBP.R).
# Must return the updated 0/1 row of W. Pure-R fallback below.
.ibp_row_kernel_R = function(Ww, Fw, TFw, a, lu, yw, inv2s) {
  fw = sum(TFw)
  for (k in seq_along(Ww)) {
    fk = TFw[k]
    r0 = yw - (fw - fk)          # residual of y_w without tree k
    r1 = r0 - Fw[k]              # residual of y_w with tree k
    z  = as.numeric(lu[k] < a[k] - (r1 * r1 - r0 * r0) * inv2s)
    if (z != Ww[k]) {
      Ww[k] = z
      nf = z * Fw[k]
      fw = fw - fk + nf
      TFw[k] = nf
    }
  }
  Ww
}


# Sweep the rows `rows` of W (global row indices).
#   Fr, TFr : one row per element of `rows`; unmasked predictions f_k(x_w) and
#             masked ones (TFr = Fr * W[rows, ]) for every current tree k
#   yv      : the "responses" for those rows
#   n_prior : the n used in the IBP predictive weights (n for training, n + n.pred in prediction)
#   guard_all_empty : training behaviour of the original code (skip the update of
#             row w if every column would be empty without w)
ibp_sweep_rows = function(W, Fr, TFr, trees, rows, yv, n_prior, eta, delta, sigma2,
                          K.new.max, K.new.par, mu_mu, sigma2_mu, guard_all_empty = TRUE,
                          row_kernel = .ibp_row_kernel_R) {

  counts_all = colSums(W)
  inv2s = 1 / (2 * sigma2)
  sdev  = sqrt(sigma2)
  log_pois = dpois(0:K.new.max, K.new.par, log = TRUE)

  for (l in seq_along(rows)) {
    w  = rows[l]
    yw = yv[l]
    Ww = W[w, ]
    counts = counts_all - Ww
    empty = which(counts == 0)
    K = length(Ww)

    do_update = if (guard_all_empty) length(empty) < K else TRUE
    if (length(empty) == K) do_update = FALSE        # would leave zero columns

    if (do_update) {
      if (length(empty) > 0) {   # trees used only by observation w are removed
        W   = W[, -empty, drop = FALSE]
        Fr  = Fr[, -empty, drop = FALSE]
        TFr = TFr[, -empty, drop = FALSE]
        trees = trees[-empty]
        counts = counts[-empty]
        Ww = Ww[-empty]
        K = length(Ww)
      }
      # log prior odds  (counts_k - eta) / (n - 1 + delta + eta - counts_k)
      a  = log(counts - eta) - log(n_prior - 1 + delta + eta - counts)
      lu = qlogis(runif(K))      # logit of uniforms: z = 1  <=>  lu < log-odds
      Ww = row_kernel(Ww, Fr[l, ], TFr[l, ], a, lu, yw, inv2s)
      W[w, ]   = Ww
      TFr[l, ] = Ww * Fr[l, ]
      counts_all = counts + Ww
    }

    #### Observation w selects new trees (computed on the log scale)
    new.mu = rnorm(K.new.max, mu_mu, sigma2_mu)
    fw = sum(TFr[l, ])
    lp = log_pois + dnorm(yw, fw + c(0, cumsum(new.mu)), sdev, log = TRUE)
    if (anyNA(lp)) {
      stop(sprintf(paste0("Non-finite quantity in the new-tree step (row %d): y_w = %s, fit = %s, sigma2 = %s, ",
                          "sigma2_mu = %s, new mu = [%s], K.new.par = %s, eta = %s, delta = %s. ",
                          "NA fitted values in %d trees."),
                   w, format(yw), format(fw), format(sigma2), format(sigma2_mu),
                   paste(format(new.mu, digits = 3), collapse = ", "), format(K.new.par),
                   format(eta), format(delta), sum(colSums(!is.finite(Fr)) > 0)), call. = FALSE)
    }
    if (all(lp == -Inf)) {
      K.new = 1        # as in the original code: no new trees if every probability underflows
    } else {
      K.new = sample.int(K.new.max + 1, 1, prob = exp(lp - max(lp)))
    }

    if (K.new > 1) {
      m  = K.new - 1
      K  = ncol(W)
      nk = K + seq_len(m)
      W   = cbind(W,   matrix(0, nrow(W), m));   W[w, nk] = 1
      Fr  = cbind(Fr,  matrix(rep(new.mu[1:m], each = nrow(Fr)), nrow(Fr), m))
      TFr = cbind(TFr, matrix(0, nrow(TFr), m)); TFr[l, nk] = new.mu[1:m]
      counts_all = c(counts_all, rep(1, m))
      trees = c(trees, lapply(new.mu[1:m], make_stump))
    }
  }
  list(W = W, Fr = Fr, TFr = TFr, trees = trees)
}


# Re-allocate observations to terminal nodes after W changed (training data).
# Same logic as the original code; F (unmasked predictions at x) is recomputed
# only for trees whose structure / leaf values actually changed.
ibp_refresh_trees = function(W, Fx, TF, trees, x, y, reset_singleton_mu = TRUE) {
  K = ncol(W)
  keep = rep(TRUE, K)
  for (p in seq_len(K)) {
    idx = which(W[, p] == 1)
    sw  = length(idx)
    if (sw == 0) { keep[p] = FALSE; next }     # (column emptied by the last rows: drop)
    tr  = trees[[p]]
    nr  = nrow(tr$tree_matrix)
    changed = FALSE
    if (sw > 1 && nr > 1) {
      tr = fill_tree_details(tr, x[idx, , drop = FALSE])
    } else if (sw > 1) {
      tr$node_indices = rep(1, sw)
      tr$tree_matrix[1, 'node_size'] = sw
    } else {                                   # sw == 1
      # original behaviour: the tree is reset to a stump with mu = 0
      tr = make_stump(if (reset_singleton_mu) 0 else Fx[idx, p])
      changed = TRUE
    }
    while (any(tr$tree_matrix[, 'node_size'] == 0)) {
      tr = clean_prune_tree(x[idx, , drop = FALSE], y[idx], tr)
      changed = TRUE
    }
    if (changed) {
      Fx[, p] = predict_tree_full(tr, x)
      TF[, p] = Fx[, p] * W[, p]
    }
    trees[[p]] = tr
  }
  if (!all(keep)) {
    W = W[, keep, drop = FALSE]; Fx = Fx[, keep, drop = FALSE]
    TF = TF[, keep, drop = FALSE]; trees = trees[keep]
  }
  list(W = W, Fx = Fx, TF = TF, trees = trees)
}


# One prediction draw (the logic of the original prediction_of_y).
# BUG FIX: the original re-computed y_hat_pred from a stale W.pred whose
# columns were no longer aligned with curr_trees once a column had been
# removed. Here y_hat_pred = rowSums of the (always aligned) fitted matrix.
predict_draw = function(trees, W, x, x.pred, eta, delta, gamma, sigma2, sigma2_mu, mu_mu,
                        K.new.max, M.pred, row_kernel = .ibp_row_kernel_R) {

  W = as.matrix(W)
  n = nrow(W); K = ncol(W); n.pred = nrow(x.pred)
  X.all = rbind(x, x.pred)
  pr = n + seq_len(n.pred)                     # rows of the test points in W / X.all

  # Initialize W.pred from the IBP predictive (same random-number order as the original)
  prob.pred = (colSums(W) - eta) / (n + delta)
  W.pred = matrix(rbinom(n.pred * K, 1, rep(prob.pred, n.pred)), n.pred, K, byrow = TRUE)
  W = rbind(W, W.pred)

  # route every point of X.all once per tree (stored in the tree while predicting)
  Fp = matrix(0, n.pred, K)
  for (k in seq_len(K)) {
    tm = trees[[k]]$tree_matrix
    if (nrow(tm) > 1) {
      ni = fill_tree_details(trees[[k]], X.all)$node_indices
      trees[[k]]$ni_all = ni
      Fp[, k] = tm[ni[pr], 'mu']
    } else Fp[, k] = tm[1, 'mu']
  }
  TFp = Fp * W.pred
  y_hat_pred = rowSums(TFp)

  K.new.par = exp( log(gamma) + lgamma(1+delta) + lgamma(n-1+delta+eta ) - lgamma(n+delta) - lgamma(delta+eta) )

  for (mm in seq_len(M.pred)) {

    s = ibp_sweep_rows(W, Fp, TFp, trees, pr, y_hat_pred, n + n.pred, eta, delta, sigma2,
                       K.new.max, K.new.par, mu_mu, sigma2_mu, guard_all_empty = FALSE,
                       row_kernel = row_kernel)
    W = s$W; Fp = s$Fr; TFp = s$TFr; trees = s$trees

    # refresh: only pruning of empty leaves and the singleton reset change predictions
    keep = rep(TRUE, ncol(W))
    for (p in seq_len(ncol(W))) {
      idx = which(W[, p] == 1); sw = length(idx)
      if (sw == 0) { keep[p] = FALSE; next }
      tr = trees[[p]]; tm = tr$tree_matrix
      if (sw == 1) {
        trees[[p]] = make_stump(0); Fp[, p] = 0; TFp[, p] = 0; next
      }
      if (nrow(tm) == 1) next
      leaf_counts = tabulate(tr$ni_all[idx], nbins = nrow(tm))
      if (any(leaf_counts[tm[, 'terminal'] == 1] == 0)) {
        tr = fill_tree_details(tr, X.all[idx, , drop = FALSE])
        while (any(tr$tree_matrix[, 'node_size'] == 0)) {
          tr = clean_prune_tree(X.all[idx, , drop = FALSE], NULL, tr)
        }
        if (nrow(tr$tree_matrix) > 1) {
          tr$ni_all = fill_tree_details(tr, X.all)$node_indices
          Fp[, p] = tr$tree_matrix[tr$ni_all[pr], 'mu']
        } else Fp[, p] = tr$tree_matrix[1, 'mu']
        TFp[, p] = Fp[, p] * W[pr, p]
        trees[[p]] = tr
      }
    }
    if (!all(keep)) {
      W = W[, keep, drop = FALSE]; Fp = Fp[, keep, drop = FALSE]
      TFp = TFp[, keep, drop = FALSE]; trees = trees[keep]
    }
    y_hat_pred = rowSums(TFp)
  }
  y_hat_pred
}


# Backward-compatible wrapper with the original signature
prediction_of_y = function(y_hat_pred, W,  curr_trees, X.pred, X.all,  eta,delta,
                           gamma, M.pred ,tree_fits_store, sigma2,K.new.max,mu_mu,sigma2_mu){
  n = nrow(as.matrix(W))
  predict_draw(curr_trees, W, X.all[seq_len(n), , drop = FALSE], X.pred, eta, delta, gamma,
               sigma2, sigma2_mu, mu_mu, K.new.max, M.pred)
}
