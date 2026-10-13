# -------------------------------------------------------------------------#
# Helpers shared by bart / cl_bart / bart.ibp / cl_bart.ibp                 #
# -------------------------------------------------------------------------#

# Input checks (NA / non-numeric inputs otherwise surface later as obscure errors)
.check_x = function(z, nm) {
  if (is.data.frame(z)) {
    bad = names(z)[!vapply(z, function(v) is.numeric(v) || is.logical(v), logical(1))]
    if (length(bad)) stop(nm, " has non-numeric columns (", paste(bad, collapse = ", "),
                          "). Convert factors to numeric/dummy columns first, e.g. model.matrix(~ . - 1, data).")
  }
  z = as.matrix(z)
  if (!is.numeric(z) && !is.logical(z)) stop(nm, " must be numeric.")
  storage.mode(z) = "double"
  if (anyNA(z)) stop(nm, " contains ", sum(is.na(z)), " missing values (NA).")
  z
}

# Stored trees: keep only the tree matrix (node_indices has length n and is recomputed
# whenever a stored tree is used for prediction)
.strip_trees = function(trees) lapply(trees, function(t) list(tree_matrix = t$tree_matrix))

.check_y = function(y, x, family) {
  if (is.factor(y)) y = as.integer(y) - 1L
  y = as.numeric(y)
  if (anyNA(y)) stop("y contains ", sum(is.na(y)), " missing values (NA).")
  if (length(y) != nrow(x)) stop("length(y) (", length(y), ") differs from nrow(x) (", nrow(x), ").")
  if (family == "probit") {
    if (!all(y %in% c(0, 1))) stop("For classification y must be coded 0/1 (or a two-level factor).")
  } else if (!is.finite(sd(y)) || sd(y) == 0) stop("y is constant: cannot be standardised.")
  y
}

# Pool the stored draws of several chains into one object
combine_chains = function(fits) {
  out = fits[[1]]
  for (nm in c("trees", "W")) out[[nm]] = do.call(c, lapply(fits, `[[`, nm))
  for (nm in c("sigma2", "ntrees", "sigma2_mu")) out[[nm]] = unlist(lapply(fits, `[[`, nm))
  for (nm in c("y_hat", "IBP.hyper", "y_predicted", "var_count_store", "s", "alpha_beta"))
    if (!is.null(out[[nm]])) out[[nm]] = do.call(rbind, lapply(fits, `[[`, nm))
  out$npost = sum(sapply(fits, `[[`, "npost"))
  out$chain = rep(seq_along(fits), sapply(fits, `[[`, "npost"))
  out
}


##########################################################################
#### Parallel lapply that works on Windows (PSOCK) and Linux/macOS (fork)
##########################################################################
ibp_lapply = function(X, FUN, ..., n_cores = 1, seed = NULL) {
  n_cores = min(n_cores, length(X))
  .args = list(...)
  .wrap = function(el, .fun, .args) do.call(.fun, c(list(el), .args))   # avoids name clashes with parLapply's own arguments
  if (n_cores <= 1) {
    if (!is.null(seed)) set.seed(seed)
    return(lapply(X, FUN, ...))
  }
  if (.Platform$OS.type == "windows" || isTRUE(getOption("bart.ibp.psock"))) {
    cl = parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    env = environment(ibp_lapply)
    fns = Filter(function(nm) is.function(get(nm, envir = env)), ls(env, all.names = TRUE))
    fns = setdiff(fns, ".ibp_row_kernel_cpp")         # compiled code cannot be shipped
    parallel::clusterExport(cl, fns, envir = env)
    if (exists(".ibp_cpp_source", envir = env)) {    # workers compile their own copy (once)
      parallel::clusterExport(cl, ".ibp_cpp_source", envir = env)
      parallel::clusterEvalQ(cl, try(.ibp_try_compile(), silent = TRUE))
    }
    parallel::clusterSetRNGStream(cl, seed)
    parallel::parLapplyLB(cl, X, .wrap, .fun = FUN, .args = .args)
  } else {
    old_kind = RNGkind()[1]
    RNGkind("L'Ecuyer-CMRG")
    on.exit(RNGkind(old_kind), add = TRUE)
    if (!is.null(seed)) set.seed(seed)
    parallel::mclapply(X, .wrap, .fun = FUN, .args = .args, mc.cores = n_cores, mc.set.seed = TRUE, mc.preschedule = FALSE)
  }
}


