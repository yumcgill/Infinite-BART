
rm(list = ls())
set.seed(1000)

source("myBART-IBP.R")      # compiles the optional C++ kernel once per session (needs Rcpp + Rtools on Windows)
source("summary_results.R")

## Friedman data
n.full = 300; p = 30
x.full = matrix(runif(n.full*p), ncol = p)
y.full = rnorm(n.full, 10*sin(pi*x.full[,1]*x.full[,2]) + 20*(x.full[,3]-0.5)^2 + 10*x.full[,4] + 5*x.full[,5], 1)
train = sample(n.full, 240)
x = x.full[train, ];  y = y.full[train]
X.pred = x.full[-train, ]; y_true_pred = y.full[-train]

## 1) One chain. Predictions are computed after the MCMC from the stored draws,
##    in parallel over draws (pred_mode = "inline" reproduces the old behaviour).
fit = bart.ibp(x, y, x.pred = X.pred, nburn = 2000, npost = 500, nthin = 2,
               n_cores = 4, seed = 1)
mean((colMeans(fit$y_predicted) - y_true_pred)^2)

## 2) Several independent chains in parallel (the most effective use of cores),
##    then pool the draws. Also useful for convergence checks (e.g. coda::gelman.diag).
fits = bart.ibp.chains(x, y, x.pred = X.pred, nburn = 2000, npost = 250, nthin = 2,
                       n_chains = 4, n_cores = 4, seed = 1)
fit.all = combine_chains(fits)
mean((colMeans(fit.all$y_predicted) - y_true_pred)^2)

## 3) Predictions for new points from an existing fit
y.new = predict_bart_ibp(fit, X.pred, n_cores = 4, seed = 2)

## 4) Calibration on the standardised scale (see notes: NOT the default, it changes the model)
K0 = 100   # roughly the initial number of trees (gamma)
fit.cal = bart.ibp(x, y, x.pred = X.pred, nburn = 2000, npost = 500, nthin = 2,
                   lambda    = qgamma(0.1, 3, rate = 1/2) / 3,                # sigma_hat = 1 on the standardised scale
                   sigma2_mu = (diff(range(scale(y))) / 2 / (2 * sqrt(K0)))^2, # Chipman et al. k = 2 on this scale
                   sparse    = FALSE,                                         # until var_count is implemented
                   n_cores = 4, seed = 1)
mean((colMeans(fit.cal$y_predicted) - y_true_pred)^2)
