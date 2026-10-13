###########################################################################
## Timing study: computational cost of BART and BART-IBP (MCMC fit only)
## for different sample sizes n and numbers of covariates p.
##
## Data: Friedman function, y = 10 sin(pi x1 x2) + 20 (x3 - 0.5)^2 + 10 x4 + 5 x5 + N(0,1),
##       x ~ U(0,1)^p.
## For each (n, p) and replicate, both samplers are run for the same number of
## iterations on the same data, and the elapsed time of each call is recorded.
##
## Usage: setwd("<this folder>"); source("timing_study.R")
## Results are appended to results/timings.csv as each run finishes (re-running
## skips the runs already done); the LaTeX table is written to latex/timing_table.tex.
###########################################################################
TIM_DIR = getwd()

cfg = list(
  n_grid = c(100, 500, 1000, 5000, 10000),
  p_grid = c(10, 50, 100, 200),
  reps   = 1,                 # replicated data sets per (n, p)
  nburn  = 100,               # iterations: nburn + npost * nthin
  npost  = 10,
  nthin  = 10,                # only npost draws are stored (keeps memory low for large n)
  bart_ntrees = 200,
  ibp_init = list(trees_per_obs = 40, rho = 0.5, total_trees = 100),
  calibrate = TRUE,           # prior calibration of Chipman et al. (2010) on the standardised scale
  sparse = FALSE,
  seed = 2026,
  code_dir = file.path(TIM_DIR, "code"),
  results_file = file.path(TIM_DIR, "results", "timings.csv"),
  table_file = file.path(TIM_DIR, "latex", "timing_table.tex")
)

## ---- load the BART / BART-IBP code
local({ old = setwd(cfg$code_dir); on.exit(setwd(old))
  suppressMessages({ source("myBART-IBP.R"); source("myBART.R") }) })

friedman = function(x) 10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-0.5)^2 + 10*x[,4] + 5*x[,5]

calib = function(y, m) {     # lambda: P(sigma < 1) = 0.9 (nu = 3); sigma_mu: k = 2 (standardised scale)
  ys = (y - mean(y)) / sd(y)
  list(lambda = qchisq(0.1, 3) / 3, sigma2_mu = (diff(range(ys)) / 2 / (2 * sqrt(m)))^2)
}

time_one = function(n, p, rep) {
  seed = cfg$seed + 1e5 * match(n, cfg$n_grid) + 100 * p + rep
  set.seed(seed)
  x = matrix(runif(n * p), n); y = rnorm(n, friedman(x))
  common = list(x = x, y = y, nburn = cfg$nburn, npost = cfg$npost, nthin = cfg$nthin,
                sparse = cfg$sparse, progress = FALSE)
  set.seed(seed + 1)
  a = c(common, list(ntrees = cfg$bart_ntrees), if (cfg$calibrate) calib(y, cfg$bart_ntrees))
  t_bart = unname(system.time(fb <- do.call(bart, a))[3])
  rm(fb); invisible(gc())
  set.seed(seed + 2)
  a = c(common, list(ibp_init = cfg$ibp_init), if (cfg$calibrate) calib(y, cfg$ibp_init$trees_per_obs))
  t_ibp = unname(system.time(fi <- do.call(bart.ibp, a))[3])
  K = mean(fi$ntrees)
  rm(fi); invisible(gc())
  data.frame(n = n, p = p, rep = rep, iterations = cfg$nburn + cfg$npost * cfg$nthin,
             time_bart = t_bart, time_bart_ibp = t_ibp, ntrees_bart_ibp = K)
}

## ---- run (cheapest settings first; already-done runs are skipped)
dir.create(dirname(cfg$results_file), showWarnings = FALSE, recursive = TRUE)
done = if (file.exists(cfg$results_file)) read.csv(cfg$results_file) else NULL
grid = expand.grid(rep = seq_len(cfg$reps), p = cfg$p_grid, n = cfg$n_grid)
grid = grid[order(grid$n * (1 + grid$p / 100)), ]
for (g in seq_len(nrow(grid))) {
  n = grid$n[g]; p = grid$p[g]; r = grid$rep[g]
  if (!is.null(done) && any(done$n == n & done$p == p & done$rep == r)) next
  res = time_one(n, p, r)
  write.table(res, cfg$results_file, sep = ",", row.names = FALSE,
              col.names = !file.exists(cfg$results_file), append = file.exists(cfg$results_file))
  cat(sprintf("[%s] n = %5d, p = %3d, rep %d: BART %6.1f s, BART-IBP %6.1f s\n",
              format(Sys.time(), "%H:%M:%S"), n, p, r, res$time_bart, res$time_bart_ibp))
}

## ---- LaTeX table: time per 100 iterations, and ratio
write_timing_table = function(d, file = cfg$table_file) {
  fmt = function(x, k) formatC(x, format = "f", digits = k)
  d$b = 100 * d$time_bart / d$iterations; d$i = 100 * d$time_bart_ibp / d$iterations
  ns = sort(unique(d$n)); ps = sort(unique(d$p))
  L = c(sprintf("\\begin{tabular}{r%s}", paste(rep("ccc", length(ps)), collapse = "")), "\\toprule",
        paste0(" & ", paste(sprintf("\\multicolumn{3}{c}{$p = %d$}", ps), collapse = " & "), " \\\\"),
        paste(sprintf("\\cmidrule(lr){%d-%d}", 3 * seq_along(ps) - 1, 3 * seq_along(ps) + 1), collapse = ""),
        paste0("$n$ & ", paste(rep("BART & BART-IBP & ratio", length(ps)), collapse = " & "), " \\\\"), "\\midrule")
  for (n in ns) {
    cells = sapply(ps, function(p) { s = d[d$n == n & d$p == p, ]
      if (!nrow(s)) return("-- & -- & --")
      sprintf("%s & %s & %s", fmt(mean(s$b), 1), fmt(mean(s$i), 1), fmt(mean(s$i) / mean(s$b), 2)) })
    L = c(L, sprintf("%s & %s \\\\", format(n, big.mark = ","), paste(cells, collapse = " & ")))
  }
  L = c(L, "\\bottomrule", "\\end{tabular}")
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  writeLines(L, file); invisible(L)
}
timings = read.csv(cfg$results_file)
write_timing_table(timings)
print(timings)
