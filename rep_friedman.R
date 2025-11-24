library(parallel)
options(mc.cores = 23)
set.seed(1000)
source("simulate_data.R")
source("myBART.R")
source("myBART-IBP.R")
set.seed(213)

fred_rep <- function(seed){
  set.seed(seed)
  # Sample size, train set
  n=300
  # Number of expl variables
  p=30
  # sample size, test set
  n.pred = 50
  
  # initial value for s vector in sparse BART
  s = rep(1/p, p)
  
  
  #####################
  #### Simulate from the plain Friedman example
  #####################
  
  # design matrix and response vectors (full: training + test)
  x = matrix(runif(n*p,0,1),ncol = p)
  y = rnorm(n, mean = (10*sin(pi * x[,1]*x[,2])+ 20*(x[,3]-0.5)^(2) +10*x[,4]+5*x[,5]), 1)
  x.full = x
  y.full = y
  n.full = length(y.full)
  n = round(n.full*0.8)
  
  # split in training and test set
  training.points = sample(1:n.full,n,replace = FALSE)
  test.points = 1:n.full
  test.points=test.points[-training.points]
  # train set
  y = y.full[training.points]
  x = x.full[training.points,]
  # test set
  X.pred = x.full[test.points, ]
  y_true_pred=y.full[test.points]
  
  # Fit BART.IBP
  bart.ibp.fit = bart.ibp(x,y,x.pred = X.pred, nburn = 1000,npost = 5000,nthin = 1,eta = -50,delta = 52,gamma = 8)
  
  # Variable importance
  vi = variable_importance(bart.ibp.fit)
  ntree = mean(bart.ibp.fit$ntrees)
  # In sample estimates of f(x) with credible regions
  in_sample_poster_summary = posterior_f_x(bart.ibp.fit)
  # fit also BART
  bart.fit = bart(x,y,x.pred = X.pred, nburn = 1000,npost = 5000,nthin = 1)
  
  # compare predictions BART and BART.IBP
  pred_bart = colMeans(bart.fit$y_predicted)
  pred_bart_ibp = colMeans(bart.ibp.fit$y_predicted)
  mse_bart = mean((pred_bart-y_true_pred)^2)
  mse_bart_ibp = mean((pred_bart_ibp-y_true_pred)^2)
  
  return(list(mse_bart=mse_bart,mse_bart_ibp=mse_bart_ibp,in_sample_poster_summary=in_sample_poster_summary,vi= vi, ntree= ntree))
}


bart_ibp_fred<-mclapply(sample(c(1:100000000),10), function(x) fred_rep(x))

