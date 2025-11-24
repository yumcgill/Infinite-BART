library(parallel)
options(mc.cores = 23) 
source("simulate_data.R")
source("myBART.R")
source("myBART-IBP.R")
set.seed(231)

bart_causal<-function(n,ss){
  set.seed(ss)
  X1 <- rnorm(n)
  X2 <- rbinom(n, 1, 0.5)
  X3 <- runif(n)
  X4 <- runif(n)
  X5 <- runif(n)
  X <- data.frame(X1, X2, X3 , X4, X5)
  
  # Propensity score (nonlinear dependence)
  pscore <- plogis(0.5 * X1 - 0.7 * X2 + 0.3 * sin(2 * pi * X3))
  T <- rbinom(n, 1, pscore)
  # True treatment effect varies with X1: ATE = 1 
  tau <- 1 + 0.5 * X1  
  
  # Nonlinear outcome model
  Y0 <- 2 + 0.3 * X1^2 - 0.5 * X2 + sin(2 * pi * X3) + rnorm(n)
  Y1 <- Y0 + tau
  
  # Observed outcome
  y <- ifelse(T == 1, Y1, Y0)
  
  
  ## estimated PS
  est.ps  <- glm(T ~ X1 + X2 + X3 + X4 + X5 , family = binomial)$fitted.values
  
  x <- data.frame(T,X,est.ps)
  
  X.pred <- rbind(data.frame(T=1,X,est.ps),data.frame(T=0,X,est.ps))
  
  prova.bart.ibp = bart.ibp(x,y,x.pred = X.pred, nburn = 1000,npost = 5000,nthin = 1)
  ate.post<- apply(prova.bart.ibp$y_predicted,1, function(x) (mean(x[1:n] - x[(n+1):(2*n)])) ) 
  
  
  
  prova.bart = bart(x,y,x.pred = X.pred, nburn = 1000,npost = 5000,nthin = 1)
  ate.post_bart<- apply(prova.bart$y_predicted,1, function(x) (mean(x[1:n] - x[(n+1):(2*n)])) ) 
  
  bart_ibp_mse <- mean((ate.post - 1)^2) 
  bart_mse <- mean((ate.post_bart - 1)^2) 
  return(list(bart_mse=bart_mse,bart_ibp_mse=bart_ibp_mse))
}
bart_ibp_causal<-mclapply(sample(c(1:100000000),10), function(x) bart_causal(100,x))


