set.seed(1000)
library(parallel)
options(mc.cores = 23) 
library(MASS)
data(Boston)
names(Boston)
source("simulate_data.R")
source("myBART.R")
source("myBART-IBP.R")
set.seed(123)

# Set design matrix and response vector (full: both train and test)
x.full = model.matrix( ~ crim + zn + indus + chas + nox + rm + age + dis + rad + tax + ptratio + black +lstat , data=Boston )
y.full = log(Boston$medv)
n.full = length(y.full)
n = round(n.full*0.8)

Boston_rep <- function(seed){
  set.seed(seed)
  # split in training and test set
  training.points = sample(1:n.full,n,replace = FALSE)
  test.points = 1:n.full
  test.points=test.points[-training.points]
  # Train set
  y = y.full[training.points]
  x = x.full[training.points,]
  # Test set
  X.pred = x.full[test.points, ]
  y_true_pred=y.full[test.points]
  
  # fit BART and BART.IBP
  bart.ibp.fit = bart.ibp(x,y,x.pred = X.pred,eta = -0.54,delta = 0.66,gamma = 16.89,nburn = 1000,npost = 5000,nthin = 1)
  bart.fit = bart(x,y,x.pred = X.pred,nburn = 1000,npost = 5000,nthin = 1)
  
  
  # Compare the predictived values of BART and BART.IBP
  pred_bart = colMeans(bart.fit$y_predicted)
  pred_bart_ibp = colMeans(bart.ibp.fit$y_predicted)
  mse_bart = mean((pred_bart-y_true_pred)^2)
  mse_bart_ibp = mean((pred_bart_ibp-y_true_pred)^2)
  return(list(mse_bart=mse_bart,mse_bart_ibp=mse_bart_ibp))
}

bart_ibp_boston<-mclapply(sample(c(1:100000000),10), function(x) Boston_rep(x))
