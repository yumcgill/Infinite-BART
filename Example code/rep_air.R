library(parallel)
options(mc.cores = 10)
source("simulate_data.R")
source("myBART.R")
source("myBART-IBP.R")
data(airquality)
airquality = airquality[complete.cases(airquality), ]
names(airquality)
set.seed(12)
# Set design matrix and response vector (full: both train and test)
x.full = model.matrix( ~ Solar.R + Wind + Temp + Month + Day, data=airquality )
y.full = log(airquality$Ozone)
n.full = length(y.full)
n = round(n.full*0.8)

air_rep <- function(seed){
  set.seed(seed)
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
  
  
  # fit BART and BART.IBP
  bart.ibp.fit = bart.ibp(x,y,x.pred = X.pred,eta = -93,delta = 100,gamma = 30,nburn = 1000,npost = 5000,nthin = 1)
  bart.fit = bart(x,y,x.pred = X.pred,nburn = 1000,npost = 5000,nthin = 1)
  
  
  # Compare the predictive values of BART and BART.IBP
  pred_bart = colMeans(bart.fit$y_predicted)
  pred_bart_ibp = colMeans(bart.ibp.fit$y_predicted)
  mse_bart = mean((pred_bart-y_true_pred)^2)
  mse_bart_ibp = mean((pred_bart_ibp-y_true_pred)^2)
  return(list(mse_bart=mse_bart,mse_bart_ibp=mse_bart_ibp))
}

bart_ibp_air<-mclapply(sample(c(1:100000000),10), function(x) air_rep(x))
