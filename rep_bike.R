library(parallel)
options(mc.cores = 23)
source("simulate_data.R")
source("myBART.R")
source("myBART-IBP.R")
motor_data <- read.csv("datasets/bike/bike.csv")

# remove missing values
motor_data = motor_data[complete.cases(motor_data), ]
set.seed(2)


x.full = model.matrix( ~ name + year + km_driven + seller_type + ex_showroom_price + owner, data=motor_data )
y.full = log(motor_data$selling_price)
n.full = length(y.full)
n = round(n.full*0.8)
motor_rep <- function(seed){
  set.seed(seed)
  training.points = sample(1:n.full,n,replace = FALSE)
  test.points = 1:n.full
  test.points=test.points[-training.points]
  # training set
  y = y.full[training.points]
  x = x.full[training.points,]
  # test set
  X.pred = x.full[test.points, ]
  y_true_pred=y.full[test.points]

  # fit BART and BART.IBP
  bart.ibp.fit = bart.ibp(x,y,x.pred = X.pred,nburn = 1000,npost = 5000,nthin = 1)
  bart.fit = bart(x,y,x.pred = X.pred,nburn = 1000,npost = 5000,nthin = 1)


  # Compare the predictived values of BART and BART.IBP
  pred_bart = colMeans(bart.fit$y_predicted)
  pred_bart_ibp = colMeans(bart.ibp.fit$y_predicted)
  mse_bart = mean((pred_bart-y_true_pred)^2)
  mse_bart_ibp = mean((pred_bart_ibp-y_true_pred)^2)
  return(list(mse_bart=mse_bart,mse_bart_ibp=mse_bart_ibp,bart=bart.fit,bart_ibp=bart.ibp.fit))
}

bart_ibp_motor<-mclapply(sample(c(1:100000000),10), function(x) motor_rep(x))
