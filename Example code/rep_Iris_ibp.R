library(parallel)
options(mc.cores = 10)
source("simulate_data.R")
source("myBART.R")
source("myBART-IBP.R")
set.seed(123)
data("iris")
df <- as.data.frame(iris)
x.full = model.matrix( ~  Species + Sepal.Width + Petal.Length+ Petal.Width, data=df )
y.full = df$Sepal.Length
n.full = length(y.full)

iris_rep_ibp<- function(seed){
  set.seed(seed)
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
  
  
  bart.ibp.fit = bart.ibp(x,y,x.pred = X.pred, nburn = 1000,npost = 5000,nthin = 1,eta = -5,delta = 10,gamma = 15)
  bart.fit = bart(x,y,x.pred = X.pred,nburn = 1000,npost = 5000,nthin = 1)
  
  x_n = x.full[training.points,-c(2:3)]
  # test set
  X.pred_n = x.full[test.points,-c(2:3) ]
  
  bart.ibp.fit2 = bart.ibp(x_n,y,x.pred = X.pred_n, nburn = 1000,npost = 5000,nthin = 1,eta = -5,delta = 10,gamma = 15)
  bart.fit2 = bart(x_n,y,x.pred = X.pred_n,nburn = 1000,npost = 5000,nthin = 1)
  
  pred_bart = colMeans(bart.fit$y_predicted)
  pred_bart_ibp = colMeans(bart.ibp.fit$y_predicted)
  mse_bart = mean((pred_bart-y_true_pred)^2)
  mse_bart_ibp = mean((pred_bart_ibp-y_true_pred)^2)
  
  pred_bart2 = colMeans(bart.fit2$y_predicted)
  pred_bart_ibp2 = colMeans(bart.ibp.fit2$y_predicted)
  mse_bart2 = mean((pred_bart2-y_true_pred)^2)
  mse_bart_ibp2 = mean((pred_bart_ibp2-y_true_pred)^2)
  return(list(mse_bart=mse_bart,mse_bart_ibp=mse_bart_ibp,mse_bart2=mse_bart2,mse_bart_ibp2=mse_bart_ibp2))
  
}

bart_ibp_iris<-mclapply(sample(c(1:100000000),10), function(x) iris_rep_ibp(x))


