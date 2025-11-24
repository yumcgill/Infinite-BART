library(parallel)
options(mc.cores = 23)
source("simulate_data.R")
source("myBART.R")
source("myBART-IBP.R")
set.seed(2151)
# Number of groups
n.groups = 5
# group size
g.size = 40
# Total sample size
n= n.groups*g.size
# Number of explanatory variables
p=10


# Generate independent explanatory variable on [0,1], with mean at j/(p+1) for j=1,..,p and variance 1
a = (1:p/(p+1))^2
b = a*(((p+1)/(1:p)) -1)

fred_clustered_rep <- function(seed){
  set.seed(seed)
  x = matrix(rep(0,n*p),ncol = p)
  for (j in 1:p){
    x[,j]= rbeta(n,a[j],b[j])
  }
  
  # Simulate response vector
  y = rep(0,n)
  # grid defining the indexes at which a new group starts
  gr=seq(from =1, to=n, by=g.size)
  for (j in 1:length(gr)){
    y[gr[j]:(gr[j]+g.size-1)] = rnorm(g.size, (10*sin(pi*x[gr[j]:(gr[j]+g.size-1),j]*x[gr[j]:(gr[j]+g.size-1),j+1])+ 20*(x[gr[j]:(gr[j]+g.size-1),j+2]-0.5)^(2) +10*x[gr[j]:(gr[j]+g.size-1),j+3]+5*x[gr[j]:(gr[j]+g.size-1),j+4]) , 1)
  }
  
  # Design matrix and response vector
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
  bart.ibp.fit = bart.ibp(x,y,x.pred = X.pred, nburn = 1000,npost = 5000,nthin = 1,eta = -50,delta = 51,gamma = 10)
  # Fit Bart
  bart.fit = bart(x,y,X.pred, nburn = 1000,npost = 5000,nthin = 1)
  
  # compare predictions BART and BART.IBP
  pred_bart = colMeans(bart.fit$y_predicted)
  pred_bart_ibp = colMeans(bart.ibp.fit$y_predicted)
  mse_bart = mean((pred_bart-y_true_pred)^2)
  mse_bart_ibp = mean((pred_bart_ibp-y_true_pred)^2)
  return(list(mse_bart=mse_bart,mse_bart_ibp=mse_bart_ibp))
}


bart_ibp_fred_clustered<-mclapply(sample(c(1:100000000),10), function(x) fred_clustered_rep(x))
fred_clustered_bart_mse<-unlist(mclapply(bart_ibp_fred_clustered, '[[', "mse_bart"))
fred_clustered_bart_ibp_mse<-unlist(mclapply(bart_ibp_fred_clustered, '[[', "mse_bart_ibp"))
saveRDS(list(fred_clustered_bart_mse=fred_clustered_bart_mse,fred_clustered_bart_ibp_mse=fred_clustered_bart_ibp_mse), "mse_fred_clustered2.rds")
#saveRDS(bart_ibp_fred_clustered, "bart_ibp_fred_clustered.rds")
