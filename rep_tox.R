library(parallel)
options(mc.cores = 10)
source("simulate_data.R")
source("myBART.R")
source("myBART-IBP.R")
aquatic_toxicity <- read.csv("datasets/toxicity/qsar_aquatic_toxicity.csv", header=FALSE, sep=";")


# Set design matrix and response vector (full: both train and test)
x.full = model.matrix( ~ V1 + V2 + V3 + V4 + V5 + V6 + as.factor(V7) + V8, data=aquatic_toxicity )
y.full = log(aquatic_toxicity$V9)
n.full = length(y.full)
n = round(n.full*0.8)

tox_rep <- function(seed){
  set.seed(seed)
# split in training and test set (random sampling)
training.points = sample(1:n.full,n,replace = FALSE)
test.points = 1:n.full
test.points=test.points[-training.points]
# train set
y = y.full[training.points]
x = x.full[training.points,]
# test set
X.pred = x.full[test.points, ]
y_true_pred=y.full[test.points]

# Fit Bart and bart.ibp
bart.ibp.fit = bart.ibp(x,y,x.pred = X.pred,nburn = 1,npost = 500,nthin = 1,eta=-0.6, delta=1.34, gamma=10.36)
bart.fit = bart(x,y,X.pred, nburn = 1000,npost = 5000,nthin = 1)


x_n = x.full[training.points,-c(8)]
# test set
X.pred_n = x.full[test.points,-c(8) ]

bart.ibp.fit2 = bart.ibp(x_n,y,x.pred = X.pred_n, nburn = 1000,npost = 5000,nthin = 1,eta=-0.6, delta=1.34, gamma=10.36)
bart.fit2 = bart(x_n,y,x.pred = X.pred_n,nburn = 1000,npost = 5000,nthin = 1)


# Compare the predictived values of BART and BART.IBP
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

bart_ibp_tox<-mclapply(sample(c(1:100000000),10), function(x) tox_rep(x))

