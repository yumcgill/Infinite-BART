
### Clean space
rm(list = ls())
### Set seed
set.seed(1000)

### Load Auxiliary functions and packages
source("simulate_data.R")
source("myBART.R")
source("myBART-IBP.R")


###################
#### Simulate the data from the model
###################

# Sample size, train set
n=100
# Number of expl variables
p=15
# sample size, test set
n.pred = 50

# initial value for s vector in sparse BART
s = rep(1/p, p)

### set IBP hyperparameters (to generate simulated data from the model) 
eta.true = -100
delta.true = 101
gamma.true = 100

# Simulate W (both training and test set) from the IBP prior
WWW.sim = prior.IBP(n, n.pred, eta.true,delta.true,gamma.true)

# W for the training set
W = WWW.sim[1:n, ]
# W in the test set
W.pred = WWW.sim[(n+1):(n+n.pred), ]

# Number of trees used in the simulated data
ntrees = ncol(W)

# Design matrix, Training set
mean_X=rnorm(p,0,10)
sd_X=rep(1,p)
x = matrix(rep(0,n*p),nrow = n)
for (j in 1:p){x[,j]=rnorm(n,mean_X[j],sd_X[j])}
## x = matrix(rnorm(n*p,0,1),nrow=n)
# Design matrix, Test set
X.pred = matrix(rnorm(n.pred*p,0,1),nrow=n.pred)
# Design matrix of both train and test sets 
X.all = rbind(x,X.pred)

# Simulate Trees from the BART prior
trees.sim = prior.trees(x,X.pred,WWW.sim, node_min_size = 5, alpha = 0.9, beta = 2, mu_mu = 0, sigma2 = 1, sigma2_mu = 1)

# Simulate Y from the model
if (dim(WWW.sim)[2]>1){
  y.sim = get_predictions_ibp( trees.sim, rbind(x,X.pred), WWW.sim, single_tree = FALSE)
} else{y.sim = get_predictions_ibp( trees.sim, rbind(x,X.pred), WWW.sim, single_tree = TRUE)  }

# Training set of Y
y_scale = y.sim[1:n]
y=y_scale
# Prediction set of Y
y_true_pred = y.sim[(n+1):(n+n.pred)]


##########
## Compare BART and BART.IBP 
##########

# Fit  BART
bart.fit = bart(x,y,X.pred,nburn = 1000,npost = 500,nthin = 1)
# Fit BART-IBP
bart.ibp.fit = bart.ibp(x,y,X.pred,nburn = 1000,npost = 500,nthin = 1)

# compare predictions BART and BART.IBP
pred_bart = colMeans(bart.fit$y_predicted)
pred_bart_ibp = colMeans(bart.ibp.fit$y_predicted)
mse_bart = mean((pred_bart-y_true_pred)^2)
mse_bart_ibp = mean((pred_bart_ibp-y_true_pred)^2)
mse_bart
mse_bart_ibp

save.image("//lancs/homes/38/battisto/My Desktop/BART_IBP_current/BART_IBP_current/results/prova.RData")

