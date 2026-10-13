

### Clean space
rm(list = ls())
### Set seed
set.seed(1000)

### Load Auxiliary functions and packages
source("simulate_data.R")
source("myBART.R")
source("myBART-IBP.R")
source("summary_results.R")

library(ggplot2)

###################
#### Simulate the data from Friedmann example
###################

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
bart.ibp.fit = bart.ibp(x,y,x.pred = X.pred, nburn = 500,npost = 100,nthin = 10)

# Variable importance
plot(variable_importance(bart.ibp.fit))

# In sample estimates of f(x) with credible regions
in_sample_poster_summary = posterior_f_x(prova.bart.ibp)

# Partial correlation of variables
# choose variable
j=2
partial_corr_j = partial(prova.bart.ibp,j,0,10,30)


# # (!!!)(TO DO) Plots 
# data <- data.frame(distance = c(10, 20, 30),
#                    value    = c(4, 5, 6),
#                    CI       = c(0.2, 0.5, 0.9))
# ggplot(data, aes(distance, value)) +
#   geom_point() +
#   geom_line() +
#   geom_errorbar(aes(ymin = value - CI, ymax = value + CI)) +
#   labs(x = "DISTANCE",
#        y = "VALUES",
#        title = "MY TITLE") +
#   theme_classic() 

# fit also BART
bart.fit = bart(x,y,x.pred = X.pred)

# compare predictions BART and BART.IBP
pred_bart = colMeans(bart.fit$y_predicted)
pred_bart_ibp = colMeans(bart.ibp.fit$y_predicted)
mse_bart = mean((pred_bart-y_true_pred)^2)
mse_bart_ibp = mean((pred_bart_ibp-y_true_pred)^2)
mse_bart
mse_bart_ibp




####################
#### Clustered Friedman example: Simulate from Friedman example in which groups use different explantory variables 
####################

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
bart.ibp.fit = bart.ibp(x,y,x.pred = X.pred, nburn = 2000,npost = 100,nthin = 10)
# Fit Bart
bart.fit = bart(x,y,X.pred, nburn = 2000,npost = 100,nthin = 10)

# compare predictions BART and BART.IBP
pred_bart = colMeans(bart.fit$y_predicted)
pred_bart_ibp = colMeans(bart.ibp.fit$y_predicted)
mse_bart = mean((pred_bart-y_true_pred)^2)
mse_bart_ibp = mean((pred_bart_ibp-y_true_pred)^2)
mse_bart
mse_bart_ibp

save.image("//lancs/homes/38/battisto/My Desktop/BART_IBP_current/BART_IBP_current/results/prova.RData")


