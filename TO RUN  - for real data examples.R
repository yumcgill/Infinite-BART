
### Clean space
rm(list = ls())

### Set seed
set.seed(1000)

### Load Auxiliary functions and packages
source("simulate_data.R")
source("myBART.R")
source("myBART-IBP.R")

#### SET YOUR WORKING DIRECTORY HERE!!!  
# setwd("C:/Users/marco/Desktop/BART_IBP_current/BART_IBP/")
#setwd("C:/Users/battisto/OneDrive - Lancaster University/BART_IBP")



##############
#### CAR example
##############

# Load Data
car_data <- read.csv("datasets/car/car_data.csv")
#View(car_data)
names(car_data)

# Set design matrix and response vector (full: both train and test) 
x.full = model.matrix( ~ Car_Name + Year + Selling_Price + Kms_Driven + Fuel_Type + Seller_Type + Transmission + Owner, data=car_data )
y.full = log(car_data$Present_Price)
n.full = length(y.full)
n = round(n.full*0.8)

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
bart.ibp.fit = bart.ibp(x,y,x.pred = X.pred)
bart.fit = bart(x,y,x.pred = X.pred)

# Compare the predictived values of BART and BART.IBP
pred_bart = colMeans(bart.fit$y_predicted)
pred_bart_ibp = colMeans(bart.ibp.fit$y_predicted)
mse_bart = mean((pred_bart-y_true_pred)^2)
mse_bart_ibp = mean((pred_bart_ibp-y_true_pred)^2)
mse_bart
mse_bart_ibp

# Save the results
save.image("results/car_results_1000.RData")

# Check traceplots of number of trees and hyperparameters of IBP
plot(bart.ibp.fit$ntrees)
plot(bart.ibp.fit$IBP.hyper[,1])
plot(bart.ibp.fit$IBP.hyper[,2])
plot(bart.ibp.fit$IBP.hyper[,3])



###############
###### Motorbikes example
###############

# Load data
motor_data <- read.csv("datasets/bike/bike.csv")
#View(motor_data)
names(motor_data)

# remove missing values
motor_data = motor_data[complete.cases(motor_data), ]

# Set design matrix and response vector (full: both train and test)
x.full = model.matrix( ~ name + year + km_driven + seller_type + ex_showroom_price + owner, data=motor_data )
y.full = log(motor_data$selling_price)
n.full = length(y.full)
n = round(n.full*0.8)

# split in training and test set
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
bart.ibp.fit = bart.ibp(x,y,x.pred = X.pred)
bart.fit = bart(x,y,x.pred = X.pred)

# Compare predictions
pred_bart = colMeans(bart.fit$y_predicted)
pred_bart_ibp = colMeans(bart.ibp.fit$y_predicted)
mse_bart = mean((pred_bart-y_true_pred)^2)
mse_bart_ibp = mean((pred_bart_ibp-y_true_pred)^2)
mse_bart
mse_bart_ibp

# Save results
save.image("results/motor_results_1000.RData")

# Check traceplots of number of trees and hyperparameters of IBP
plot(bart.ibp.fit$ntrees)
plot(bart.ibp.fit$IBP.hyper[,1])
plot(bart.ibp.fit$IBP.hyper[,2])
plot(bart.ibp.fit$IBP.hyper[,3])



######################
# 1974 Motor Trend US magazine and comprises 32 observations on 11 variables
######################

# Load data
data(mtcars)
#View(car_data)
names(mtcars)

# Set design matrix and response vector (full: both train and test)
x.full = model.matrix( ~ cyl + disp + hp + drat + wt + qsec + vs + am + gear + carb , data=mtcars )
y.full = log(mtcars$mpg)
n.full = length(y.full)
n = round(n.full*0.8)

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

# Fit BART and BART.IBP
bart.ibp.fit = bart.ibp(x,y,x.pred = X.pred)
bart.fit = bart(x,y,x.pred = X.pred)

# Compare predictions
pred_bart = colMeans(bart.fit$y_predicted)
pred_bart_ibp = colMeans(bart.ibp.fit$y_predicted)
mse_bart = mean((pred_bart-y_true_pred)^2)
mse_bart_ibp = mean((pred_bart_ibp-y_true_pred)^2)
mse_bart
mse_bart_ibp

# Save the results
save.image("results/mtcars_1000.RData")

# Check traceplots of number of trees and hyperparameters of IBP
plot(bart.ibp.fit$ntrees)
plot(bart.ibp.fit$IBP.hyper[,1])
plot(bart.ibp.fit$IBP.hyper[,2])
plot(bart.ibp.fit$IBP.hyper[,3])



##################
### Boston housing 506 observations on 14 variables
##################

# load data
library(MASS)
data(Boston)
names(Boston)

# Set design matrix and response vector (full: both train and test)
x.full = model.matrix( ~ crim + zn + indus + chas + nox + rm + age + dis + rad + tax + ptratio + black +lstat , data=Boston )
y.full = log(Boston$medv)
n.full = length(y.full)
n = round(n.full*0.8)

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

# Fit BART and BART.IBP
bart.ibp.fit = bart.ibp(x,y,x.pred = X.pred)
bart.fit = bart(x,y,x.pred = X.pred)

# compare predictions
pred_bart = colMeans(bart.fit$y_predicted)
pred_bart_ibp = colMeans(bart.ibp.fit$y_predicted)
mse_bart = mean((pred_bart-y_true_pred)^2)
mse_bart_ibp = mean((pred_bart_ibp-y_true_pred)^2)
mse_bart
mse_bart_ibp

# Save results
save.image("results/Boston_1000.RData")

# Check traceplots of number of trees and hyperparameters of IBP
plot(bart.ibp.fit$ntrees)
plot(bart.ibp.fit$IBP.hyper[,1])
plot(bart.ibp.fit$IBP.hyper[,2])
plot(bart.ibp.fit$IBP.hyper[,3])



################
#### Airquality 153 observations on 6 variables
###############

# Load data
data(airquality)
airquality = airquality[complete.cases(airquality), ]
names(airquality)

# Set design matrix and response vector (full: both train and test)
x.full = model.matrix( ~ Solar.R + Wind + Temp + Month + Day, data=airquality )
y.full = log(airquality$Ozone)
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

# Fit BART and BART.IBP
bart.ibp.fit = bart.ibp(x,y,x.pred = X.pred)
bart.fit = bart(x,y,x.pred = X.pred)

# compare prediction
pred_bart = colMeans(bart.fit$y_predicted)
pred_bart_ibp = colMeans(bart.ibp.fit$y_predicted)
mse_bart = mean((pred_bart-y_true_pred)^2)
mse_bart_ibp = mean((pred_bart_ibp-y_true_pred)^2)
mse_bart
mse_bart_ibp

# Save results
save.image("results/airquality_1000.RData")

# Check traceplots of number of trees and hyperparameters of IBP
plot(bart.ibp.fit$ntrees)
plot(bart.ibp.fit$IBP.hyper[,1])
plot(bart.ibp.fit$IBP.hyper[,2])
plot(bart.ibp.fit$IBP.hyper[,3])



################
#### Aquatic Toxicity - https://archive.ics.uci.edu/dataset/505/qsar+aquatic+toxicity 
###############

# Load data
aquatic_toxicity <- read.csv("datasets/toxicity/qsar_aquatic_toxicity.csv", header=FALSE, sep=";")
# View(qsar_aquatic_toxicity)

# Set design matrix and response vector (full: both train and test) 
x.full = model.matrix( ~ V1 + V2 + V3 + V4 + V5 + V6 + V7 + V8, data=aquatic_toxicity )
y.full = log(aquatic_toxicity$V9)
n.full = length(y.full)
n = round(n.full*0.8)

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
bart.ibp.fit = bart.ibp(x,y,x.pred = X.pred)
bart.fit = bart(x,y,x.pred = X.pred)

# Compare the predictived values of BART and BART.IBP
pred_bart = colMeans(bart.fit$y_predicted)
pred_bart_ibp = colMeans(bart.ibp.fit$y_predicted)
mse_bart = mean((pred_bart-y_true_pred)^2)
mse_bart_ibp = mean((pred_bart_ibp-y_true_pred)^2)
mse_bart
mse_bart_ibp

# Save results
save.image("results/toxicity_1000.RData")

# Check traceplots of number of trees and hyperparameters of IBP
plot(bart.ibp.fit$ntrees)
plot(bart.ibp.fit$IBP.hyper[,1])
plot(bart.ibp.fit$IBP.hyper[,2])
plot(bart.ibp.fit$IBP.hyper[,3])

#### (!!!) For clustered example, let us try to remove V7 (possibly also V3 and/or V8)

