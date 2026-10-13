######################
# 1974 Motor Trend US magazine and comprises 32 observations on 11 variables
######################
source("simulate_data.R")
source("myBART.R")
source("myBART-IBP.R")
# Load data
data(mtcars)
#View(car_data)
names(mtcars)
set.seed(21)
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
bart.ibp.fit_mtcar = bart.ibp(x,y,x.pred = X.pred,nburn = 1,npost = 6000,nthin = 1)
saveRDS(bart.ibp.fit_mtcar, "bart.ibp.fit_mtcar1.rds")