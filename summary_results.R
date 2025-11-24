
source("aux_tree.R")


##################
### Variable Importance function
##################

# Description: Given the trees from the MCMC, compute the estimated variable importance
# the input 'bart' is the output of a bart or bart.ibp function
variable_importance = function(bart.ibp){
  
  # tree is a list of lists. The t-th element in the list is the list of trees of iteration t (after removing burnin and thinning) 
  trees = bart.ibp$trees
  # number of mcmc (after burnin and including only thinned iterations)
  store_size = length(trees)
  # Design matrix and number of explan. variables
  x = bart.ibp$x
  p = dim(x)[2]
  # matrix to store results
  variable_presence = matrix(rep(0,store_size*p), nrow = p)
  
  for (t in 1:store_size){
    # number of trees active at iteration t
    ntrees = length(trees[[t]])
    # we count how many times each of the p expl.var. is used among the ntrees active at iteration t
    counts = rep(0,p)
    for (j in 1:ntrees){
      variables_used = as.numeric(names(table(trees[[t]][[j]]$tree_matrix[,5])))
      counts[variables_used]=counts[variables_used] + 1
    }
    variable_presence[,t] =counts/ntrees
  }
  return(rowMeans(variable_presence)) 
}



########################
## Point estimate and Credible intervals f(x) (In-sample)
########################

# Estimate of f(x) for x in-sample
# bart is the output of bart.ibp
posterior_f_x = function(bart.ibp) {
  
  # get the trees, out of the MCMC. It is a list (for each iteration) of lists (for each tree)
  trees = bart.ibp$trees
  store_size = length(trees)
  # get the design matrix (needed only if in-sample)
  x = bart.ibp$x 
  # get the W matrix
  W = bart.ibp$W
  n = dim(x)[1]
  f_x = matrix(rep(0, store_size* n), ncol = store_size)
  for (t in 1:store_size){
    f_x[ ,t] = get_predictions_ibp(trees[[t]], x, W[[t]], single_tree = FALSE) 
  }
  point_estimate = rowMeans(f_x)
  credible = matrix(rep(0,2*n),nrow =2)
  for(j in 1:n){
    credible[,j] = quantile(f_x[j,], probs = c(0.025,0.975))
  }
  
  return(list(point_estimate,credible))
} 



##########################
##### Partial Correlation Function
##########################


### Partial function
# variable is the variable that we want to find the correlation
# lower, upper, length.out: lower and upper bounds of the grid, length of the grid.
partial = function(bart.ibp, variable, lower = NULL, upper = NULL, length.out = 10){
  
  # Create a grid of values
  if (is.null(lower)){lower = 0}
  if (is.null(upper)){upper = 1}
  grid = seq(from = lower, to = upper, length.out = length.out)
  
  # Extract design matrix x, trees and W from the bart.ibp object 
  x = bart.ibp$x
  n = dim(x)[1]
  p = dim(x)[2]
  trees = bart.ibp$trees
  store_size = length(trees)
  W = bart.ibp$W
  
  # space to store results. For each point in grid, we compute the point estimate and
  # 95% credible ball for f(x_s)
  f_x_s_point = rep(0, length.out)
  f_x_s_credible  = matrix(rep(0, 2*length.out), nrow = 2)
  
  # loop over the grid point
  for ( i in 1:length.out){
    
    # create a new design matrix in which variable 'variable' is replaced by a specific
    # value in grid, at this iteration value 'i'
    x_s = x
    x_s[ , variable] = rep(grid[i],n)
    
    # compute f(x_s)^{(t)} for each iteration t of the MCMC (after removing burnin and thinning)
    f_x_s_i = rep(0, store_size)
    # This step below computes f^{*}_{k}(x_s) (mid of page 276 in CGM(10) paper (the k index is out t, number of iterations))
    for (t in 1:store_size){
      f_x_s_i[t] = mean(get_predictions_ibp(trees[[t]], x_s, W[[t]], single_tree = FALSE))  
    }
    # This step compute the point estimates and credible region for f(x_s), for each x_s in the grid
    f_x_s_point[i] = mean(f_x_s_i)
    f_x_s_credible[,i] = quantile(f_x_s_i, probs = c(0.025,0.975))
  }
 
  return(list(f_x_s_point,f_x_s_credible))
}


############################
####### Casual Effects - CATE and CATT
############################

### Casual effects
# treatment variable is the index (among 1:p) of the binary treatment assignment variable
cate_catt = function(bart.ibp, treat_var){
  
  x = bart.ibp$x
  n = dim(x)[1]
  p = dim(x)[2]
  
  trees = bart.ibp$trees
  store_size = length(trees)

  # indexes for treatment equal to one
  indexes = which(x[,treat_var]==1)
  
  cate = rep(0,2)
  catt = rep(0,2)
  
  # loop over the grid point
  for ( i in 1:2){
    # set treatment either to 0 or 1
    x_s = x
    x_s[ , treat_var] = rep(i-1,n)
    
    # compute f(x_s)^{(t)} for each iteration t of the MCMC (after removing burnin and thinning)
    cate_i = rep(0, store_size)
    catt_i = rep(0, store_size)
    for (t in 1:store_size){
      cate_i[t] =mean( get_predictions(trees[[t]], x_s, single_tree = FALSE) ) 
      catt_i[t] =mean( get_predictions(trees[[t]], x_s, single_tree = FALSE)[indexes] )
    }
    cate[i]=mean(cate_i)
    catt[i]=mean(catt_i)
  }

  return(c(cate[2]-cate[1],catt[2]-catt[1]))
}




