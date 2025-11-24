

### Load Auxiliary functions and packages
source("aux_tree.R")
source("modify_tree.R")
source("predict_tree.R")
source("parameters_quantities_tree.R")
source("IBP_functions.R")
source("simulate_data.R")
source("predict_y.R")

bart.ibp = function(x,y,x.pred = NULL,sparse = TRUE,node_min_size = 5,alpha = 0.95,beta = 2,
          nu = 3,lambda = NULL,mu_mu = 0, sigma2 = 1,sigma2_mu = 1,lambda_cov = 0.4,nu_cov = 2,
          nburn = 20000,npost = 1000,nthin = 10,penalise_num_cov = TRUE,
          eta.alpha = 0.05,eta.beta = 0.01,delta.alpha = 0.1,delta.beta = 0.01,gamma.alpha = 0.05,gamma.beta = 0.01, # prior hyperparmaters
          sigma.eta= 1,sigma.delta= 1,K.new.max = 5,M.pred = 10,eta = -100,delta = 101,gamma = 100 # MH variances, MCMC par and initial values
        ) {

  
  # Scale the response target variable
  y_mean = mean(y)
  y_sd = sd(y)
  y_scale = (y - y_mean)/y_sd
  
  # Don't scale the response variable
  # y_scale=y
  # y_mean=0
  # y_sd=1
  
  
  n = length(y_scale)
  p = ncol(x)
  s = rep(1/p, p)
  
## TO DO: DOUBLE CHECK HERE y_sd should be the sd of the already standardized variables, so 1???  
  if (is.null(lambda)){
    q = 0.9
    lambda = qgamma(1-q,nu, rate =1/2)*(y_sd^2)/nu
  }
  
  # Initialize W (training set n, only) from the IBP prior
  W = prior.IBP(n,0,eta,delta,gamma)
  
  # Initial value of number of trees
  ntrees = ncol(W)
   
  #####################
  ##### Set storage space for MCMC and initialize the remaining parameters
  ####################
  
  # Initialize Trees: Create a list of trees for the initial stump
  curr_trees = create_stump_ibp(num_trees = ntrees, y = y_scale, X = x, W)
  new_trees = curr_trees
  
  # Initialize Fitted values (basically to zero)
  y_hat = get_predictions_ibp(curr_trees, x, W, single_tree = ntrees == 1)

  # Total number of iterations  
  TotIter = nburn + npost*nthin 
  
  # Set Storage space 
  # (Remark: if interest is only in a specific summary, better to store less and save memory;
  # currently I am saving everything so that I don't have to re-run again.)
  store_size = npost
  tree_store = vector('list', store_size)
  sigma2_store = rep(NA, store_size)
  y_hat_store = matrix(NA, ncol = length(y_scale), nrow = store_size)
  var_count = rep(0, ncol(x))
  var_count_store = matrix(0, ncol = ncol(x), nrow = store_size)
  s_prob_store = matrix(0, ncol = ncol(x), nrow = store_size)
  tree_fits_store = matrix(0, ncol = ntrees, nrow = length(y))
  IBP.hyper.store= matrix(rep(0,store_size*3),ncol=3)
  IBP.accept.store= matrix(rep(0,TotIter*2),ncol=2)
  ntrees.store = rep(NA, store_size)
  W_store = vector('list',store_size)
  # initialize counters to check the acceptance rates of eta and delta
  accept.eta=0
  accept.delta=0
  # Additional storage if doing out of sample prediction
  if (is.null(x.pred) == FALSE){
    X.all = rbind(x,x.pred)
    n.pred = dim(x.pred)[1]
    y_hat_pred = rep(0,n.pred)
    y_pred_store = matrix(0, ncol = n.pred, nrow = store_size)
    # storage space for predicted values
    y_hat_pred_store = matrix(NA, ncol = n.pred, nrow = store_size)
    
  } else{ 
    y_pred_store = 0}
  
  # Set up a progress bar and counter for running time
  pb = utils::txtProgressBar(min = 1, max = TotIter,style = 3, width = 60, title = 'Running IBP-BART...')
  
  
  ####################  
  ####### MCMC Sampler 
  ####################
  
  for (i in 1:TotIter) {
   
    utils::setTxtProgressBar(pb, i)
    
    ###################
    #####  Store values and do predictions for y.pred
    ###################
    
    # Start storing results only after burnin and for non-thinned iterations
    if((i > nburn) & ((i - nburn) %% nthin) == 0) {
      
      ###  Store everything 
      curr = (i - nburn)/nthin
      tree_store[[curr]] = curr_trees
      W_store[[curr]] = W
      sigma2_store[curr] = sigma2
      y_hat_store[curr,] = y_hat
      var_count_store[curr,] = var_count
      s_prob_store[curr,] = s
      IBP.hyper.store[curr,]=c(eta,delta,gamma)
      ntrees.store[curr]=ntrees
      
      ### Prediction Step
      if (is.null(x.pred) == FALSE){
      # do prediction of out-of-sample y
      y_hat_pred = prediction_of_y( y_hat_pred,  W,  curr_trees, x.pred,X.all,  
                        eta,delta,gamma, M.pred, tree_fits_store,sigma2,K.new.max,mu_mu,sigma2_mu)
      # Store the predicted values y_hat
      y_hat_pred_store[curr,] = y_hat_pred*y_sd + y_mean
      }
    }
    
    
    ################
    ### UPDATE THE TREES
    ################
    
    for (j in 1:ntrees) {
      
      if (length(curr_trees[[j]]$node_indices)>node_min_size)     {    
        
        current_partial_residuals = y_scale - y_hat + tree_fits_store[,j]
        
        # Choose the type of new tree proposed in MH among: grow/change/prune/swap
        type = sample_move(curr_trees[[j]], i, nburn)
        
        # Generate a NEW TREE based on the current tree and chosen proposal move
        # (Remark: to update tree j, we use only obs assigned to that tree, by using the j-th column of W)
        new_trees[[j]] = update_tree(y = y_scale[which(W[,j]==1)], X = x[which(W[,j]==1),], type = type, curr_tree = curr_trees[[j]], node_min_size = node_min_size, s = s)
        
        # CURRENT TREE: compute the log of the marginalized likelihood + log of the tree prior
        l_old = tree_full_conditional(curr_trees[[j]], current_partial_residuals[W[,j]==1], sigma2, sigma2_mu) + get_tree_prior(curr_trees[[j]], alpha, beta)
        
        # NEW TREE: compute the log of the marginalized likelihood + log of the tree prior
        l_new = tree_full_conditional(new_trees[[j]], current_partial_residuals[W[,j]==1], sigma2, sigma2_mu) +  get_tree_prior(new_trees[[j]], alpha, beta)
        
        # Compute the MH acceptance probability
        a = exp(l_new - l_old)
        
        # Accept MH move 
        if(a > runif(1)) {
          curr_trees[[j]] = new_trees[[j]]
        }
        
        # Update mu whether tree accepted or not
        curr_trees[[j]] = simulate_mu(curr_trees[[j]], current_partial_residuals[W[,j]==1], sigma2, sigma2_mu)
        
        # Updating BART predictions
        current_fit = get_predictions_ibp(curr_trees[[j]], x,W[,j], single_tree = TRUE)  # predict using only tree j (gives 0 for obs not having feature j)
        y_hat = y_hat - tree_fits_store[,j] # subtract the old fit
        y_hat = y_hat + current_fit # add the new fit
        tree_fits_store[,j] = current_fit # update the new fit
        
      } # end if node_min_size. If the j-th tree doesn't have enough observations, we don't update it
      
    } # End loop through trees
    
    
    ##################
    ###  Update the Matrix W
    ##################
    
    # Update the parameter of the number of new features K
    K.new.par = exp( log(gamma) + lgamma(1+delta) + lgamma(n-1+delta+eta ) - lgamma(n+delta) - lgamma(delta+eta) )
    
    # Update each row of W
    for (w in 1:n){
      
      ### Update row w in W
      
      # (Remark: this ugly condition ncol(W)>1 is because R turns 1-column matrices into vectors,
      #  I should have fixed that (by keeping W always as a matrix) but need to double check it.)
      if (ncol(W)>1){
        
        # compute m_k, the number of data points associated to tree k, for each k
        counts = colSums( W[-w, ])
        
        # remove empty columns (columns of all zeros) and corresponding trees (these are trees used only by obs w)
        empty.columns= which(counts == 0)
        
        # this is just to avoid that all columns could possibly be dropped (should almost never occur)
        if (length(empty.columns) != dim(W)[2]){
           if (length(empty.columns)>0 ){
            W = W[,-empty.columns,drop = FALSE]  
            counts = counts[-empty.columns]
            tree_fits_store = tree_fits_store[,-empty.columns,drop = FALSE]
            curr_trees = curr_trees[-empty.columns]
          }
          
          ## Observation w selects trees among the currently existing ones. 
          #  We sample K Bernoulli r.v.s, with the probabilities of 0 and 1 proportional to prob0 and prob1
          for (k in 1:dim(W)[2]){
            # compute the fitted value, y_hat computed without using k-th tree
            fit.no.k = tree_fits_store[w,-k] %*% W[w,-k] 
            # compute the probability of observation w not using tree k
            prob0 = (n - 1 + delta - counts[k] + eta ) * dnorm( y_scale[w], fit.no.k, sqrt(sigma2)) + 10^(-100)
            # Compute the possible additional fit added by tree k, i.e. predicted (mean) value of y_w if we use the k-th tree
            add.fit = get_predictions_ibp(curr_trees[[k]], x,1, single_tree = TRUE)[w]
            # compute the probability of observation w joining tree k
            prob1 = (counts[k] - eta) * dnorm( y_scale[w], fit.no.k +  add.fit , sqrt(sigma2)) + 10^(-100)
            # Sample W_{w,k} 
            W[w,k]=sample(c(0,1),1, prob = c(prob0,prob1) )
            
            # Update the fit store matrix 
            if (W[w,k]==1){
              tree_fits_store[w,k] = add.fit
            } else{ tree_fits_store[w,k] = 0}
          } }
      } else{  ## if W has only one column
        
       counts = sum(W[-w])
       if (counts>0){
          # compute the probability of observation w not using tree k
          prob0 = (n - 1 + delta - counts + eta ) * dnorm( y_scale[w], 0 , sqrt(sigma2)) + 10^(-100)
          # Get the predicted (mean) value of y_w if we use the k-th tree
          add.fit = get_predictions_ibp(curr_trees[[1]], x,1, single_tree = TRUE)[w]
          # compute the probability of observation w joining tree k
          prob1 = (counts - eta) * dnorm( y_scale[w],  add.fit,  sqrt(sigma2)) + 10^(-100)
          # Sample W_{w,k} 
          W[w,1]=sample(c(0,1),1, prob = c(prob0,prob1) )
          # Update the fit store matrix
          if (W[w,1]==1){
            tree_fits_store[w,1] = add.fit
          } else{ tree_fits_store[w,1] = 0}  
        }
        
      }
      
      
      #### Observation w selects new trees
      
      prob.new= rep(0, K.new.max + 1)
      
      # (Comment: Can we integrate this step out?)
      new.mu = rnorm(K.new.max , mu_mu ,  sigma2_mu )
      new.mu.cum = cumsum(new.mu)
      fit.no.k = tree_fits_store[w,] %*% W[w,]
      # Probability of not generating any new tree
      prob.new[1] =  dpois( 0, K.new.par ) * dnorm( y_scale[w], fit.no.k , sqrt(sigma2))
      # Probability of generating k additional new trees
      for (k in 1:K.new.max){
        prob.new[k+1] = dpois(k, K.new.par )* dnorm( y_scale[w], fit.no.k + new.mu.cum[k], sqrt(sigma2))
      }
      # select how many new trees (the number of new trees is K.new-1, i.e. if K.new=1, no new trees are added)
      if (max(prob.new)>10^(-100)){
      prob.new=prob.new/sum(prob.new)
      K.new = sample(1:(K.new.max+1),1, prob = prob.new)
      } else{K.new = 0}
      
      # Add new columns to W and new trees
      if (K.new >1){
        W.new = matrix(rep(0,(K.new-1)*n), nrow=n)
        W.new[w,]= rep(1,K.new-1)
        W = cbind(W,W.new)
        
        # Create and Add the new trees to curr_trees and tree_fits_store  
        #  (Comment: Double check this part)
        trees.new = create_stump_ibp(num_trees = K.new - 1 , y = y_scale[w], X = x[w,], W.new)
        tree_fits_store.new = matrix(rep(0,(K.new-1)*n), nrow=n)
        for (k in 1:(K.new-1)){
          trees.new[[k]][[1]][1,7]=new.mu[k]
          if (is.null(dim(W.new))){
            tree_fits_store.new[,k] = get_predictions_ibp(trees.new[[k]], x, W.new, single_tree = TRUE)
          } else {
            tree_fits_store.new[,k] = get_predictions_ibp(trees.new[[k]], x, W.new[,k], single_tree = TRUE)
          }
        }
        # Attached the newly created trees to the list with the others
        curr_trees = c(curr_trees,trees.new)
        tree_fits_store = cbind(tree_fits_store,tree_fits_store.new)
      }
    } # end loop in n for W
    
    
    
    ### Update allocations of observations to terminal nodes and fitted values
    
    curr_trees_p = curr_trees
    
    for (p in 1:dim(W)[2]){
      if (sum(W[,p])>1 & nrow(curr_trees[[p]]$tree_matrix)>1 ){
        curr_trees[[p]] = fill_tree_details(curr_trees[[p]], x[which(W[,p]==1),])
      } 
      if (sum(W[,p])>1 & nrow(curr_trees[[p]]$tree_matrix)==1 ){
        curr_trees[[p]]$node_indices = rep(1,sum(W[,p]))
        curr_trees[[p]]$tree_matrix[1,'node_size'] = sum(W[,p])
      } 
      if (sum(W[,p])==1){
        new_trees= create_stump_ibp(1,y[which(W[,p]==1)],x[which(W[,p]==1),],1)
        curr_trees[[p]]= new_trees[[1]]
      }
      
      # If a leaf node becomes empty, we need to prune the tree
      while( length(which(as.numeric(curr_trees[[p]]$tree_matrix[,'node_size']) ==0 )) > 0  )   {
        curr_trees[[p]]=clean_prune_tree(x[which(W[,p]==1),],y[which(W[,p]==1)],curr_trees[[p]])
      }
      
      # (Comment: this should be reduntant. We already update this above.)
      tree_fits_store[,p] = get_predictions_ibp(curr_trees[[p]], x,W[,p], single_tree = TRUE)
      
    }
    
    # update number of trees
    ntrees = dim(W)[2]
    
    
    ################
    ### Update other BART parameters
    ################
    
    # Update y_hat
    y_hat = rowSums(tree_fits_store)
    
    # y_hat = get_predictions_ibp(curr_trees, x, single_tree = ntrees == 1)
    sum_of_squares = sum((y_scale - y_hat)^2)
    
    # Update sigma2 (variance of the residuals)
    sigma2 = update_sigma2(sum_of_squares, n = length(y_scale), nu, lambda)
    
    #  TO DO !!!
    # var_count =
    
    # Update s = (s_1, ..., s_p), where s_p is the probability that predictor p is used to create new terminal nodes
    if (sparse == 'TRUE' & i > floor(TotIter*0.1)){
      s = update_s(var_count, p, 1)
    }
    
    ##############
    ### Update IBP hyperparameters (eta,delta,gamma)
    ##############
    
    #up.eta = update_eta(W,accept.eta,eta,delta,gamma, sigma.eta,sigma.delta, eta.alpha,eta.beta,delta.alpha,delta.beta)
    up.etadelta = update_eta_delta_block(W,accept.eta,eta,delta,gamma, sigma.eta,sigma.delta, eta.alpha,eta.beta,delta.alpha,delta.beta)
    eta = up.etadelta[1]
    accept.eta = up.etadelta[3]
    
    #up.delta = update_delta(W,accept.delta,eta,delta,gamma, sigma.eta,sigma.delta, delta.alpha ,delta.beta)
    delta = up.etadelta[2]
    accept.delta = up.etadelta[3]
    
    gamma = update_gamma(W,eta,delta,gamma, gamma.alpha,gamma.beta)
    
    IBP.accept.store[i,]= c(accept.eta,accept.delta)
    
    ### This is some adaptation of the MH variances
    if ( i%%100 == 0 & i < nburn){
      accept_rate.eta = (IBP.accept.store[i,1] - IBP.accept.store[i-99,1])/100
      accept_rate.delta  = (IBP.accept.store[i,2] - IBP.accept.store[i-99,2])/100
      if (accept_rate.eta >0.5){sigma.eta=sigma.eta*1.1}
      if (accept_rate.eta <0.1){sigma.eta=sigma.eta*0.9}
      if (accept_rate.delta >0.5){sigma.delta=sigma.delta*1.1}
      if (accept_rate.delta <0.1){sigma.delta=sigma.delta*0.9}
        }
    
    print(c(i,ntrees) )
    
  } # End MCMC iterations loop
  
  cat('\n') # Make sure progress bar ends on a new line
  
  
  if (is.null(x.pred) == TRUE){
    y_pred_store = NULL
    y_hat_pred_store = NULL
    
  }
     return(list(trees = tree_store,   # 1
                sigma2 = sigma2_store*y_sd^2, # 2
                y_hat = y_hat_store*y_sd + y_mean, # 3
                npost = npost,  # 4
                nburn = nburn,  # 5
                nthin = nthin,  # 6
                ntrees = ntrees.store,  # 7
                x = x,  # 8 
                y_mean = y_mean,  # 9 
                y_sd = y_sd,  # 10
                var_count_store = var_count_store,  # 11
                s = s_prob_store,  # 12
                IBP.hyper = IBP.hyper.store,  # 13
                y_fitted = y_pred_store,  # 14
                y_predicted = y_hat_pred_store,  # 15
                IBP.accept.store = IBP.accept.store,  # 16
                W = W_store  # 17
                ))
  
} # End main function
  
  
