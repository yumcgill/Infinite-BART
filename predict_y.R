
prediction_of_y = function(y_hat_pred, W,  curr_trees,  X.pred, X.all,  eta,delta,
                           gamma, M.pred ,tree_fits_store, sigma2,K.new.max,mu_mu,sigma2_mu){
  
  n = dim(W)[1]

  # number of observations to predict
  n.pred = length(y_hat_pred)
  # current number of features
  K.curr = length(curr_trees)
  
  # compute trees counts and probabilities
  feature.counts = colSums(W)
  prob.pred = (feature.counts - eta)/(n + delta)

  # Initialize W.pred (W_{(n+1):(n+1+N)} from the prior, 
  #  given the current values of the hyperparameters and W_{1:n}
  W.pred  = matrix( rep(0,n.pred* dim(W)[2]),nrow=n.pred)
  for (ll in 1:n.pred){
    W.pred[ll,] = rbinom( dim(W)[2] ,1,prob.pred)
  }
  # W matrix 1:N (n rows in sample + n.pred rows out sample)
  W = rbind(W,W.pred)
  
  # Initialize y.pred
  if (dim(W)[2]>1){
  y_hat_pred = get_predictions_ibp(curr_trees, X.pred,W.pred, single_tree = FALSE)
  } else{  y_hat_pred = get_predictions_ibp(curr_trees, X.pred,W.pred, single_tree = TRUE) }
  
  tree_fits_store_pred = matrix(0, ncol = ncol(tree_fits_store), nrow = n.pred)
  for (j in 1:K.curr){
    tree_fits_store_pred[,j] = get_predictions_ibp(curr_trees[[j]], X.pred,W.pred[,j], single_tree = TRUE)
  }
  
  tree_fits_store = rbind( tree_fits_store , tree_fits_store_pred)
  
  # Update the parameter of the number of new features
  K.new.par = exp( log(gamma) + lgamma(1+delta) + lgamma(n-1+delta+eta ) - lgamma(n+delta) - lgamma(delta+eta) )
  
  # For the number of prediction steps, alternate between re-sampling W.pred and y.pred
  for (mm in 1:M.pred){
  
    # Resample W.pred  (given y_hat_pred, and the rest)
  
      for (w in (n+1):(n+n.pred)){
      
        ##### Update row w in W
        # compute the number of customers trying each dish
        if (ncol(W)>1){
          # compute m_k, the number of data points associated to tree k, for each k
          counts = colSums( W[-w, ])
          
          # remove empty columns (columns of all zeros) and corresponding trees
          empty.columns= which(counts == 0)
          if (length(empty.columns)>0 ){
            W = W[,-empty.columns,drop = FALSE]  
            counts = counts[-empty.columns]
            tree_fits_store = tree_fits_store[,-empty.columns,drop = FALSE]
            curr_trees = curr_trees[-empty.columns]
          }
        } else{counts = sum(W)}
      
        #### Observation w selects trees among the currently existing ones
        if (ncol(W)>1){
          for (k in 1:dim(W)[2]){
            # compute the fitted value, y_hat computed without using k-th tree
            fit.no.k = tree_fits_store[w,-k] %*% W[w,-k] 
            # compute the probability of observation w not using tree k
            prob0 = (n + n.pred - 1 + delta - counts[k] + eta ) * dnorm( y_hat_pred[w-n], fit.no.k, sqrt(sigma2)) + 10^(-100)
            # compute additional fit
            add.fit = get_predictions_ibp(curr_trees[[k]], X.pred,1, single_tree = TRUE)[w-n]
            # compute the probability of observation w joining tree k
            prob1 = (counts[k] - eta) * dnorm( y_hat_pred[w-n], fit.no.k + add.fit , sqrt(sigma2)) + 10^(-100) 
            # Sample W_{w,k} 
            W[w,k]=sample(c(0,1),1, prob = c(prob0,prob1) )
            if (W[w,k]==1){
              tree_fits_store[w,k] = add.fit
            } else{ tree_fits_store[w,k] = 0}
          }
        } else{
          # compute the probability of observation w not using tree k
          prob0 = (n + n.pred - 1 + delta - counts + eta ) * dnorm( y_hat_pred[w-n], 0 , sqrt(sigma2)) + 10^(-100)
          # Get the predicted (mean) value of y_w if we use the k-th tree
          add.fit = get_predictions_ibp(curr_trees[[1]], X.pred,1, single_tree = TRUE)[w-n]
          # compute the probability of observation w joining tree k
          prob1 = (counts - eta) * dnorm( y_hat_pred[w-n],  add.fit, sqrt(sigma2)) + 10^(-100)
          # Sample W_{w,k} 
          W[w,1]=sample(c(0,1),1, prob = c(prob0,prob1) )
          if (W[w,1]==1){
            tree_fits_store[w,1] = add.fit
          } else{ tree_fits_store[w,1] = 0}
        }
      
        #### Observation w selects new trees
        prob.new= rep(0, K.new.max + 1)
        ### CAN WE INTEGRATE THIS OUT ????
        new.mu = rnorm(K.new.max , mu_mu ,  sigma2_mu )
        new.mu.cum = cumsum(new.mu)
        fit.no.k = tree_fits_store[w,] %*% W[w,]
        # Probability of not generating any new tree
        prob.new[1] =  dpois( 0, K.new.par ) * dnorm( y_hat_pred[w-n], fit.no.k , sqrt(sigma2))
        # Probability of generating k new trees
        for (k in 1:K.new.max){
          prob.new[k+1] = dpois(k, K.new.par )* dnorm( y_hat_pred[w-n], fit.no.k + new.mu.cum[k], sqrt(sigma2))
        }
        # select how many new trees (the number of new trees is K.new-1, i.e. if K.new=1, no new trees are added)
        if (max(prob.new)>10^(-100)){
        prob.new=prob.new/sum(prob.new)
        K.new = sample(1:(K.new.max+1),1, prob = prob.new)
        } else{K.new = 0}
      
        # Add new columns to W and new trees
        if (K.new >1){
          W.new = matrix(rep(0,(K.new-1)*(n+n.pred)) , nrow=(n+n.pred))
          W.new[w,]= rep(1,K.new-1)
          W = cbind(W,W.new)
          W.pred = W[(n+1):(n+n.pred),]
        
          # Add new trees  
          trees.new = create_stump_ibp(num_trees = K.new - 1 , y = y_hat_pred[w-n] , X = X.pred[ w-n ,], W.new)
          tree_fits_store.new = matrix(rep(0,(K.new-1)*(n+n.pred)), nrow=(n+n.pred))
          for (k in 1:(K.new-1)){
            trees.new[[k]][[1]][1,7]=new.mu[k]
            if (is.null(dim(W.new))){
              tree_fits_store.new[,k] = get_predictions_ibp(trees.new[[k]], X.all , W.new, single_tree = TRUE)
            } else {
              tree_fits_store.new[,k] = get_predictions_ibp(trees.new[[k]], X.all , W.new[,k], single_tree = TRUE)
            }
          }
          curr_trees = c(curr_trees,trees.new)
          tree_fits_store = cbind(tree_fits_store,tree_fits_store.new)
        
        }
      } # end loop n for W
    
    
      ### update allocations of observations to terminal nodes and fitted values
    
      curr_trees_p = curr_trees
    
      for (p in 1:dim(W)[2]){
        if (sum(W[,p])>1 & nrow(curr_trees[[p]]$tree_matrix)>1 ){
          curr_trees[[p]] = fill_tree_details(curr_trees[[p]], X.all[which(W[,p]==1),])
        } 
        if (sum(W[,p])>1 & nrow(curr_trees[[p]]$tree_matrix)==1 ){
          curr_trees[[p]]$node_indices = rep(1,sum(W[,p]))
          curr_trees[[p]]$tree_matrix[1,'node_size'] = sum(W[,p])
        } 
        if (sum(W[,p])==1){
          new_trees= create_stump_ibp(1,c(y,y_hat_pred)[which(W[,p]==1)],X.all[which(W[,p]==1),],1)
          curr_trees[[p]]= new_trees[[1]]
        }
      
        ######## IF A LEAF NODE BECOMES EMPTY I PRUNE THE TREE
        while( length(which(as.numeric(curr_trees[[p]]$tree_matrix[,'node_size']) ==0 )) > 0  )   {
          curr_trees[[p]]=clean_prune_tree(X.all[which(W[,p]==1),],c(y,y_hat_pred)[which(W[,p]==1)],curr_trees[[p]])
        }
      
        tree_fits_store[,p] = get_predictions_ibp(curr_trees[[p]], X.all,W[,p], single_tree = TRUE)
      }
      y_hat = rowSums(tree_fits_store)
    
    
      # update number of trees
      ntrees = dim(W)[2]
  
    
    
    
    
    # Re-sample y.pred  (given W.pred, and the rest)
      if (dim(W)[2]>1) {
           y_hat_pred = get_predictions_ibp(curr_trees, X.pred,W.pred, single_tree = FALSE)
      } else{ y_hat_pred = get_predictions_ibp(curr_trees, X.pred,W.pred, single_tree = TRUE)}
    }  

  return(y_hat_pred)
}


