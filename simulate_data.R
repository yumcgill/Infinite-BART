

##################
### Draw a value from the IBP prior
##################

prior.IBP <- function(n,n.pred, eta,delta,gamma){
  
  NN = n + n.pred
  K = rpois(1,gamma)
  if (K==0){K=1}
  
  W =matrix( rep(0, NN*K),nrow=NN)
  
  W[1,]= rep(1,K)
  
  W[2,] = rbinom(K,1,(1-eta)/(delta+1))
  K.new=rpois(1, exp( log(gamma) + lgamma(1+delta) + lgamma(2-1+delta+eta) - lgamma(2+delta) - lgamma(delta+eta) )  )
  W.new= matrix(rep(0,K.new*NN), nrow=NN)
  W.new[2,]=rep(1,K.new)
  W = cbind(W,W.new)
  
  for (i in 3:NN){
    K = ncol(W)
    if (K>1 ){
      counts = colSums(W[1:(i-1),])
    } else{ counts = sum(W)}
    
    prob = (counts - eta) / (i - 1 + delta )
    W[i, ]= rbinom(K,1,prob)
    
    K.new=rpois(1, exp( log(gamma) + lgamma(1+delta) + lgamma(i-1+delta+eta) - lgamma(i+delta) - lgamma(delta+eta) )  )
    W.new= matrix(rep(0,K.new*NN), nrow=NN)
    W.new[i,]=rep(1,K.new)
    W = cbind(W,W.new)
  }
  return(W)
}




###############
##### Function to generate trees from the BART prior
###############

prior.trees <- function( x , X.pred , WWW.sim , node_min_size = 5, 
                         alpha = 0.9, beta = 2, mu_mu = 0, sigma2 = 1,sigma2_mu = 1){
   
   X=rbind(x,X.pred)
    
   num_trees= ncol(WWW.sim)
   
   # Create holder for trees
   all_trees = vector('list', length = num_trees)
   
   # Loop through trees
   for (j in 1:num_trees) {
     
     ### which observations are allocated to tree j?
     obs_in_tree_j = which( WWW.sim[,j] == 1)
     
     ### number of observations in tree j
     nn_tree_j = length(obs_in_tree_j)
     
     ### Derive a design matrix with only observations allocated to tree j
     X_tree_j = X[ obs_in_tree_j, ]
     
     #  Set up each tree to have two elements in the list as described above
     all_trees[[j]] = vector('list', length = 2)
     # Give the elements names
     names(all_trees[[j]]) = c('tree_matrix',
                               'node_indices')
     # Create the two elements: first is a matrix
     all_trees[[j]][[1]] = matrix(NA, ncol = 8, nrow = 1)
     
     # Second is the assignment to node indices
     all_trees[[j]][[2]]= rep(1, nn_tree_j )
     
     # Create column names
     colnames(all_trees[[j]][[1]]) = c('terminal','child_left','child_right','parent','split_variable','split_value','mu','node_size')
     
     # Set values for stump
    all_trees[[j]][[1]][1,] = c(1, NA, NA, NA, NA, NA, NA , nn_tree_j )
  
    terminals = 1
    depth = 0
    
    while ( length(terminals)>0) {
      
      ## if a node among the possible splitting terminal node doesn't have enough observations to generate two children nodes
      ## with at least node_min_size nodes, we don't split it
      if ( all_trees[[j]]$tree_matrix[terminals[1],'node_size'] < 2*node_min_size    ){
      
          all_trees[[j]]$tree_matrix[terminals[1],'terminal']=1
      
          } 
      else{
    
      # shall we split the first candidate terminal node?
      split = rbinom( 1 , 1 , alpha * (1+depth[1])^(-beta) )
      
      if (split == 1){
      # choose a split variable
      split_variable = sample(1:ncol(X), 1, prob = s)
      
      # which observations are currently in the splitting node
      obs_to_partition =  which( all_trees[[j]]$node_indices == terminals[1])
      
      # find all possible values that we could use for a split
      possible_split_value = sort(X_tree_j [ obs_to_partition , split_variable]) 
      
      # choose a split value
      split_value =  sample(  possible_split_value[ (node_min_size+1): (length(possible_split_value)-node_min_size+1)  ]  ,1 )

     
      
      # Add two extra rows to the tree in question
      all_trees[[j]]$tree_matrix = rbind(all_trees[[j]]$tree_matrix,
                                   c(1, NA, NA, NA, NA, NA, NA, NA), # Make sure they're both terminal
                                   c(1, NA, NA, NA, NA, NA, NA, NA))
      
      # new number of nodes in tree j
      nn_nodes_j = nrow(all_trees[[j]]$tree_matrix)
      
      curr_parent = all_trees[[j]]$tree_matrix[terminals[1], 'parent'] # Make sure to keep the current parent in there. Will be NA if at the root node
      
      all_trees[[j]]$tree_matrix[terminals[1],1:6] = c(0, # Now not temrinal
                                                  nn_nodes_j - 1, # child_left is penultimate row
                                                  nn_nodes_j,  # child_right is penultimate row
                                                  curr_parent,
                                                  split_variable,
                                                  split_value)
      
      #  Fill in the parents of these two nodes
      all_trees[[j]]$tree_matrix[ nn_nodes_j ,'parent'] = terminals[1]
      all_trees[[j]]$tree_matrix[ nn_nodes_j - 1 ,'parent'] = terminals[1]
      
      # Now call the fill function on this tree
      all_trees[[j]] = fill_tree_details(all_trees[[j]], X_tree_j)
      
      ### Add the new children nodes in the set of possible splitting nodes
      terminals=c(terminals, nn_nodes_j -1 , nn_nodes_j )
      ## the depth of the new nodes is the depth of the split variable +1
      depth =c(depth, depth[1]+1, depth[1]+1 )
         
      } else{
        
        # if the draw from the bernoulli is equal to 0, we simply remove the node from the splitting nodes
        # we don't split the node, hence we remove it among the possible splitting terminal nodes
        all_trees[[j]]$tree_matrix[terminals[1],'terminal'] = 1
      }
      
      }
      
      ### Remove the node from the possible splitting nodes
      terminals = terminals[-1]
      depth = depth[-1]
    }
   
    # for each terminal node, draw a value of mu from the prior
    term = which(all_trees[[j]]$tree_matrix[,'terminal']==1)
    all_trees[[j]]$tree_matrix[ term , 'mu' ] =rnorm( length(term), mu_mu , sigma2_mu)
    
  
   } # End of loop through trees
   
   return(all_trees)
   
}








