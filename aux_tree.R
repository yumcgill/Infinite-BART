# -------------------------------------------------------------------------#
# Description: this script contains auxiliar functions needed to update    #
# the trees with details and to map the predicted values to each obs       #
# -------------------------------------------------------------------------#

# 1. fill_tree_details: takes a tree matrix and returns the number of obs in each node in it and the indices of each observation in each terminal node
# 2.a get_predictions: gets the predicted values from a current set of trees for plain BART
# 2.b get_predictions_ibp: gets the predicted values from a current set of trees for BART-IBP
# 3. get_children: it's a function that takes a node and, if the node is terminal, returns the node. If not, returns the children and calls the function again on the children
# 4. resample: an auxiliar function
# 5. get_ancestors: get the ancestors of all terminal nodes in a tree
# 6. update_s: full conditional of the vector of splitting probability.
# 7. get_number_distinct_cov: given a tree, it returns the number of distinct covariates used to create its structure



# Fill_tree_details -------------------------------------------------------
# OPTIMISED: instead of scanning all n observations once per node
# (O(n * #nodes), with repeated `node_indices == parent` comparisons), the
# observations are pushed down the tree once: each internal node splits the
# index set of its parent. Cost is O(n * depth). Output is identical.
# (Same assumption as the original code: a parent row precedes its children.)

fill_tree_details = function(curr_tree, X) {

  tree_matrix = curr_tree$tree_matrix
  n  = nrow(X)
  nr = nrow(tree_matrix)

  node_indices = rep.int(1L, n)
  node_size    = numeric(nr)
  node_size[1] = n

  if (nr > 1) {
    terminal = tree_matrix[, 'terminal']
    cl  = tree_matrix[, 'child_left']
    cr  = tree_matrix[, 'child_right']
    sv  = tree_matrix[, 'split_variable']
    sval = tree_matrix[, 'split_value']

    members = vector('list', nr)
    members[[1]] = seq_len(n)
    for (i in which(terminal == 0)) {          # internal nodes, in row order
      idx = members[[i]]
      go_left = X[idx, sv[i]] < sval[i]
      members[[cl[i]]] = idx[go_left]
      members[[cr[i]]] = idx[!go_left]
    }
    node_size = as.numeric(lengths(members))
    for (i in which(terminal == 1)) {
      if (length(members[[i]])) node_indices[members[[i]]] = i
    }
  }

  tree_matrix[, 'node_size'] = node_size
  return(list(tree_matrix = tree_matrix,
              node_indices = node_indices))

} # End of function


# NEW: prediction of ONE tree at every row of X (no masking by W) ----------
# This is the basic building block used everywhere below: it routes every
# row of X to a leaf and returns the leaf value.
predict_tree_full = function(tree, X) {
  tm = tree$tree_matrix
  if (nrow(tm) == 1) return(rep.int(tm[1, 'mu'], nrow(X)))
  as.numeric(tm[fill_tree_details(tree, X)$node_indices, 'mu'])
}


# Get predictions ---------------------------------------------------------
# OPTIMISED: loop instead of recursion (the recursion copied the list of
# trees at every level: O(K^2) copies).

get_predictions = function(trees, X, single_tree = FALSE) {

  # Stop nesting problems in case of multiple trees
  if(is.null(names(trees)) & (length(trees) == 1)) trees = trees[[1]]
  if (!is.null(names(trees))) return(predict_tree_full(trees, X))

  predictions = numeric(nrow(X))
  for (k in seq_along(trees)) predictions = predictions + predict_tree_full(trees[[k]], X)
  return(predictions)
}

#### for BART-IBP
# W: vector (single tree) or n x K matrix (one column per tree)
get_predictions_ibp = function(trees, X, W, single_tree = FALSE) {

  if(is.null(names(trees)) & (length(trees) == 1)) trees = trees[[1]]

  if (!is.null(names(trees))) {           # a single tree
    w = if (is.null(dim(W))) W else W[, 1]
    return(predict_tree_full(trees, X) * w)
  }

  W = as.matrix(W)
  predictions = numeric(nrow(X))
  for (k in seq_along(trees)) {
    predictions = predictions + predict_tree_full(trees[[k]], X) * W[, k]
  }
  return(predictions)
}




# get_children ------------------------------------------------------------

get_children = function(tree_mat, parent) {
  # Create a holder for the children
  all_children = NULL
  if(as.numeric(tree_mat[parent,'terminal']) == 1) {
    # If the node is terminal return the list so far
    return(c(all_children, parent))
  } else {
    # If not get the current children
    curr_child_left = as.numeric(tree_mat[parent, 'child_left'])
    curr_child_right = as.numeric(tree_mat[parent, 'child_right'])
    # Return the children and also the children of the children recursively
    return(c(all_children,
             get_children(tree_mat,curr_child_left),
             get_children(tree_mat,curr_child_right)))
  }
}






# Sample function ----------------------------------------------------------

resample <- function(x, ...) x[sample.int(length(x), size=1), ...]

update_s = function(var_count, p, alpha_s){
  s_ = rdirichlet(1, alpha_s/p + var_count)
  return(s_)
}

get_number_distinct_cov <- function(tree){
  
  # Select the rows that correspond to internal nodes
  which_terminal = which(tree$tree_matrix[,'terminal'] == 0)
  # Get the covariates that are used to define the splitting rules
  num_distinct_cov = length(unique(tree$tree_matrix[which_terminal,'split_variable']))
  
  return(num_distinct_cov)
}

sample_move = function(curr_tree, i, nburn){
  
  if (nrow(curr_tree$tree_matrix) == 1 || i < max(floor(0.1*nburn), 10)) {
    type = 'grow'
  } else {
    type = sample(c('grow', 'prune', 'change'), 1)
  }
  return(type)
}