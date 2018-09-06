#' @title Convergence criterion
#' @description Apply the Gelman-Rubin convergence criterion. If the chains have converged than the value returned should be close to 1 (conventionally less than 1.1). This is a one-way relation such that getting a value close to 1 does not guarantee convergence, but having a value much larger than 1 indicates lack of convergence.
#'
#' @param multichain: list.
#'      List of chains built with 'add_chain'. Chains should be pruned before added
#'      to the multichain.
#'
#'
#' @return The value of the Gelman-Rubin test for each parameter.
#'
#'


convergence_criterion = function(multichain){
    # Make sure there is more than one chain:
    if (length(multichain$chains) <= 1){
        return(NULL)
    }
    
    # Make sure the chains have been pruned
    if (all_pruned(multichain)==FALSE){
        stop("The chains should be pruned before calculating convergence")
    }
    
    # Iterate over the MCMC set to find the minimum number of steps in the set
    chain_set = list()
    min_accepts = Inf
    for (chain in multichain$chains){
        chain_set = append(chain_set, list(chain$positions))
        if (dim(chain$positions)[1] < min_accepts){
            min_accepts = dim(chain$positions)[1]
        }
    }
    
    # Make sure we have some steps
    if (min_accepts ==0){
        return(NULL)
    }
    
    # Truncate the chains to make them all the length of the one with the fewest accepts
    for (i in 1:length(chain_set)){
        chain = chain_set[[i]]
        chain_set[[i]] = chain[(dim(chain)[1]-min_accepts+1):dim(chain)[1],]
    }
    
    # Run the calculations on the chain set
    W = within_chain_variances(chain_set)
    var = parameter_variance_estimates(chain_set)
    
    return(sqrt(var/W))
}



check_chain_lengths = function(chain_set){
  # Checks to make sure there is more than one chain in the set, 
  # and that all chains are the same length.
  if (length(chain_set)<2){
    stop("To calculate the convergence criterion, there must be more than one chain")
  }
  n = length(chain_set[[1]])
  for (i in 1:length(chain_set)){
    if (length(chain_set[[i]]) != n){
      stop("To calculate the within-chain variances, all chains must be the same length")
    }
  }
}

within_chain_variances = function(chain_set){
  # Calculate the vector of average within-chain variances, W.
  # Takes a list of chains, each expressed as an array of positions.
  # Note that all chains should be the same length.
  
  # Make sure all chains are the same length
  check_chain_lengths(chain_set)
  
  #Calculate each within-chain variance (s_j^2):
  chain_variances = list()
  for (chain in chain_set){
    chain_variances = append(chain_variances, list(diag(var(chain))))
  }
  
  # Calculate the average within-chain variance (W):
  W = colMeans(matrix(unlist(chain_variances), nrow=length(chain_set), byrow=TRUE))
  
  return(W)
}

between_chain_variances = function(chain_set){
  # Caluclate the vector of between-chain variances, B.
  # Takes a list of chains (each expressed as an array of positions).
  # Note that all chains should be the same length.
  
  # Make sure all chains are the same length:
  check_chain_lengths(chain_set)
  
  # Calculate each within-chain average (i.e. psi_{.j})
  chain_averages = list()
  if (is.null(dim(chain_set[[1]]))){
    for (chain in chain_set){
      chain_averages = append(chain_averages, list(mean(chain)))
    }
    n = length(chain_set[[1]])
  } else {
    for (chain in chain_set){
      chain_averages = append(chain_averages, list(colMeans(chain)))
    }
    n = length(chain_set)
  }
  
  # Caluclate the between-chain variance (B):
  # n = length(chain_set)
  B = n * diag(var(matrix(unlist(chain_averages), nrow = length(chain_set), byrow=TRUE)))
  
  return(B)
}

parameter_variance_estimates = function(chain_set){
  # Calculate the best estimate of the variances of the parameters given
  # the chains that have been run, that is,
  
  #  \widehat{\mathrm{var}}^+ (\psi|y) = \frac{n-1}{n} W + \frac{1}{n} B
  
  # Takes a list of chains (each expressed as an array of positions).
  # Note that all chains should be the same length.
  
  if (is.null(dim(chain_set[[1]]))){
    n = as.numeric(length(chain_set[[1]]))
  } else {
    n =  as.numeric(dim(chain_set[[1]])[1])
  }
  
  W = within_chain_variances(chain_set)
  B = between_chain_variances(chain_set)
  
  return ((((n-1)/n) * W) + ((1/n) * B))
}


# -----TESTS-------
# The following are tests of the convergence criterion calculations using
# the simple data given below.

test_data = list(c(1,2,3), c(4,5,6))

test_within_chain_variances = function(){
  # The average of the first chain is 2; the second is 5
  # The variance of the first chain is 1/2 * 2 = 1
  # The variance of the second chain is also 1/2 * 2 = 1
  # The average of the two variances is therefore 1.0
  if (within_chain_variances(test_data) != 1){
   stop("Failed to correctly calculate within-chain variance!") 
  } else{
    return(TRUE)
  }
}

test_between_chain_variances = function(){
  # The average of the first chain is 2; the second is 5
  # The average of the two averages is thus 3.5
  # The variance between the chains is thus
  # = 3/1 * [(3.5-2)^2 + (3.5-5)^2]
  # = 3 * [9/4 + 9/4]
  # = 54/4 = 13.5
  if (between_chain_variances(test_data) != 13.5){
    stop("Failed to correctly calculate between-chain variance!") 
  } else {
    return(TRUE)
  }
}

test_parameter_variance_estimates = function(){
  # For the test data, W is 1.0 and B is 13.5 (see tests above)
  # The variance estimate is thus
  # = (2/3 * 1.0) + (1/3 * 13.5)
  # = 2/3 + 18/4
  # = 8/12 + 54/12
  # = 62/12
  if (parameter_variance_estimates(test_data) != 62/12){
    stop("Failed to correctly calculate parameter variance estimates!")
  } else {
    return(TRUE)
  }
}

test_convergence_criterion = function(){
  # For the test data, the estimated variance is 62/12, and the
  # within-chain variance is 1.0 (see other tests above)
  # The convergence criterion is therefore
  # sqrt( (62/12) / 1.0)
  if (sqrt(62/12.) != convergence_criterion(test_data)){
    stop("Failed to correctly calculate the convergence criterion!")
  } else{
    return(TRUE)
  }
}








