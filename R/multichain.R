#' @title Build the multichain
#' @description "add_chain" adds a chain to the multichain list. "all_pruned" checks
#'    if all chains have been pruned already.
#'
#' @param chain: List of output values with entries as explained below.
#' @param num_estimate: int.
#'    Number of parameters to estimate.
#' @param estimate_idx: list of int.
#'    Indices of parameters to estimate in the model's full parameter list.
#' @param initial_values: list of float.
#'    Starting values for parameters to estimate, taken from the parameters'
#'    nominal values in the model or explicitly specified in 'options'.
#' @param initial_position: list of float.
#'    Starting position of the MCMC walk in parameter space (log10 of 'initial_values').
#' @param position: list of float.
#'    Current position of MCMC walk in parameter space, i.e. the most
#'    recently accepted move.
#' @param test_position: list of float.
#'    Proposed MCMC mmove.
#' @param acceptance: int.
#'    Number of accepted moves.
#' @param T: float.
#'    Current value of the simulated annealing temperature.
#' @param T_decay: float.
#'    Constant for exponential decay of 'T', automatically calculated such
#'    that T will decay from 'options$T_init' down to 1 over the first
#'    'options$anneal_length' setps.
#' @param sig_value: float.
#'    Current value of sigma, the scaling factor for the proposal distribution.
#'    The MCMC algorithm dynamically tunes this to maintain the aaceptance
#'    rate specified in 'options$accept_rate_target'.
#' @param iter: int.
#'    Current MCMC step number.
#' @param start_iter: int.
#'    Starting MCMC step number.
#' @param ode_options: list.
#'    Options for the ODE integrator, currently just 'rtol' for relative
#'    tolerance and 'atol' for absolute tolerance.
#' @param initial_prior: float.
#'    Starting prior value, i.e. the value at 'initial_position'.
#' @param initial_likelihood: float.
#'    Starting likelihood value, i.e. the value at 'initial_position'.
#' @param initial_posterior: float.
#'    Starting posterior value, i.e. the value at 'initial_position'.
#' @param accept_prior: float.
#'    Current prior value i.e. the value at 'position'.
#' @param accept_likelihood: float.
#'    Current likelihood value i.e. the value at 'position'.
#' @param accept_posterior: float.
#'    Current posterior value i.e. the value at 'position'.
#' @param test_prior: float.
#'    Prior value at 'test_position'.
#' @param test_likelihood: float.
#'    Likelihood value at 'test_position'.
#' @param test_posterior: float.
#'    Posterior value at 'test_position'.
#' @param hessian: array of float.
#'    Current hessian of the posterior landscape. Size is
#'    'num_estimate' x 'num_estimate'.
#' @param positions: array of float.
#'    Trace of all proposed moves. Size is 'num_estimate' x 'nsteps'.
#' @param priors: array of float.
#'    Trace of all priors corresponding to 'positions'. Length is 'nsteps'.
#' @param likelihoods: array of float.
#'    Trace of all likelihoods corresponding to 'positions'. Length is 'nsteps'.
#' @param posteriors: array of float.
#'    Trace of all posteriors corresponding to 'positions'. Length is 'nsteps'.
#' @param alphas: array of float. Trace of 'alpha' parameter and calculated values. Length is 'nsteps'.
#' @param sigmas: array of float. Trace of 'sigma' parameter and calculated values. Length is 'nsteps'.
#' @param delta_posteriors: array of float. Trace of 'delta_posterior' parameter and calculated values. Length is 'nsteps'.
#' @param ts: array of float. Trace of 'T' parameter and calculated values. Length is 'nsteps'.
#' @param accepts: logical array.
#'    Trace of wheter each proposed move was accepted or not.
#'    Length is 'nsteps'.
#' @param rejects: logical array. Trace of wheter each proposed move was rejected or not. Length is 'nsteps'.
#' @param hessians: array of float.
#'    Trace of all hessians. Size is 'num_estimate' x 'num_estimate' x 'num_hessians'
#'    where 'num_hessians' is the actual number of hessians to be calculated.

#' @param options: list with entries as explained below.
#'    Options set -- defines the problem and sets some
#'    parameters to control the MCMC algorithm.
#' @param model: List of model parameters - to estimate.
#'    The parameter objects must each have a
#'    'value' attribute containing the parameter's numerical value.
#' @param estimate_params: list.
#'    List of parameters to estimate, all of which must also be listed
#'    in 'options$model$parameters'.
#' @param initial_values: list of float, optional.
#'    Starting values for parameters to estimate. If omitted, will use
#'    the nominal values from 'options$model$parameters'
#' @param step_fn: callable f(output), optional.
#'    User callback, called on every MCMC iteration.
#' @param likelihood_fn: callable f(output, position).
#'    User likelihood function.
#' @param prior_fn: callable f(output, position), optional.
#'    User prior function. If omitted, a flat prior will be used.
#' @param nsteps: int.
#'    Number of MCMC iterations to perform.
#' @param use_hessian: logical, optional.
#'    Wheter to use the Hessian to guide the walk. Defaults to FALSE.
#' @param rtol: float or list of float, optional.
#'    Relative tolerance for ode solver.
#' @param atol: float or list of float, optional.
#'    Absolute tolerance for ode solver.
#' @param norm_step_size: float, optional.
#'    MCMC step size. Defaults to a reasonable value.
#' @param hessian_period: int, optional.
#'    Number of MCMC steps between Hessian recalculations. Defaults
#'    to a reasonable but fairly large value, as Hessian calculation is expensive.
#' @param hessian_scale: float, optional.
#'    Scaling factor used in generating Hessian-guided steps. Defaults to a
#'    reasonable value.
#' @param sigma_adj_interval: int, optional.
#'    How often to adjust 'output$sig_value' while annealing to meet
#'    'accept_rate_target'. Defaults to a reasonable value.
#' @param anneal_length: int, optional.
#'    Length of initial "burn-in" annealing period. Defaults to 10% of
#'    'nsteps', or if 'use_hessian' is TRUE, to 'hessian_period' (i.e.
#'    anneal until first hessian is calculated)
#' @param T_init: float, optional.
#'    Initial temperature for annealing. Defaults to a resonable value.
#' @param accept_rate_target: float, optional.
#'    Desired acceptance rate during annealing. Defaults to a reasonable value.
#'    See also 'sigma_adj_interval' above.
#' @param sigma_max: float, optional.
#'    Maximum value for 'output$sig_value'. Defaults to a resonable value.
#' @param sigma_min: float, optional.
#'    Minimum value for 'output$sig_value'. Defaults to a resonable value.
#' @param sigma_step: float, optional.
#'    Increment for 'output$sig_value' adjustments. Defaults to a resonable
#'    value. To eliminate adaptive step size, set sigma_step to 1.
#' @param thermo_temp: float in the range [0,1], optional.
#'    Temperature for thermodynamic integration support. Used to scale
#'    likelihood when calculating the posterior value. Defaults to 1,
#'    i.e. no effect.

#' @param multichain: list.
#'    List of chains to be built by "add_chain".
#'
#'
#' @return multichain: list of chains.
#'


add_chain = function(multichain, chain){
  # Add an MCMC chain to the set.
  multichain$chains = append(multichain$chains, list(chain))
  return(multichain)
}

### prune_all_chains does not work due to structure of R lists.
# need to prune each chain individually and then add to the multichain with add_chain
prune_all_chains = function(multichain, options, burn, thin=1){
  # Iterates over all the chains and prunes each one with the specified arguments.
  
  for (chain in multichain$chains){
    chain = prune(chain, options, burn ,thin)
  }
  
  # If any chains are empty after pruning (i.e. there were no accepts)
  # then remove them from the list
  for (chain in multichain$chains){
    if (length(chain$positions)==0){
      cat("WARNING: Chain was empty after pruning and is being removed")
      multichain$chains[which(multichain$chains$positions==chain$positions)]=NULL
    }
  }
  
  return(multichain)
}

all_pruned = function(multichain){
  # Indicates whether all chains have been pruned already.
  
  if (length(multichain$chains)==0){
    stop("There are no chains.")
  }
  for (chain in multichain$chains){
    if(chain$pruned == FALSE){
      return(FALSE)
    } else {
      return(TRUE)
    }
  }
}



pool_chains = function(multichain){
  # Pool the chains into a single set of pooled positions stored along with the MCMCSet.
  if (length(multichain$chains)==0){
    stop("There are no chains.")
  }
  
  # First, count the total number of steps after pruning and make sure
  # all chains have been pruned.
  total_positions = 0
  for (chain in multichain$chains){
    if(chain$pruned==FALSE){
      stop("The chains have not yet been pruned.")
    } else {
      total_positions = total_positions + dim(chain$positions)[1]
    }
  }
  # Allocate enough space for the pooled positions
  multichain$pooled_positions = matrix(nrow=total_positions, ncol=multichain$chains[[1]]$num_estimate, byrow=TRUE)
  
  # Iterate again, filling in the pooled positions
  start_index=0
  for (chain in multichain$chains){
    last_index = start_index + dim(chain$positions)[1]
    multichain$pooled_positions[start_index:last_index,] = chain$positions
    start_index = last_index
  }
  return(multichain)
}

get_sample_position = function(multichain){
  if(length(multichain$chains)==0){
    stop("There are no chains.")
  }
  
  if(is.null(multichain$pooled_positions)){
    stop("Cannot get a sample position until the chains have been pooled.")
  }
  
  if (length(multichain$pooled_positions)==0){
    stop("There are no positions in the combined pool of positions.")
  }
  
  rand_index = sample(1, dim(multichain$pooled_positions)[1])
  return(multichain$pooled_positions[rand_index])
}


initialize_and_pool = function(multichain, chains, options, burn, thin=1){
  # Adds the chains to the multichain and prunes and pools them
  for (chain in chains){
    chain = prune(chain, options, burn, thin=1)
    multichain = add_chain(multichain)
  }
  pool_chains(mulichain)
  return(multichain)
}

maximum_likelihood = function(multichain){
  # Return the maximum log likelihood (minimum negative log likelihood)
  # from the set of chains, along with the position giving the maximum  likelihood.
  if (length(multichain$chains)==0){
    stop("There are no chains")
  }
  
  max_likelihood = Inf
  max_likelihood_position = NULL
  for (chain in multichain$chains){
    if (length(chain$likelihoods)>0){
      chain_max_likelihood_index = which.min(chain$likelihoods)
      chain_max_likelihood = chain$likelihoods[chain_max_likelihood_index]
      if (chain_max_likelihood < max_likelihood){
        max_likelihood = chain_max_likelihood
        max_likelihood_position = chain$positions[chain_max_likelihood_index,]
      }
    }
  }
  # Check if there are no positions
  if (is.null(max_likelihood_position)){
    stop ("The maximum likelihood could not be determined because there are no accepted positions")
  }
  return(list("max_likelihood" = max_likelihood, "max_likelihood_position"=max_likelihood_position))
}

maximum_posterior = function(multichain){
  # Return the maximum log posterior (minimum negative log posterior)
  # from the set of chains, along with the position giving the maximum  posterior.
  if (length(multichain$chains)==0){
    stop("There are no chains")
  }
  
  max_posterior = Inf
    max_posterior_position = NULL
  for (chain in multichain$chains){
    if (length(chain$posteriors)>0){
      chain_max_posterior_index = which.min(chain$posteriors)
      chain_max_posterior = chain$posteriors[chain_max_posterior_index]
      if (chain_max_posterior < max_posterior){
        max_posterior = chain_max_posterior
        max_posterior_position = chain$positions[chain_max_posterior_index,]
      }
    }
  }
  # Check if there are no positions
  if (is.null(max_posterior_position)){
    stop ("The maximum posterior could not be determined because there are no accepted positions")
  }
  return(list("max_posterior" = max_posterior, "max_posterior_position"=max_posterior_position))
}


