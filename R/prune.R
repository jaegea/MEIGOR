#' @title Truncates the chain to the thinned, mixed, accepted positions.
#' 
#' @description Note that if the chain has already been pruned then 
#' this function does nothing.
#' 
#' Side Effects: 
#' 
#' - After this method is called, output$positions is (destructively) set 
#' to the mixed, accepted positions. 
#' - The priors, likelihoods, and posteriors, sigmas, alphas, etc. 
#' are modified to only record values for the mixed, accepted positions. 
#' - The ``output$pruned`` field is set to TRUE, indicating that the walk 
#' has been irreversibly pruned. 
#' - The indices of the thinned, mixed, accepted steps are recorded in 
#' the field ``output$thinned_accept_steps``. This can be useful for 
#' plotting walks from multiple chains. 
#' - The values of the ``burn`` and ``thin`` parameters are recorded 
#' in the options object as ``output$options$burn`` and ``output$options$thin``.
#' 
#' 
#' @param options: list with entries as explained below. 
#'    Options set -- defines the problem and sets some
#'    parameters to control the MCMC algorithm.
#' @param model: List of model parameters - to estimate.. 
#'    The parameter objects must each have a
#'    'value' attribute containing the parameter's numerical value.
#' @param estimate_params: list. 
#'    List of parameters to estimate, all of which must also be listed
#'    in 'options$model$parameters'.
#' @param initail_values: list of float, optional. 
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
#'    Wheter to use the Hessian to giude the walk. Defaults to false.
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
#'    'nsteps', or if 'use_hessian' is true, to 'hessian_period' (i.e. 
#'    anneal until first hessian is calculated)
#' @param T_init: float, optional. 
#'    Initial temperature for annealing. Defaults to a resonable values.
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
#' 
#' 
#' 
#' @param output: List of output values with entries as explained below.. 
#' @param num_estimate: int. 
#'    Number of parameters to estimate.
#' @param estimate_idx: list of int. 
#'    Indices of parameters to estimate in the model's full parameter list.
#' @param initial_values : list of float. 
#'    Starting values for parameters to estimate, taken from the parameters'
#'    nominal values in the model or explicitly specified in `options`.
#' @param initial_position : list of float. 
#'    Starting position of the MCMC walk in parameter space
#'    (log10 of 'initial_values')
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
#'    Constant for exponential decay 'T', automatically calculated such 
#'    that T will decay from 'options$T_init' down to 1 over the first
#'    'options$anneal_length' setps.
#' @param sig_value: float. 
#'    Current value of sigma, the scaling factor for the proposal distribution.
#'    Te MCMC algorithm dynamically tunes this to maintain the aaceptance
#'    rate specified in 'options$accept_rate_target'.
#' @param iter: int. 
#'    Current MCMC step number.
#' @param start_iter: int. 
#'    Starting MCMC step number.
#' @param ode_options: list (dict in python). 
#'    Options for the ODE integrator, currently just 'rtol' for relative
#'    tolerance and 'atol' for absolute tolerance.
#' @param initial_prior: float. 
#'    Starting prior value, i.e. the value at 'initial_position'.
#' @param initial_likelihood: float. 
#'    Starting likelihood value, i.e. the values at 'initial_position'.
#' @param initial_posterior: float. 
#'    Starting posterior value, i.e. the values at 'initial_position'.
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
#'    Trace of all priors corresponding to positions'. Length is 'nsteps'.
#' @param likelihoods: array of float. 
#'    Trace of all likelihoods corresponding to positions'. Length is 'nsteps'.
#' @param posteriors: array of float. 
#'    Trace of all posteriors corresponding to positions'. Length is 'nsteps'.
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
#' 
#' 
#' 
#' @param burn: int.
#'    An integer specifying the number of steps to cut off
#'    from the beginning of the walk.
#' @param thin: int.
#'    An integer specifying how to thin the accepted steps of the walk.
#'    If 1, it returns every step; if 2, every other step; if 5, every fifth step, etc.
#' 
#'    
#' @return The updated output list, aftr the chain is pruned.
 
prune = function(output, options, burn, thin=1){
  thinned_accepts = get_mixed_accepts(output, options, burn, thin)$thinned_accepts
  thinned_accept_steps = get_mixed_accepts(output, options, burn, thin)$thinned_accept_steps
  
  output$positions = thinned_accepts
  output$options$burn = burn
  output$options$thin = thin
  output$pruned = TRUE
  output$thinned_accept_steps = thinned_accept_steps
  
  output$delta_posteriors = output$delta_posteriors[thinned_accept_steps]
  output$ts = output$ts[thinned_accept_steps]
  output$priors = output$priors[thinned_accept_steps]
  output$likelihoods = output$likelihoods[thinned_accept_steps]
  output$posteriors = output$posteriors[thinned_accept_steps]
  
  output$alphas = output$alphas[thinned_accept_steps]
  output$sigmas = output$sigmas[thinned_accept_steps]
  output$accepts = output$accepts[thinned_accept_steps]
  output$rejects = output$rejects[thinned_accept_steps]
  
  #output$delta_posteriors = output$delta_posteriors[-length(output$delta_posteriors)]
  #output$ts = output$ts[-length(output$ts)]
  #output$priors = output$priors[-length(output$priors)]
  #output$likelihoods = output$likelihoods[-length(output$likelihoods)]
  #output$posteriors = output$posteriors[-length(output$posteriors)]
  
  #output$alphas = output$alphas[-length(output$alphas)]
  #output$sigmas = output$sigmas[-length(output$sigmas)]
  #output$accepts = output$accepts[-length(output$accepts)]
  #output$rejects = output$rejects[-length(output$rejects)]

  return(output)
}


 


