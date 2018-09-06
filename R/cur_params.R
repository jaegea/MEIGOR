#' @title Current values of all model parameters
#' 
#' @description For a given set of values for the parameters to be estimated, 
#'    this method returns an array containing the actual (not log-transformed)
#'    values of all model parameters, not just those to be estimated,
#'    in the same order as specified in the model. This is helpful 
#'    when simulating the model at a given position in parameter space.
#' 
#' @param options list with entries as explained below. 
#'    Options set -- defines the problem and sets some
#'    parameters to control the MCMC algorithm.
#' model: List of model parameters - to estimate. 
#'    The parameter objects must each have a
#'    'value' attribute containing the parameter's numerical value.
#' estimate_params: list. 
#'    List of parameters to estimate, all of which must also be listed
#'    in 'options$model$parameters'.
#' initial_values: list of float, optional. 
#'    Starting values for parameters to estimate. If omitted, will use
#'    the nominal values from 'options$model$parameters'
#' step_fn: callable f(output), optional. 
#'    User callback, called on every MCMC iteration.
#' likelihood_fn: callable f(output, position). 
#'    User likelihood function.
#' prior_fn: callable f(output, position), optional. 
#'    User prior function. If omitted, a flat prior will be used.
#' nsteps: int. 
#'    Number of MCMC iterations to perform.
#' use_hessian: logical, optional. 
#'    Wheter to use the Hessian to guide the walk. Defaults to FALSE.
#' rtol: float or list of float, optional. 
#'    Relative tolerance for ode solver.
#' atol: float or list of float, optional. 
#'    Absolute tolerance for ode solver.    
#' norm_step_size: float, optional. 
#'    MCMC step size. Defaults to a reasonable value.
#' hessian_period: int, optional. 
#'    Number of MCMC steps between Hessian recalculations. Defaults
#'    to a reasonable but fairly large value, as Hessian calculation is expensive.
#' hessian_scale: float, optional. 
#'    Scaling factor used in generating Hessian-guided steps. Defaults to a
#'    reasonable value.
#' sigma_adj_interval: int, optional. 
#'    How often to adjust 'output$sig_value' while annealing to meet
#'    'accept_rate_target'. Defaults to a reasonable value.
#' anneal_length: int, optional. 
#'    Length of initial "burn-in" annealing period. Defaults to 10% of
#'    'nsteps', or if 'use_hessian' is TRUE, to 'hessian_period' (i.e. 
#'    anneal until first hessian is calculated)
#' T_init: float, optional. 
#'    Initial temperature for annealing. Defaults to a resonable value.
#' accept_rate_target: float, optional. 
#'    Desired acceptance rate during annealing. Defaults to a reasonable value.
#'    See also 'sigma_adj_interval' above.
#' sigma_max: float, optional. 
#'    Maximum value for 'output$sig_value'. Defaults to a resonable value.
#' sigma_min: float, optional. 
#'    Minimum value for 'output$sig_value'. Defaults to a resonable value. 
#' sigma_step: float, optional. 
#'    Increment for 'output$sig_value' adjustments. Defaults to a resonable 
#'    value. To eliminate adaptive step size, set sigma_step to 1.
#' thermo_temp: float in the range [0,1], optional. 
#'    Temperature for thermodynamic integration support. Used to scale
#'    likelihood when calculating the posterior value. Defaults to 1,
#'    i.e. no effect.
#' 
#' 
#' @param output List of output values with entries as explained below. 
#' num_estimate: int. 
#'    Number of parameters to estimate.
#' estimate_idx: list of int.  
#'    Indices of parameters to estimate in the model's full parameter list.
#' initial_values: list of float.  
#'    Starting values for parameters to estimate, taken from the parameters'
#'    nominal values in the model or explicitly specified in 'options'.
#' initial_position: list of float. 
#'    Starting position of the MCMC walk in parameter space (log10 of 'initial_values').
#' position: list of float. 
#'    Current position of MCMC walk in parameter space, i.e. the most
#'    recently accepted move.
#' test_position: list of float. 
#'    Proposed MCMC mmove.
#' acceptance: int. 
#'    Number of accepted moves.
#' T: float.  
#'    Current value of the simulated annealing temperature.
#' T_decay: float.  
#'    Constant for exponential decay of 'T', automatically calculated such 
#'    that T will decay from 'options$T_init' down to 1 over the first
#'    'options$anneal_length' setps.
#' sig_value: float.  
#'    Current value of sigma, the scaling factor for the proposal distribution.
#'    The MCMC algorithm dynamically tunes this to maintain the aaceptance
#'    rate specified in 'options$accept_rate_target'.
#' iter: int.  
#'    Current MCMC step number.
#' start_iter: int.  
#'    Starting MCMC step number.
#' ode_options: list.  
#'    Options for the ODE integrator, currently just 'rtol' for relative
#'    tolerance and 'atol' for absolute tolerance.
#' initial_prior: float. 
#'    Starting prior value, i.e. the value at 'initial_position'.
#' initial_likelihood: float.  
#'    Starting likelihood value, i.e. the value at 'initial_position'.
#' initial_posterior: float.  
#'    Starting posterior value, i.e. the value at 'initial_position'.
#' accept_prior: float.  
#'    Current prior value i.e. the value at 'position'.
#' accept_likelihood: float.  
#'    Current likelihood value i.e. the value at 'position'.
#' accept_posterior: float.  
#'    Current posterior value i.e. the value at 'position'.
#' test_prior: float.  
#'    Prior value at 'test_position'.
#' test_likelihood: float.  
#'    Likelihood value at 'test_position'.
#' test_posterior: float.  
#'    Posterior value at 'test_position'.
#' hessian: array of float. 
#'    Current hessian of the posterior landscape. Size is 
#'    'num_estimate' x 'num_estimate'.
#' positions: array of float.  
#'    Trace of all proposed moves. Size is 'num_estimate' x 'nsteps'.
#' priors: array of float. 
#'    Trace of all priors corresponding to 'positions'. Length is 'nsteps'.
#' likelihoods: array of float. 
#'    Trace of all likelihoods corresponding to 'positions'. Length is 'nsteps'.
#' posteriors: array of float. 
#'    Trace of all posteriors corresponding to 'positions'. Length is 'nsteps'.
#' alphas: array of float. Trace of 'alpha' parameter and calculated values. Length is 'nsteps'.
#' sigmas: array of float. Trace of 'sigma' parameter and calculated values. Length is 'nsteps'.
#' delta_posteriors: array of float. Trace of 'delta_posterior' parameter and calculated values. Length is 'nsteps'.
#' ts: array of float. Trace of 'T' parameter and calculated values. Length is 'nsteps'.
#' accepts: logical array. 
#'    Trace of wheter each proposed move was accepted or not.
#'    Length is 'nsteps'.
#' rejects: logical array. Trace of wheter each proposed move was rejected or not. Length is 'nsteps'.
#' hessians: array of float. 
#'    Trace of all hessians. Size is 'num_estimate' x 'num_estimate' x 'num_hessians'
#'    where 'num_hessians' is the actual number of hessians to be calculated.

#' 
#' 
#' 
#' 
#' @param position list of float, optional. 
#'    log10 of the values of the parameters being estimated.
#'    If omitted, 'output$position' (the most recent accepted output move)
#'    will be used. The model's nominal values will be used for all parameters
#'    *not* being estimated, regardless.
#'    
#' @return A list of the values of all model parameters.

cur_params = function(output,options,position=NULL){
  if(is.null(position)){
    position = output$position
  }
  values = as.numeric(options$model$parameters$value)
  values[output$estimate_idx] = 10^position
  return(values)
}
